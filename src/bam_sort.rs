//! Pure-Rust external (disk-spilling) k-way merge sort for BAM files, keyed on
//! a two-byte aux tag (e.g. the cell-barcode tag `CB`).
//!
//! This replaces the runtime dependency on the external `samtools sort -t CB`
//! binary. It is net-new infrastructure with **no** counterpart in velocyto.py
//! (the Python pipeline shelled out to `samtools`), so the project's
//! one-Rust-function-per-Python-function / no-helpers translation rule does not
//! apply here — the algorithm genuinely needs a spill phase, a merge phase, and
//! a few internal helpers. The whole module is gated behind the `bam` feature.
//!
//! Correctness contract: the downstream `ExInCounter::count` re-sorts the reads
//! within each cell batch itself, so only **contiguity by the full tag value**
//! must be preserved — not byte-for-byte within-group order. That lets us use an
//! unstable in-chunk sort and an arbitrary-but-deterministic position for
//! records that lack the tag (they trail all tagged records).

use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;
use std::fs::File;
use std::num::NonZero;
use std::path::Path;

use anyhow::Context;
use rayon::slice::ParallelSliceMut;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::value::Array;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::RecordBuf;

use crate::constants::BAM_COMPRESSION;

/// Sort key for one record. Variant order is load-bearing: `derive(Ord)` orders
/// `Present` before `Missing`, and `Present(a)` vs `Present(b)` compares the raw
/// bytes lexicographically — exactly the grouping `samtools sort -t TAG` produces,
/// with untagged records gathered into a single trailing group.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum SortKey {
    Present(Vec<u8>),
    Missing,
}

/// Extract the sort key from a record. Keys on the raw tag bytes (no UTF-8
/// assumption); a `Character`-typed tag is treated as a one-byte key.
fn extract_key(record: &RecordBuf, tag: &Tag) -> SortKey {
    match record.data().get(tag) {
        Some(Value::String(s)) => {
            let bytes: &[u8] = s.as_ref();
            SortKey::Present(bytes.to_vec())
        }
        Some(Value::Character(c)) => SortKey::Present(vec![*c]),
        _ => SortKey::Missing,
    }
}

/// Conservative payload estimate for a record and its cached sort key.
fn record_footprint(record: &RecordBuf, tag: &Tag) -> usize {
    let name = record.name().map_or(0, |n| n.len());
    let seq = record.sequence().len();
    let cigar = record.cigar().as_ref().len().saturating_mul(16);
    let aux = record.data().iter().fold(0usize, |n, (_, value)| {
        let payload = match value {
            Value::String(s) | Value::Hex(s) => s.len(),
            Value::Array(a) => match a {
                Array::Int8(v) => v.len(),
                Array::UInt8(v) => v.len(),
                Array::Int16(v) => v.len().saturating_mul(2),
                Array::UInt16(v) => v.len().saturating_mul(2),
                Array::Int32(v) => v.len().saturating_mul(4),
                Array::UInt32(v) => v.len().saturating_mul(4),
                Array::Float(v) => v.len().saturating_mul(4),
            },
            _ => 8,
        };
        n.saturating_add(payload).saturating_add(32)
    });
    let key = match record.data().get(tag) {
        Some(Value::String(s)) => s.len(),
        Some(Value::Character(_)) => 1,
        _ => 0,
    };
    512usize
        .saturating_add(name)
        .saturating_add(seq.saturating_mul(2))
        .saturating_add(cigar)
        .saturating_add(aux)
        .saturating_add(key)
}

/// Build a BAM writer whose BGZF stream uses compression level `BAM_COMPRESSION`
/// (7), matching the old `samtools sort -l 7`, with `workers` threads compressing
/// blocks in parallel (mirrors samtools `-@`). The noodles BAM writer builder
/// cannot set the level directly, so we wrap a pre-configured multithreaded BGZF
/// writer; `bam::io::Writer` accepts any `Write`, and `MultithreadedWriter` is one.
/// Finalize with [`finish_bam`] (not `try_finish`, which only exists on the
/// single-threaded BGZF writer).
fn level7_bam_writer<W: std::io::Write + Send + 'static>(
    w: W,
    workers: usize,
) -> bam::io::Writer<bgzf::io::MultithreadedWriter<W>> {
    let level = bgzf::io::writer::CompressionLevel::new(BAM_COMPRESSION as u8)
        .expect("BAM_COMPRESSION must be a valid BGZF compression level");
    let workers = NonZero::new(workers.max(1)).expect("workers floored at 1");
    let inner = bgzf::io::multithreaded_writer::Builder::default()
        .set_compression_level(level)
        .set_worker_count(workers)
        .build_from_writer(w);
    bam::io::Writer::from(inner)
}

/// Flush and finalize a multithreaded BAM writer: shuts down the compression
/// workers and appends the BGZF EOF block, surfacing any deferred write error.
fn finish_bam<W: std::io::Write + Send + 'static>(
    writer: bam::io::Writer<bgzf::io::MultithreadedWriter<W>>,
) -> std::io::Result<()> {
    writer.into_inner().finish().map(|_| ())
}

/// A heap entry for the k-way merge. Ordered by `(key, source)` only — the
/// `record` payload is ignored for ordering (and `RecordBuf` is not `Ord`). The
/// `source` tiebreaker keeps the merge total-ordered and deterministic.
struct HeapEntry {
    key: SortKey,
    source: usize,
    record: RecordBuf,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.source == other.source
    }
}
impl Eq for HeapEntry {}
impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key
            .cmp(&other.key)
            .then(self.source.cmp(&other.source))
    }
}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Sort `input` BAM by the two-byte `tag` (e.g. "CB"), writing the result to
/// `output` so that all records sharing a tag value are contiguous.
///
/// Uses an external (disk-spilling) k-way merge so inputs larger than RAM sort
/// correctly. `mem_budget_mb` bounds the resident memory of each in-memory chunk
/// (mirrors samtools `-m`); `threads` controls the rayon sort parallelism. Temp
/// spill files are removed when this call returns; the pipeline wrapper keeps
/// interruption state in a separate validated work directory.
pub fn sort_bam_by_tag(
    input: &Path,
    output: &Path,
    tag: &str,
    mem_budget_mb: usize,
    threads: usize,
) -> anyhow::Result<()> {
    let out_parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let tempdir = tempfile::Builder::new()
        .prefix("cellsort_tmp_")
        .tempdir_in(out_parent)
        .with_context(|| format!("creating spill temp dir in {}", out_parent.display()))?;
    sort_bam_by_tag_with_spill_dir(input, output, tag, mem_budget_mb, threads, tempdir.path())
}

fn sort_bam_by_tag_with_spill_dir(
    input: &Path,
    output: &Path,
    tag: &str,
    mem_budget_mb: usize,
    threads: usize,
    spill_dir: &Path,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        tag.len() == 2,
        "BAM sort tag must be exactly two bytes, got {tag:?}"
    );
    let tag_bytes = tag.as_bytes();
    let tag = Tag::from([tag_bytes[0], tag_bytes[1]]);

    // Always admit at least one record per chunk before checking the budget, so
    // a budget smaller than a single record still makes progress (many tiny
    // chunks) rather than spinning. Floor at 1 byte for the same reason.
    let budget_bytes = mem_budget_mb.saturating_mul(1024 * 1024).max(1);

    let mut reader = File::open(input)
        .map(bam::io::Reader::new)
        .with_context(|| format!("opening BAM {}", input.display()))?;
    let header = reader
        .read_header()
        .with_context(|| format!("reading BAM header {}", input.display()))?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build()
        .context("building rayon thread pool for BAM sort")?;

    // ── Spill phase ───────────────────────────────────────────────────────────
    let mut record = RecordBuf::default();
    let mut chunk: Vec<RecordBuf> = Vec::new();
    let mut spilled: Vec<std::path::PathBuf> = Vec::new();
    let mut eof = false;

    loop {
        chunk.clear();
        let mut chunk_bytes: usize = 0;
        loop {
            let n = reader
                .read_record_buf(&header, &mut record)
                .with_context(|| format!("reading record from {}", input.display()))?;
            if n == 0 {
                eof = true;
                break;
            }
            chunk_bytes = chunk_bytes.saturating_add(record_footprint(&record, &tag));
            chunk.push(std::mem::take(&mut record));
            if chunk_bytes >= budget_bytes {
                break;
            }
        }

        // Fast path: the entire input fit one in-memory chunk. Sort and stream
        // straight to the output — no temp files, no merge.
        if spilled.is_empty() && eof {
            pool.install(|| chunk.par_sort_by_cached_key(|r| extract_key(r, &tag)));
            let mut writer = level7_bam_writer(
                File::create(output).with_context(|| format!("creating {}", output.display()))?,
                threads,
            );
            writer
                .write_header(&header)
                .with_context(|| "writing BAM header")?;
            for r in &chunk {
                writer
                    .write_alignment_record(&header, r)
                    .with_context(|| "writing BAM record")?;
            }
            finish_bam(writer).with_context(|| "finalizing BAM")?;
            return Ok(());
        }

        // Otherwise spill this chunk (skip an empty trailing chunk at EOF).
        if !chunk.is_empty() {
            pool.install(|| chunk.par_sort_by_cached_key(|r| extract_key(r, &tag)));
            std::fs::create_dir_all(spill_dir)
                .with_context(|| format!("creating spill directory {}", spill_dir.display()))?;
            let path = spill_dir.join(format!("chunk_{}.bam", spilled.len()));
            let mut writer = level7_bam_writer(
                File::create(&path)
                    .with_context(|| format!("creating spill chunk {}", path.display()))?,
                threads,
            );
            writer
                .write_header(&header)
                .with_context(|| "writing spill chunk header")?;
            for r in &chunk {
                writer
                    .write_alignment_record(&header, r)
                    .with_context(|| "writing spill chunk record")?;
            }
            finish_bam(writer).with_context(|| "finalizing spill chunk")?;
            File::open(&path)?.sync_all()?;
            spilled.push(path);
        }

        if eof {
            break;
        }
    }

    drop(chunk);
    drop(reader);
    let mut pass = 0;
    const FAN_IN: usize = 32;
    while spilled.len() > FAN_IN {
        let mut next = Vec::new();
        for (i, runs) in spilled.chunks(FAN_IN).enumerate() {
            let path = spill_dir.join(format!("merge_{pass}_{i}.bam"));
            merge_runs(runs, &path, &header, &tag, threads)?;
            File::open(&path)?.sync_all()?;
            next.push(path);
        }
        spilled = next;
        pass += 1;
    }
    merge_runs(&spilled, output, &header, &tag, threads)?;
    Ok(())
}

fn merge_runs(
    spilled: &[std::path::PathBuf],
    output: &Path,
    header: &noodles_sam::Header,
    tag: &Tag,
    threads: usize,
) -> anyhow::Result<()> {
    // ── Merge phase (k-way, min-heap) ──────────────────────────────────────────
    let mut readers: Vec<bam::io::Reader<bgzf::io::Reader<File>>> =
        Vec::with_capacity(spilled.len());
    for path in spilled {
        let mut r = File::open(path)
            .map(bam::io::Reader::new)
            .with_context(|| format!("reopening spill chunk {}", path.display()))?;
        // Each chunk carries a copy of the source header; skip past it and reuse
        // the original `header` for decoding records.
        r.read_header()
            .with_context(|| format!("reading spill chunk header {}", path.display()))?;
        readers.push(r);
    }

    let mut heap: BinaryHeap<Reverse<HeapEntry>> = BinaryHeap::new();
    let mut rec = RecordBuf::default();
    for (i, r) in readers.iter_mut().enumerate() {
        let n = r
            .read_record_buf(&header, &mut rec)
            .with_context(|| "priming merge heap")?;
        if n != 0 {
            let key = extract_key(&rec, &tag);
            heap.push(Reverse(HeapEntry {
                key,
                source: i,
                record: std::mem::take(&mut rec),
            }));
        }
    }

    let mut writer = level7_bam_writer(
        File::create(output).with_context(|| format!("creating {}", output.display()))?,
        threads,
    );
    writer
        .write_header(&header)
        .with_context(|| "writing merged BAM header")?;
    while let Some(Reverse(entry)) = heap.pop() {
        writer
            .write_alignment_record(&header, &entry.record)
            .with_context(|| "writing merged BAM record")?;
        let i = entry.source;
        let n = readers[i]
            .read_record_buf(&header, &mut rec)
            .with_context(|| "advancing merge reader")?;
        if n != 0 {
            let key = extract_key(&rec, &tag);
            heap.push(Reverse(HeapEntry {
                key,
                source: i,
                record: std::mem::take(&mut rec),
            }));
        }
    }
    finish_bam(writer).with_context(|| "finalizing merged BAM")?;

    Ok(())
}

/// Reuse a completed sort or restart safely from the source BAM.
pub fn prepare_sorted_bam(
    input: &Path,
    output: &Path,
    tag: &str,
    memory_mb: usize,
    threads: usize,
) -> anyhow::Result<()> {
    use crate::count_checkpoint::InputIdentity;

    #[derive(serde::Serialize, serde::Deserialize, PartialEq, Eq)]
    struct Receipt {
        version: u32,
        input: InputIdentity,
        output: InputIdentity,
        tag: String,
    }

    #[derive(serde::Serialize, serde::Deserialize)]
    struct SortState {
        version: u32,
        input: InputIdentity,
        tag: String,
        output_bytes: Option<u64>,
        output_modified_ns: Option<u128>,
    }

    let receipt_path = output.with_extension("bam.sort.json");
    let lock = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(output.with_extension("bam.sort.lock"))?;
    lock.try_lock()
        .context("another process is sorting this BAM")?;
    let input_identity = InputIdentity::read(input)?;
    let work_dir = output.with_extension("bam.sort-work");
    let state_path = work_dir.join("sort-state.json");
    let state_partial = work_dir.join("sort-state.partial");
    let staged = work_dir.join("sorted.bam");
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));

    if output.exists() {
        if receipt_path.exists() {
            let receipt: Receipt = read_json_bounded(&receipt_path)?;
            anyhow::ensure!(
                receipt
                    == Receipt {
                        version: 1,
                        input: input_identity.clone(),
                        output: InputIdentity::read(output)?,
                        tag: tag.into()
                    },
                "sorted BAM provenance differs; retained unchanged"
            );
            cleanup_owned_sort_work(&work_dir, &input_identity, tag);
            return Ok(());
        }

        let state: SortState = read_json_bounded(&state_path).with_context(|| format!(
            "Unverified sorted BAM {}; retained unchanged. Move it aside before sorting, or resume from an existing counting checkpoint using write-output",
            output.display()
        ))?;
        anyhow::ensure!(
            state.version == 1
                && state.input == input_identity
                && state.tag == tag
                && state.output_bytes == Some(std::fs::metadata(output)?.len())
                && state.output_modified_ns == Some(InputIdentity::read(output)?.modified_ns),
            "sorted BAM provenance differs; retained unchanged"
        );
        write_sort_receipt(&receipt_path, &input_identity, output, tag)?;
        cleanup_owned_sort_work(&work_dir, &input_identity, tag);
        return Ok(());
    }

    if work_dir.exists() {
        if !std::fs::symlink_metadata(&work_dir)?.file_type().is_dir() {
            anyhow::bail!(
                "sort work path exists and is not a directory: {}",
                work_dir.display()
            );
        }
        if state_path.exists() {
            let state: SortState = read_json_bounded(&state_path)?;
            anyhow::ensure!(
                state.version == 1 && state.input == input_identity && state.tag == tag,
                "sort work directory belongs to different inputs; retained unchanged"
            );
            if let (Some(bytes), Some(modified_ns)) = (state.output_bytes, state.output_modified_ns)
            {
                anyhow::ensure!(
                    matches_file_stamp(&staged, bytes, modified_ns)?,
                    "completed staged BAM is missing or changed; retained unchanged"
                );
                std::fs::rename(&staged, output)?;
                File::open(parent)?.sync_all()?;
                write_sort_receipt(&receipt_path, &input_identity, output, tag)?;
                cleanup_owned_sort_work(&work_dir, &input_identity, tag);
                return Ok(());
            }
        } else if std::fs::read_dir(&work_dir)?.next().is_some() {
            anyhow::bail!(
                "unmarked sort work directory exists; retained unchanged: {}",
                work_dir.display()
            );
        }
        std::fs::remove_dir_all(&work_dir)?;
    }
    if receipt_path.exists() {
        std::fs::remove_file(&receipt_path)?;
    }

    std::fs::create_dir(&work_dir)?;
    let mut state = SortState {
        version: 1,
        input: input_identity.clone(),
        tag: tag.into(),
        output_bytes: None,
        output_modified_ns: None,
    };
    write_sort_state(&state_path, &state_partial, &state)?;
    if let Err(error) = sort_bam_by_tag_with_spill_dir(
        input,
        &staged,
        tag,
        memory_mb,
        threads,
        &work_dir.join("spills"),
    ) {
        let _ = std::fs::remove_dir_all(&work_dir);
        return Err(error).context("sorting BAM by barcode");
    }
    File::open(&staged)?.sync_all()?;
    let staged_identity = InputIdentity::read(&staged)?;
    state.output_bytes = Some(staged_identity.bytes);
    state.output_modified_ns = Some(staged_identity.modified_ns);
    write_sort_state(&state_path, &state_partial, &state)?;
    std::fs::rename(&staged, output)?;
    File::open(parent)?.sync_all()?;
    write_sort_receipt(&receipt_path, &input_identity, output, tag)?;
    cleanup_owned_sort_work(&work_dir, &input_identity, tag);
    Ok(())
}

fn read_json_bounded<T: serde::de::DeserializeOwned>(path: &Path) -> anyhow::Result<T> {
    let file = File::open(path)?;
    anyhow::ensure!(file.metadata()?.len() <= 1_048_576, "sort state too large");
    Ok(serde_json::from_reader(file)?)
}

fn matches_file_stamp(path: &Path, bytes: u64, modified_ns: u128) -> anyhow::Result<bool> {
    match crate::count_checkpoint::InputIdentity::read(path) {
        Ok(identity) => Ok(identity.bytes == bytes && identity.modified_ns == modified_ns),
        Err(error)
            if error
                .downcast_ref::<std::io::Error>()
                .is_some_and(|e| e.kind() == std::io::ErrorKind::NotFound) =>
        {
            Ok(false)
        }
        Err(error) => Err(error),
    }
}

fn write_sort_state(
    path: &Path,
    partial: &Path,
    state: &impl serde::Serialize,
) -> anyhow::Result<()> {
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(partial)?;
    serde_json::to_writer(&mut file, state)?;
    file.sync_all()?;
    std::fs::rename(partial, path)?;
    File::open(path.parent().unwrap_or(Path::new(".")))?.sync_all()?;
    Ok(())
}

fn write_sort_receipt(
    path: &Path,
    input: &crate::count_checkpoint::InputIdentity,
    output: &Path,
    tag: &str,
) -> anyhow::Result<()> {
    #[derive(serde::Serialize)]
    struct Receipt<'a> {
        version: u32,
        input: &'a crate::count_checkpoint::InputIdentity,
        output: crate::count_checkpoint::InputIdentity,
        tag: &'a str,
    }

    let parent = path.parent().unwrap_or(Path::new("."));
    let partial = path.with_extension(format!("partial-{:016x}", rand::random::<u64>()));
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&partial)?;
    serde_json::to_writer(
        &mut file,
        &Receipt {
            version: 1,
            input,
            output: crate::count_checkpoint::InputIdentity::read(output)?,
            tag,
        },
    )?;
    file.sync_all()?;
    std::fs::rename(partial, path)?;
    File::open(parent)?.sync_all()?;
    Ok(())
}

fn cleanup_owned_sort_work(
    work_dir: &Path,
    input: &crate::count_checkpoint::InputIdentity,
    tag: &str,
) {
    #[derive(serde::Deserialize)]
    struct SortState {
        version: u32,
        input: crate::count_checkpoint::InputIdentity,
        tag: String,
    }
    let state_path = work_dir.join("sort-state.json");
    if let Ok(state) = read_json_bounded::<SortState>(&state_path) {
        if state.version == 1 && state.input == *input && state.tag == tag {
            if let Err(error) = std::fs::remove_dir_all(work_dir) {
                log::warn!(
                    "Could not remove completed sort work {}: {error}",
                    work_dir.display()
                );
            }
        }
    }
}

#[cfg(all(test, feature = "bam"))]
mod tests {
    use super::*;
    use noodles_sam as sam;

    #[test]
    fn failed_merge_cleans_spill_files() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        write_bam(
            &input,
            &[("r0", Some("b")), ("r1", Some("a")), ("r2", Some("c"))],
        );
        let output = dir.path().join("is-a-directory");
        std::fs::create_dir(&output).unwrap();
        assert!(sort_bam_by_tag(&input, &output, "CB", 0, 1).is_err());
        assert!(std::fs::read_dir(dir.path())
            .unwrap()
            .flatten()
            .all(|e| !e.file_name().to_string_lossy().starts_with("cellsort_tmp_")));
    }

    #[test]
    fn publishes_a_completed_staged_sort_after_interruption() {
        use crate::count_checkpoint::InputIdentity;

        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        let work = output.with_extension("bam.sort-work");
        let staged = work.join("sorted.bam");
        write_bam(&input, &[("r0", Some("b")), ("r1", Some("a"))]);
        std::fs::create_dir_all(&work).unwrap();
        sort_bam_by_tag(&input, &staged, "CB", 0, 1).unwrap();
        let input_id = InputIdentity::read(&input).unwrap();
        let staged_id = InputIdentity::read(&staged).unwrap();
        std::fs::write(
            work.join("sort-state.json"),
            serde_json::to_vec(&serde_json::json!({
                "version": 1,
                "input": input_id,
                "tag": "CB",
                "output_bytes": staged_id.bytes,
                "output_modified_ns": staged_id.modified_ns
            }))
            .unwrap(),
        )
        .unwrap();

        prepare_sorted_bam(&input, &output, "CB", 1, 1).unwrap();

        assert_eq!(read_keys(&output).len(), 2);
        assert!(output.with_extension("bam.sort.json").exists());
        assert!(!work.exists());
    }

    #[test]
    fn retries_cleanly_after_an_interrupted_sort() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        let work = output.with_extension("bam.sort-work");
        write_bam(&input, &[("r0", Some("b")), ("r1", Some("a"))]);
        std::fs::create_dir_all(work.join("spills")).unwrap();
        std::fs::write(work.join("spills/chunk_0.bam"), b"incomplete").unwrap();
        std::fs::write(
            work.join("sort-state.json"),
            serde_json::to_vec(&serde_json::json!({
                "version": 1,
                "input": crate::count_checkpoint::InputIdentity::read(&input).unwrap(),
                "tag": "CB",
                "output_bytes": null,
                "output_modified_ns": null
            }))
            .unwrap(),
        )
        .unwrap();

        prepare_sorted_bam(&input, &output, "CB", 1, 1).unwrap();

        assert_eq!(read_keys(&output).len(), 2);
        assert!(!work.exists());
    }

    #[test]
    fn sort_budget_estimate_includes_cigar_and_aux_payloads() {
        use noodles_sam::alignment::record::cigar::{op::Kind, Op};

        let tag = cb_tag();
        let mut small = RecordBuf::default();
        small.data_mut().insert(tag, Value::String("CB".into()));
        small.cigar_mut().as_mut().push(Op::new(Kind::Match, 10));

        let mut large = RecordBuf::default();
        large
            .data_mut()
            .insert(tag, Value::String("C".repeat(4096).into()));
        large.cigar_mut().as_mut().extend([
            Op::new(Kind::Match, 10),
            Op::new(Kind::Skip, 20),
            Op::new(Kind::Match, 10),
        ]);

        assert!(record_footprint(&large, &tag) >= record_footprint(&small, &tag) + 4096);
    }

    #[test]
    fn many_runs_merge_with_bounded_fan_in() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        let names: Vec<_> = (0..70).map(|i| format!("read{i}")).collect();
        let records: Vec<_> = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), Some(if i % 2 == 0 { "b" } else { "a" })))
            .collect();
        write_bam(&input, &records);
        sort_bam_by_tag(&input, &output, "CB", 0, 1).unwrap();
        let keys = read_keys(&output);
        assert_eq!(keys.len(), 70);
        assert_contiguous(&keys);
    }

    #[test]
    fn sort_receipt_reuses_only_unchanged_complete_outputs() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        write_bam(&input, &[("r", Some("b"))]);
        prepare_sorted_bam(&input, &output, "CB", 1, 1).unwrap();
        let modified = std::fs::metadata(&output).unwrap().modified().unwrap();
        prepare_sorted_bam(&input, &output, "CB", 1, 1).unwrap();
        assert_eq!(
            std::fs::metadata(&output).unwrap().modified().unwrap(),
            modified
        );
        std::fs::OpenOptions::new()
            .write(true)
            .open(&output)
            .unwrap()
            .set_len(12)
            .unwrap();
        assert!(prepare_sorted_bam(&input, &output, "CB", 1, 1).is_err());
        assert_eq!(std::fs::metadata(&output).unwrap().len(), 12);
    }

    fn cb_tag() -> Tag {
        Tag::from([b'C', b'B'])
    }

    /// Write a BAM at `path` with one unmapped record per `(name, cb)`. `cb`
    /// of `None` writes a record with no CB tag.
    fn write_bam(path: &Path, records: &[(&str, Option<&str>)]) {
        let header = sam::Header::default();
        let mut writer = level7_bam_writer(File::create(path).unwrap(), 1);
        writer.write_header(&header).unwrap();
        let tag = cb_tag();
        for (name, cb) in records {
            let mut r = RecordBuf::default();
            *r.name_mut() = Some(name.as_bytes().to_vec().into());
            if let Some(cb) = cb {
                r.data_mut()
                    .insert(tag, Value::String(cb.to_string().into()));
            }
            writer.write_alignment_record(&header, &r).unwrap();
        }
        finish_bam(writer).unwrap();
    }

    /// Read back the CB key of every record, in file order.
    fn read_keys(path: &Path) -> Vec<SortKey> {
        let mut reader = File::open(path).map(bam::io::Reader::new).unwrap();
        let header = reader.read_header().unwrap();
        let tag = cb_tag();
        let mut out = Vec::new();
        let mut rec = RecordBuf::default();
        while reader.read_record_buf(&header, &mut rec).unwrap() != 0 {
            out.push(extract_key(&rec, &tag));
        }
        out
    }

    /// Assert no key value reappears after a different value intervened.
    fn assert_contiguous(keys: &[SortKey]) {
        let mut closed: Vec<SortKey> = Vec::new();
        let mut prev: Option<&SortKey> = None;
        for k in keys {
            if prev != Some(k) {
                if let Some(p) = prev {
                    closed.push(p.clone());
                }
                assert!(
                    !closed.contains(k),
                    "key reappeared non-contiguously after a different key"
                );
            }
            prev = Some(k);
        }
    }

    #[test]
    fn sorts_small_bam_contiguous_by_cb() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        write_bam(
            &input,
            &[
                ("r0", Some("AAA")),
                ("r1", Some("BBB")),
                ("r2", Some("AAA")),
                ("r3", Some("CCC")),
                ("r4", Some("BBB")),
                ("r5", Some("AAA")),
            ],
        );
        // Large budget → single-chunk fast path.
        sort_bam_by_tag(&input, &output, "CB", 1024, 2).unwrap();
        let keys = read_keys(&output);
        assert_eq!(keys.len(), 6);
        assert_contiguous(&keys);
    }

    #[test]
    fn spills_and_merges_when_budget_tiny() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        write_bam(
            &input,
            &[
                ("r0", Some("AAA")),
                ("r1", Some("BBB")),
                ("r2", Some("AAA")),
                ("r3", Some("CCC")),
                ("r4", Some("BBB")),
                ("r5", Some("AAA")),
            ],
        );
        // Budget 0 → at least one chunk per record → forces the merge path.
        sort_bam_by_tag(&input, &output, "CB", 0, 1).unwrap();
        let keys = read_keys(&output);
        assert_eq!(keys.len(), 6, "record count must be preserved");
        assert_contiguous(&keys);
    }

    #[test]
    fn records_missing_cb_form_trailing_group() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        write_bam(
            &input,
            &[
                ("r0", Some("BBB")),
                ("r1", None),
                ("r2", Some("AAA")),
                ("r3", None),
                ("r4", Some("BBB")),
            ],
        );
        sort_bam_by_tag(&input, &output, "CB", 0, 1).unwrap();
        let keys = read_keys(&output);
        assert_eq!(keys.len(), 5);
        assert_contiguous(&keys);
        // The two missing-tag records must be the final two.
        assert_eq!(keys[3], SortKey::Missing);
        assert_eq!(keys[4], SortKey::Missing);
    }

    #[test]
    fn empty_bam_roundtrips() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        let output = dir.path().join("out.bam");
        write_bam(&input, &[]);
        sort_bam_by_tag(&input, &output, "CB", 1024, 1).unwrap();
        assert!(read_keys(&output).is_empty());
    }

    #[test]
    fn key_ordering_total() {
        let a = SortKey::Present(b"AAA".to_vec());
        let b = SortKey::Present(b"BBB".to_vec());
        assert!(a < b);
        assert!(b < SortKey::Missing);
        assert!(a < SortKey::Missing);
        assert_eq!(a, SortKey::Present(b"AAA".to_vec()));
    }

    /// Single pass over a BAM: per-CB record counts, total records, and whether
    /// every CB value forms a single contiguous run (the only property the
    /// downstream counter relies on).
    fn cb_group_counts(path: &Path) -> (std::collections::HashMap<SortKey, u64>, u64, bool) {
        use std::collections::{HashMap, HashSet};
        let mut reader = File::open(path).map(bam::io::Reader::new).unwrap();
        let header = reader.read_header().unwrap();
        let tag = cb_tag();
        let mut counts: HashMap<SortKey, u64> = HashMap::new();
        let mut total: u64 = 0;
        let mut closed: HashSet<SortKey> = HashSet::new();
        let mut current: Option<SortKey> = None;
        let mut contiguous = true;
        let mut rec = RecordBuf::default();
        while reader.read_record_buf(&header, &mut rec).unwrap() != 0 {
            let key = extract_key(&rec, &tag);
            *counts.entry(key.clone()).or_insert(0) += 1;
            total += 1;
            if current.as_ref() != Some(&key) {
                if let Some(prev) = current.take() {
                    closed.insert(prev);
                }
                if closed.contains(&key) {
                    contiguous = false;
                }
                current = Some(key);
            }
        }
        (counts, total, contiguous)
    }

    /// End-to-end parity against a real `samtools sort -t CB` output. Compares
    /// the pure-Rust sort of the 1.8 GB fixture to the checked-in samtools
    /// ground truth: identical per-CB read counts, identical total, and CB
    /// contiguity in both. Ignored by default (needs the large fixtures and
    /// minutes of runtime); run with:
    ///   cargo test --release --features bam -- --ignored --nocapture e2e_matches_samtools
    #[test]
    #[ignore = "requires the large real BAM fixtures in tests/; run with --ignored"]
    fn e2e_matches_samtools_grouping() {
        let manifest = env!("CARGO_MANIFEST_DIR");
        let input = Path::new(manifest).join("tests/p395aWT1_velocyto_sample_tag_04.bam");
        let truth =
            Path::new(manifest).join("tests/cellsorted_p395aWT1_velocyto_sample_tag_04.bam");
        if !input.exists() || !truth.exists() {
            eprintln!("skipping e2e: fixtures not present ({})", input.display());
            return;
        }

        let dir = tempfile::tempdir().unwrap();
        let mine = dir.path().join("mine_cellsorted.bam");
        // ~1 GB per-chunk budget on a 1.8 GB BAM → forces the spill+merge path.
        // Time the sort alone (the read-back below dominates total test time and
        // is single-threaded, so it would mask the multithreaded-writer speedup).
        let t0 = std::time::Instant::now();
        sort_bam_by_tag(&input, &mine, "CB", 1024, 16).unwrap();
        eprintln!(
            "pure-Rust sort (16 BGZF workers) wall time: {:.1} s",
            t0.elapsed().as_secs_f64()
        );

        let (mine_counts, mine_total, mine_contig) = cb_group_counts(&mine);
        let (truth_counts, truth_total, truth_contig) = cb_group_counts(&truth);

        eprintln!(
            "rust: {mine_total} reads / {} barcodes (contiguous={mine_contig})",
            mine_counts.len()
        );
        eprintln!(
            "samtools: {truth_total} reads / {} barcodes (contiguous={truth_contig})",
            truth_counts.len()
        );

        assert!(mine_contig, "pure-Rust output is not contiguous by CB");
        assert!(
            truth_contig,
            "samtools ground truth is not contiguous by CB"
        );
        assert_eq!(
            mine_total, truth_total,
            "total record count differs (rust {mine_total} vs samtools {truth_total})"
        );
        assert_eq!(
            mine_counts, truth_counts,
            "per-CB read counts differ between rust and samtools"
        );
    }
}
