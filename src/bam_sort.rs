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

/// Rough resident-memory estimate for one record, used only to bound a chunk to
/// the `-m` budget. The dominant variable cost is the sequence plus its quality
/// scores (both ~one byte per base); the constant covers struct/Vec overhead,
/// name, cigar, and aux data. Order-of-magnitude accuracy is sufficient.
fn record_footprint(record: &RecordBuf) -> usize {
    let name = record.name().map_or(0, |n| n.len());
    let seq = record.sequence().len();
    256 + name + 2 * seq
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
/// Spill files live persistently in `output`'s parent directory, including on failure.
pub fn sort_bam_by_tag(
    input: &Path,
    output: &Path,
    tag: &str,
    mem_budget_mb: usize,
    threads: usize,
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

    let out_parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));

    // ── Spill phase ───────────────────────────────────────────────────────────
    let mut record = RecordBuf::default();
    let mut chunk: Vec<RecordBuf> = Vec::new();
    let mut spilled: Vec<std::path::PathBuf> = Vec::new();
    let mut tmpdir: Option<std::path::PathBuf> = None;
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
            chunk_bytes += record_footprint(&record);
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
            if tmpdir.is_none() {
                tmpdir = Some(
                    tempfile::Builder::new()
                        .prefix("cellsort_tmp_")
                        .tempdir_in(out_parent)
                        .with_context(|| {
                            format!("creating spill temp dir in {}", out_parent.display())
                        })?
                        .keep(),
                );
            }
            let path = tmpdir
                .as_ref()
                .unwrap()
                .join(format!("chunk_{}.bam", spilled.len()));
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
            let manifest = tmpdir.as_ref().unwrap().join("sort-runs.json");
            let partial = manifest.with_extension("partial");
            let mut receipt = File::create(&partial)?;
            serde_json::to_writer(&mut receipt, &spilled)?;
            receipt.sync_all()?;
            std::fs::rename(partial, manifest)?;
            File::open(tmpdir.as_ref().unwrap())?.sync_all()?;
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
            let path = tmpdir
                .as_ref()
                .unwrap()
                .join(format!("merge_{pass}_{i}.bam"));
            merge_runs(runs, &path, &header, &tag, threads)?;
            File::open(&path)?.sync_all()?;
            next.push(path);
        }
        spilled = next;
        pass += 1;
    }
    merge_runs(&spilled, output, &header, &tag, threads)?;
    log::info!(
        "Sort spill files retained at {}",
        tmpdir.as_ref().unwrap().display()
    );
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

/// Reuse only a published, input-bound sort; incomplete files are retained.
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
    if output.exists() {
        let receipt_file = File::open(&receipt_path).with_context(|| format!(
            "Unverified sorted BAM {}; retained unchanged. Move it aside before sorting, or resume from an existing counting checkpoint using write-output", output.display()))?;
        anyhow::ensure!(
            receipt_file.metadata()?.len() <= 1_048_576,
            "sort receipt too large"
        );
        let receipt: Receipt = serde_json::from_reader(receipt_file)?;
        anyhow::ensure!(
            receipt
                == Receipt {
                    version: 1,
                    input: input_identity,
                    output: InputIdentity::read(output)?,
                    tag: tag.into()
                },
            "sorted BAM provenance differs; retained unchanged"
        );
        return Ok(());
    }
    let name = output
        .file_name()
        .context("sorted BAM needs a filename")?
        .to_string_lossy();
    let partial = output.with_file_name(format!("{name}.partial-{:016x}", rand::random::<u64>()));
    sort_bam_by_tag(input, &partial, tag, memory_mb, threads).with_context(|| {
        format!(
            "Sort failed; intermediate files and {} retained",
            partial.display()
        )
    })?;
    File::open(&partial)?.sync_all()?;
    std::fs::rename(&partial, output)?;
    let receipt = Receipt {
        version: 1,
        input: input_identity,
        output: InputIdentity::read(output)?,
        tag: tag.into(),
    };
    let partial_receipt =
        receipt_path.with_extension(format!("partial-{:016x}", rand::random::<u64>()));
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&partial_receipt)?;
    serde_json::to_writer(&mut file, &receipt)?;
    file.sync_all()?;
    std::fs::rename(partial_receipt, receipt_path)?;
    File::open(
        output
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new(".")),
    )?
    .sync_all()?;
    Ok(())
}

#[cfg(all(test, feature = "bam"))]
mod tests {
    use super::*;
    use noodles_sam as sam;

    #[test]
    fn failed_merge_keeps_spill_files() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.bam");
        write_bam(
            &input,
            &[("r0", Some("b")), ("r1", Some("a")), ("r2", Some("c"))],
        );
        let output = dir.path().join("is-a-directory");
        std::fs::create_dir(&output).unwrap();
        assert!(sort_bam_by_tag(&input, &output, "CB", 0, 1).is_err());
        let retained: Vec<_> = std::fs::read_dir(dir.path())
            .unwrap()
            .flatten()
            .filter(|e| e.file_name().to_string_lossy().starts_with("cellsort_tmp_"))
            .collect();
        assert_eq!(
            retained.len(),
            1,
            "failed sorts must retain their spill directory"
        );
        assert!(std::fs::read_dir(retained[0].path()).unwrap().count() >= 3);
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
