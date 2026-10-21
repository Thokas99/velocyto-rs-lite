//! Translated from velocyto/commands/_run.py and velocyto/commands/run.py

use std::collections::HashSet;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Command, Stdio};

use anyhow::{bail, Context};
use clap::Args;
use rand::Rng;

use crate::constants::BAM_COMPRESSION;
use crate::counter::ExInCounter;
use crate::logic::logic_from_name;

#[derive(Args)]
pub struct CountingOptions {
    /// Persistent count checkpoint directory (default: OUTPUT/SAMPLE.counts)
    #[arg(long, conflicts_with = "resume_from")]
    pub checkpoint_dir: Option<std::path::PathBuf>,
    /// Resume counting from this checkpoint (inputs and settings must match)
    #[arg(long)]
    pub resume_from: Option<std::path::PathBuf>,
    /// Maximum raw barcodes per count batch; groups are never split
    #[arg(long, default_value_t = 100)]
    pub cell_batch_size: usize,
    /// MiB budget per count read/evidence/scratch component (not total job RSS)
    #[arg(long, default_value_t = 512)]
    pub count_memory_mb: usize,
}
impl Default for CountingOptions {
    fn default() -> Self {
        Self {
            checkpoint_dir: None,
            resume_from: None,
            cell_batch_size: 100,
            count_memory_mb: 512,
        }
    }
}

/// Arguments for the generic `run` subcommand.
///
/// Accepts one or more BAM files and a genome annotation GTF. All other
/// options give fine-grained control over barcodes, logic, samtools
/// parallelism, and output format.
#[derive(Args)]
pub struct RunArgs {
    #[command(flatten)]
    pub counting: CountingOptions,
    /// Enable verbose (debug) logging
    #[arg(short, long)]
    pub verbose: bool,

    /// One or more BAM files (position-sorted)
    #[arg(required = true)]
    pub bamfile: Vec<String>,
    /// Genome annotation GTF file
    pub gtffile: String,
    /// Valid barcodes file (plain text or .gz), one barcode per line
    #[arg(short = 'b', long)]
    pub bcfile: Option<String>,
    /// Output folder (created if absent; defaults to `<bam-dir>/velocyto`)
    #[arg(short = 'o', long)]
    pub outputfolder: Option<String>,
    /// Sample name used as the output loom filename stem
    #[arg(short = 'e', long)]
    pub sampleid: Option<String>,
    /// CSV metadata table; rows = samples, cols = entries
    #[arg(short = 's', long)]
    pub metadatatable: Option<String>,
    /// GTF file containing genomic intervals to mask (e.g. repeats)
    #[arg(short = 'm', long)]
    pub mask: Option<String>,
    /// Treat each input BAM as a single cell (SmartSeq2-style)
    #[arg(short = 'c', long, default_value_t = false)]
    pub onefilepercell: bool,
    /// Molecule-filtering logic class name (default: Default)
    #[arg(short = 'l', long, default_value = "Default")]
    pub logic: String,
    /// Do not use UMI information (read-count mode)
    #[arg(short = 'U', long, default_value_t = false)]
    pub without_umi: bool,
    /// Extend UMI identity to include genomic position context
    /// (`no`, `chr`, `Gene`, `Cluster`, `all`)
    #[arg(short = 'u', long, default_value = "no")]
    pub umi_extension: String,
    /// Consider non-unique mappings (not recommended)
    #[arg(short = 'M', long, default_value_t = false)]
    pub multimap: bool,
    /// Number of threads for the built-in Rust barcode sort
    #[arg(long, default_value_t = 16)]
    pub samtools_threads: usize,
    /// MiB per sort chunk for the built-in Rust barcode sort
    #[arg(long, default_value_t = 2048)]
    pub samtools_memory: usize,
    /// dtype for loom layer arrays: "uint32" (default, lossless) or "uint16" (smaller, saturates at 65535)
    #[arg(short = 't', long, default_value = "uint32")]
    pub dtype: String,
    /// Output file format: "h5ad" (default, AnnData), "loom", or "both"
    #[arg(long, default_value = "h5ad")]
    pub output_format: String,
    /// Debug dump: save a molecular mapping report every N cells (0 = disabled)
    #[arg(short = 'd', long, default_value = "0")]
    pub dump: String,
    /// BAM tag for cell barcode (overrides auto-detection; e.g. `CB` or `XC`)
    #[arg(long)]
    pub cb_tag: Option<String>,
    /// BAM tag for UMI barcode (overrides auto-detection; e.g. `UB` or `XM`)
    #[arg(long)]
    pub ub_tag: Option<String>,
    /// BAM tag carrying the sample identity for demultiplexing a multi-sample
    /// BAM in place (e.g. BD Rhapsody `ST`). When set, the cell identity becomes
    /// `(sample, barcode)`, CellIDs are formatted `{sampleid}_{sample}:{bc}`,
    /// and a `SampleID` column attribute is added. Omit for single-sample BAMs.
    #[arg(long)]
    pub sample_tag: Option<String>,
}

/// Runs the velocity analysis for generic BAM input.
///
/// Thin wrapper that forwards all [`RunArgs`] fields to [`run_inner`].
pub fn run(args: RunArgs) -> anyhow::Result<()> {
    run_inner_with_options(
        &args.bamfile,
        &args.gtffile,
        args.bcfile.as_deref(),
        args.outputfolder.as_deref(),
        args.sampleid.as_deref(),
        args.mask.as_deref(),
        args.onefilepercell,
        &args.logic,
        args.without_umi,
        &args.umi_extension,
        args.multimap,
        args.samtools_threads,
        args.samtools_memory,
        &args.dump,
        &args.dtype,
        &[],
        args.cb_tag.as_deref(),
        args.ub_tag.as_deref(),
        args.sample_tag.as_deref(),
        &args.output_format,
        &args.counting,
    )
}

/// Generate a random alphanumeric identifier of the given length.
///
/// Characters are drawn from `A-Z0-9`. Used to produce unique suffixes for
/// auto-generated sample IDs so that repeated runs do not overwrite each
/// other's output files.
///
/// Translated from `_run.id_generator` in Python, which used
/// `random.choice(string.ascii_uppercase + string.digits)`.
pub fn id_generator(size: usize) -> String {
    let mut rng = rand::thread_rng();
    (0..size)
        .map(|_| {
            let idx = rng.gen_range(0..36usize);
            if idx < 26 {
                (b'A' + idx as u8) as char
            } else {
                (b'0' + (idx - 26) as u8) as char
            }
        })
        .collect()
}

/// Core velocity analysis pipeline shared by all subcommands.
///
/// Translated from `_run._run` in Python. All higher-level entry points
/// (`run10x`, `run_dropest`, `run_smartseq2`, `run`) call this function.
///
/// # Pipeline steps
/// 1. **Resolve inputs** — validate BAM/barcode combinations; auto-generate
///    `sampleid` and `outputfolder` when not provided.
/// 2. **Parse barcodes** — read the whitelist from `bcfile` (plain or `.gz`);
///    strip GEM-group suffixes; derive `gem_grp` (`-N` when all barcodes share
///    one group, `x` when mixed, empty string when no whitelist).
/// 3. **ExInCounter setup** — instantiate the counter with the chosen `logic`
///    and barcode set.
/// 4. **Peek BAM** — inspect the first 1000 records to discover cell-barcode
///    and UMI tag names.
/// 5. **Launch samtools sort** (async) — sort each BAM by the cell-barcode tag
///    (`cellsorted_<bam>`) while annotation loading proceeds in parallel.
/// 6. **Load GTF annotations** — parse transcript models and, optionally, a
///    repeat-mask GTF.
/// 7. **Mark up introns** — first-pass scan of all BAMs to validate intron
///    intervals against actual read splicing patterns.
/// 8. **Wait for samtools** — block until all sort processes finish.
/// 9. **Count molecules** — second-pass scan of the cell-sorted BAMs to
///    accumulate per-cell, per-gene spliced/unspliced/ambiguous counts.
/// 10. **Write loom** — call [`ExInCounter::dump_loom`] to produce the output
///     `.loom` file.
///
/// # Differences from Python
/// - `additional_ca` float column attributes (tSNE, clusters) are accepted but
///   silently ignored — hdf5-pure-rs does not expose a float-array attribute
///   API.
/// - Memory and CPU availability are read directly from `/proc/meminfo` and
///   `/proc/cpuinfo` rather than via `subprocess` / `multiprocessing`.
/// - The `test` / pickle-debug shortcut present in the Python code is omitted.
#[cfg(feature = "bam")]
pub fn run_inner(
    bamfile: &[String],
    gtffile: &str,
    bcfile: Option<&str>,
    outputfolder: Option<&str>,
    sampleid: Option<&str>,
    repmask: Option<&str>,
    onefilepercell: bool,
    logic: &str,
    without_umi: bool,
    umi_extension: &str,
    multimap: bool,
    samtools_threads: usize,
    samtools_memory: usize,
    dump: &str,
    loom_numeric_dtype: &str,
    additional_ca: &[(&str, Vec<f32>)],
    cb_tag: Option<&str>,
    ub_tag: Option<&str>,
    sample_tag: Option<&str>,
    output_format: &str,
) -> anyhow::Result<()> {
    run_inner_with_options(
        bamfile,
        gtffile,
        bcfile,
        outputfolder,
        sampleid,
        repmask,
        onefilepercell,
        logic,
        without_umi,
        umi_extension,
        multimap,
        samtools_threads,
        samtools_memory,
        dump,
        loom_numeric_dtype,
        additional_ca,
        cb_tag,
        ub_tag,
        sample_tag,
        output_format,
        &CountingOptions::default(),
    )
}

#[cfg(feature = "bam")]
fn run_inner_with_options(
    bamfile: &[String],
    gtffile: &str,
    bcfile: Option<&str>,
    outputfolder: Option<&str>,
    sampleid: Option<&str>,
    repmask: Option<&str>,
    onefilepercell: bool,
    logic: &str,
    without_umi: bool,
    umi_extension: &str,
    multimap: bool,
    samtools_threads: usize,
    samtools_memory: usize,
    dump: &str,
    loom_numeric_dtype: &str,
    additional_ca: &[(&str, Vec<f32>)],
    cb_tag: Option<&str>,
    ub_tag: Option<&str>,
    sample_tag: Option<&str>,
    output_format: &str,
    counting: &CountingOptions,
) -> anyhow::Result<()> {
    // ── Resolve inputs ────────────────────────────────────────────────────────
    validate_loom_dtype(loom_numeric_dtype)?;
    validate_output_format(output_format)?;

    anyhow::ensure!(!bamfile.is_empty(), "at least one BAM is required");
    anyhow::ensure!(
        counting.cell_batch_size > 0 && counting.count_memory_mb > 0,
        "counting limits must be positive"
    );
    let count_bytes = counting
        .count_memory_mb
        .checked_mul(1024 * 1024)
        .context("count budget overflow")?;
    let previous = counting
        .resume_from
        .as_ref()
        .map(|p| crate::count_checkpoint::Checkpoint::open(p))
        .transpose()?;
    let sampleid = sampleid.or_else(|| {
        previous
            .as_ref()
            .map(|c| c.manifest.metadata.sample.as_str())
    });
    let multi = bamfile.len() > 1;

    if onefilepercell && multi && bcfile.is_some() {
        bail!("Inputs incompatibility. --bcfile/-b used together with --onefilepercell/-c.");
    }

    let sampleid: String = match sampleid {
        Some(s) => s.to_string(),
        None => {
            let stem0 = Path::new(&bamfile[0])
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();
            if multi && !onefilepercell {
                let full_name: String = bamfile
                    .iter()
                    .map(|f| {
                        Path::new(f)
                            .file_stem()
                            .unwrap_or_default()
                            .to_string_lossy()
                            .to_string()
                    })
                    .collect::<Vec<_>>()
                    .join("_");
                if full_name.len() > 50 {
                    format!("multi_input_{}_{}", stem0, id_generator(5))
                } else {
                    format!("multi_input_{}_and_others_{}", full_name, id_generator(5))
                }
            } else if multi && onefilepercell {
                format!("onefilepercell_{}_and_others_{}", stem0, id_generator(5))
            } else {
                format!("{}_{}", stem0, id_generator(5))
            }
        }
    };

    let outputfolder: String = match outputfolder {
        Some(o) => o.to_string(),
        None => {
            let parent = Path::new(&bamfile[0]).parent().unwrap_or(Path::new("."));
            parent.join("velocyto").to_string_lossy().to_string()
        }
    };
    if !Path::new(&outputfolder).exists() {
        fs::create_dir_all(&outputfolder)
            .with_context(|| format!("Creating output folder {outputfolder}"))?;
    }

    // ── Barcodes ──────────────────────────────────────────────────────────────
    let (valid_bcset, gem_grp) = if let Some(bc_path) = bcfile {
        let content = if bc_path.ends_with(".gz") {
            let f = fs::File::open(bc_path).with_context(|| format!("Opening {bc_path}"))?;
            let mut gz = flate2::read::GzDecoder::new(f);
            let mut s = String::new();
            std::io::Read::read_to_string(&mut gz, &mut s)?;
            s
        } else {
            fs::read_to_string(bc_path).with_context(|| format!("Reading {bc_path}"))?
        };
        let bcs: Vec<String> = content.split_whitespace().map(|s| s.to_string()).collect();
        let gem_grp_val = {
            let suffixes: HashSet<&str> =
                bcs.iter().filter_map(|bc| bc.split('-').last()).collect();
            if suffixes.len() == 1 {
                format!("-{}", bcs[0].split('-').last().unwrap_or("1"))
            } else {
                "x".to_string()
            }
        };
        let bcset: HashSet<String> = bcs
            .iter()
            .map(|bc| bc.splitn(2, '-').next().unwrap_or(bc).to_string())
            .collect();
        log::info!("Read {} cell barcodes from {bc_path}", bcset.len());
        (Some(bcset), gem_grp_val)
    } else {
        (None, String::new())
    };

    // ── ExInCounter setup ─────────────────────────────────────────────────────
    let umi_ext = if without_umi {
        "without_umi"
    } else {
        umi_extension
    };
    let logic_box = logic_from_name(logic);

    let mut exincounter = ExInCounter::new(
        sampleid.clone(),
        logic_box,
        valid_bcset,
        umi_ext,
        onefilepercell,
        dump,
        outputfolder.clone(),
        loom_numeric_dtype.to_string(),
    )?;

    // ── Memory / thread heuristic ─────────────────────────────────────────────
    let mb_available = read_mem_available_mb();
    let threads_to_use = samtools_threads.max(1).min(num_cpus());
    let mb_to_use = (samtools_memory as u64)
        .min((mb_available / 2).max(1))
        .max(1) as usize;

    // ── Peek BAM for tag names ────────────────────────────────────────────────
    // When both tags are supplied explicitly, skip the 1000-read peek entirely.
    // When only one (or neither) is supplied, peek auto-detects both, then the
    // explicit values override whichever the caller specified.
    let tagname = if onefilepercell && without_umi {
        "NOTAG".to_string()
    } else if onefilepercell {
        if ub_tag.is_none() {
            exincounter.peek_umi_only(&bamfile[0], 1000)?;
        }
        "NOTAG".to_string()
    } else {
        if cb_tag.is_none() || ub_tag.is_none() {
            exincounter.peek(&bamfile[0], 1000)?;
        }
        cb_tag.unwrap_or(&exincounter.cellbarcode_str).to_string()
    };
    if let Some(cb) = cb_tag {
        exincounter.cellbarcode_str = cb.to_string();
    }
    if let Some(ub) = ub_tag {
        exincounter.umibarcode_str = ub.to_string();
    }
    if let Some(st) = sample_tag {
        log::info!("Sample demultiplexing enabled on BAM tag '{st}'");
        exincounter.sample_tag = Some(st.to_string());
    }

    // ── Cell-sorted BAM paths ─────────────────────────────────────────────────
    let bamfile_cellsorted: Vec<String> = if multi && onefilepercell {
        bamfile.to_vec()
    } else if onefilepercell {
        vec![bamfile[0].clone()]
    } else {
        bamfile
            .iter()
            .map(|bmf| {
                let parent = Path::new(bmf).parent().unwrap_or(Path::new("."));
                let basename = Path::new(bmf).file_name().unwrap_or_default();
                parent
                    .join(format!("cellsorted_{}", basename.to_string_lossy()))
                    .to_string_lossy()
                    .to_string()
            })
            .collect()
    };

    // Sort sequentially before annotation allocation; persistent Rust spills survive errors.
    if tagname != "NOTAG" {
        for (input, output) in bamfile.iter().zip(&bamfile_cellsorted) {
            crate::bam_sort::prepare_sorted_bam(
                Path::new(input),
                Path::new(output),
                &tagname,
                mb_to_use,
                threads_to_use,
            )?;
        }
    }

    // ── Load annotations after sorting ────────────────────────────────
    log::info!("Loading annotation from {gtffile}");
    exincounter.read_transcriptmodels(gtffile)?;

    if let Some(mask_path) = repmask {
        log::info!("Loading repeat mask from {mask_path}");
        exincounter.read_repeats(mask_path, 0)?;
    }

    // ── Mark up introns ───────────────────────────────────────────────────────
    log::info!("Scanning BAM(s) to validate intron intervals");
    for bmf in bamfile.iter() {
        exincounter.mark_up_introns(bmf, multimap)?;
    }

    // ── Durable counts and bounded output ────────────────────────────────────
    use crate::count_checkpoint::{Checkpoint, InputIdentity, Metadata};
    let identities = bamfile_cellsorted
        .iter()
        .map(|p| InputIdentity::read(Path::new(p)))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let mut provenance_paths: Vec<&str> = bamfile.iter().map(String::as_str).collect();
    provenance_paths.push(gtffile);
    provenance_paths.extend(bcfile);
    provenance_paths.extend(repmask);
    let provenance = provenance_paths
        .iter()
        .map(|p| InputIdentity::read(Path::new(p)))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let settings = serde_json::to_string(&(
        env!("CARGO_PKG_VERSION"),
        logic,
        multimap,
        umi_ext,
        onefilepercell,
        &exincounter.cellbarcode_str,
        &exincounter.umibarcode_str,
        sample_tag,
        provenance,
    ))?;
    let gem = if exincounter.filter_mode {
        gem_grp
    } else {
        String::new()
    };
    let metadata = Metadata::from_counter(&exincounter, gem, identities, settings)?;
    let root = counting
        .resume_from
        .clone()
        .or_else(|| counting.checkpoint_dir.clone())
        .unwrap_or_else(|| Path::new(&outputfolder).join(format!("{sampleid}.counts")));
    let mut checkpoint = if counting.resume_from.is_some() {
        Checkpoint::resume(&root, &metadata)?
    } else {
        Checkpoint::create(&root, metadata)?
    };
    log::info!("Persistent counting checkpoint: {}", root.display());
    for (source, bam) in bamfile_cellsorted.iter().enumerate() {
        if checkpoint.manifest.sources[source].complete {
            continue;
        }
        let offset = checkpoint.manifest.sources[source].offset;
        exincounter.count_stream(
            bam,
            multimap,
            counting.cell_batch_size,
            offset,
            count_bytes,
            |layers, keys, next| {
                checkpoint.commit_batch(source, next, &layers, &keys)?;
                log::info!(
                    "Committed count batch {} ({} cells)",
                    checkpoint.manifest.batches.len(),
                    keys.len()
                );
                Ok(())
            },
        )?;
        checkpoint.finish_source(source)?;
    }
    // Annotation and barcode buffers are not needed by checkpoint-backed export.
    drop(exincounter);
    crate::commands::write_output::export(
        &checkpoint,
        Path::new(&outputfolder),
        &sampleid,
        output_format,
    )?;

    let _ = additional_ca; // float col attrs not yet supported by hdf5-pure-rust string-only API

    log::info!("Terminated successfully!");
    Ok(())
}

/// Stub that errors when the crate is built without the `bam` feature.
///
/// Recompile with `cargo build --features bam` to enable BAM support.
#[cfg(not(feature = "bam"))]
pub fn run_inner(
    _bamfile: &[String],
    _gtffile: &str,
    _bcfile: Option<&str>,
    _outputfolder: Option<&str>,
    _sampleid: Option<&str>,
    _repmask: Option<&str>,
    _onefilepercell: bool,
    _logic: &str,
    _without_umi: bool,
    _umi_extension: &str,
    _multimap: bool,
    _samtools_threads: usize,
    _samtools_memory: usize,
    _dump: &str,
    _loom_numeric_dtype: &str,
    _additional_ca: &[(&str, Vec<f32>)],
    _cb_tag: Option<&str>,
    _ub_tag: Option<&str>,
    _sample_tag: Option<&str>,
    _output_format: &str,
) -> anyhow::Result<()> {
    anyhow::bail!("BAM support required. Recompile with: cargo build --features bam")
}

#[cfg(not(feature = "bam"))]
fn run_inner_with_options(
    bamfile: &[String],
    gtffile: &str,
    bcfile: Option<&str>,
    outputfolder: Option<&str>,
    sampleid: Option<&str>,
    repmask: Option<&str>,
    onefilepercell: bool,
    logic: &str,
    without_umi: bool,
    umi_extension: &str,
    multimap: bool,
    samtools_threads: usize,
    samtools_memory: usize,
    dump: &str,
    loom_numeric_dtype: &str,
    additional_ca: &[(&str, Vec<f32>)],
    cb_tag: Option<&str>,
    ub_tag: Option<&str>,
    sample_tag: Option<&str>,
    output_format: &str,
    counting: &CountingOptions,
) -> anyhow::Result<()> {
    anyhow::bail!("BAM support required. Recompile with: cargo build --features bam")
}

/// Read the `MemAvailable` field from `/proc/meminfo` and return it in MB.
///
/// Falls back to 32 000 MB (≈ 32 GB) when `/proc/meminfo` is unavailable
/// (non-Linux systems) or the field cannot be parsed.
#[cfg(feature = "bam")]
fn read_mem_available_mb() -> u64 {
    let host = fs::read_to_string("/proc/meminfo")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|line| line.starts_with("MemAvailable:"))?
                .split_whitespace()
                .nth(1)?
                .parse::<u64>()
                .ok()
        })
        .map(|kb| kb / 1024)
        .unwrap_or(32_000);
    let env_number = |name| std::env::var(name).ok().and_then(|s| s.parse::<u64>().ok());
    let slurm = env_number("SLURM_MEM_PER_NODE").or_else(|| {
        env_number("SLURM_MEM_PER_CPU")?.checked_mul(
            env_number("SLURM_CPUS_PER_TASK").or_else(|| env_number("SLURM_CPUS_ON_NODE"))?,
        )
    });
    let mut available = host.min(slurm.unwrap_or(host));
    if let Ok(groups) = fs::read_to_string("/proc/self/cgroup") {
        if let Some(relative) = groups.lines().find_map(|line| line.strip_prefix("0::")) {
            let mount = Path::new("/sys/fs/cgroup");
            let mut current = mount.join(relative.trim_start_matches('/'));
            loop {
                let read_number = |name| {
                    fs::read_to_string(current.join(name))
                        .ok()?
                        .trim()
                        .parse::<u64>()
                        .ok()
                };
                if let Some(used) = read_number("memory.current") {
                    for name in ["memory.max", "memory.high"] {
                        if let Some(limit) = read_number(name) {
                            available = available.min(limit.saturating_sub(used) / (1024 * 1024));
                        }
                    }
                }
                if current == mount || !current.pop() {
                    break;
                }
            }
        }
    }
    available
}

/// Count logical CPU cores by counting `processor` entries in `/proc/cpuinfo`.
///
/// Falls back to 1 when `/proc/cpuinfo` is unavailable (non-Linux systems).
#[cfg(feature = "bam")]
fn num_cpus() -> usize {
    let affinity = std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1);
    ["SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"]
        .iter()
        .filter_map(|name| std::env::var(name).ok()?.parse::<usize>().ok())
        .fold(affinity, usize::min)
        .max(1)
}

/// Validate the `--dtype` value for loom layer datasets.
///
/// Only `"uint32"` (default, lossless) and `"uint16"` (narrower, saturates at
/// 65535) are supported. Any other value — including typos like `"unit32"` —
/// is rejected rather than silently coerced.
#[cfg(feature = "bam")]
fn validate_loom_dtype(dtype: &str) -> anyhow::Result<()> {
    if matches!(dtype, "uint16" | "uint32") {
        Ok(())
    } else {
        bail!("Invalid --dtype '{dtype}'. Supported values: \"uint32\" (default) or \"uint16\".")
    }
}

/// Validate the `--output-format` value.
///
/// Only `"h5ad"` (default, AnnData), `"loom"`, and `"both"` are supported. Any
/// other value is rejected rather than silently coerced.
#[cfg(feature = "bam")]
fn validate_output_format(format: &str) -> anyhow::Result<()> {
    if matches!(format, "h5ad" | "loom" | "both") {
        Ok(())
    } else {
        bail!(
            "Invalid --output-format '{format}'. Supported values: \"h5ad\" (default), \
             \"loom\", or \"both\"."
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // ── id_generator ──────────────────────────────────────────────────────────

    #[test]
    fn id_generator_length() {
        for size in [5, 6, 10] {
            assert_eq!(id_generator(size).len(), size);
        }
    }

    #[test]
    fn id_generator_charset() {
        let id = id_generator(1000);
        assert!(id
            .chars()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit()));
    }

    #[test]
    fn id_generator_zero_length() {
        assert_eq!(id_generator(0), "");
    }

    // ── read_mem_available_mb ─────────────────────────────────────────────────

    #[test]
    #[cfg(feature = "bam")]
    fn read_mem_available_mb_returns_nonzero_on_linux() {
        // /proc/meminfo exists on Linux; fallback is 32000 which is also valid
        let mb = read_mem_available_mb();
        assert!(mb > 0, "expected positive MB value, got {mb}");
    }

    // ── num_cpus ──────────────────────────────────────────────────────────────

    #[test]
    #[cfg(feature = "bam")]
    fn num_cpus_returns_at_least_one() {
        assert!(num_cpus() >= 1);
    }

    // ── barcode parsing (gem_grp logic extracted for testing) ─────────────────

    fn parse_gem_grp(bcs: &[&str]) -> String {
        let suffixes: std::collections::HashSet<&str> =
            bcs.iter().filter_map(|bc| bc.split('-').last()).collect();
        if suffixes.len() == 1 {
            format!("-{}", bcs[0].split('-').last().unwrap_or("1"))
        } else {
            "x".to_string()
        }
    }

    fn parse_bcset(bcs: &[&str]) -> std::collections::HashSet<String> {
        bcs.iter()
            .map(|bc| bc.splitn(2, '-').next().unwrap_or(bc).to_string())
            .collect()
    }

    #[test]
    fn gem_grp_single_suffix_extracts_dash_suffix() {
        let bcs = ["ACGT-1", "TTTT-1", "GGGG-1"];
        assert_eq!(parse_gem_grp(&bcs), "-1");
    }

    #[test]
    fn gem_grp_multiple_suffixes_returns_x() {
        let bcs = ["ACGT-1", "TTTT-2"];
        assert_eq!(parse_gem_grp(&bcs), "x");
    }

    #[test]
    fn gem_grp_no_suffix_returns_x() {
        let bcs = ["ACGTACGT", "TTTTGGGG"];
        assert_eq!(parse_gem_grp(&bcs), "x");
    }

    #[test]
    fn bcset_strips_gem_group_suffix() {
        let bcs = ["ACGT-1", "TTTT-1"];
        let set = parse_bcset(&bcs);
        assert!(set.contains("ACGT"));
        assert!(set.contains("TTTT"));
        assert!(!set.iter().any(|s| s.contains('-')));
    }

    #[test]
    fn bcset_no_suffix_keeps_full_barcode() {
        let bcs = ["ACGTACGT"];
        let set = parse_bcset(&bcs);
        assert!(set.contains("ACGTACGT"));
    }

    // ── sampleid auto-generation ──────────────────────────────────────────────

    fn derive_sampleid(bamfile: &[&str], multi: bool, onefilepercell: bool) -> String {
        let stem0 = Path::new(bamfile[0])
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        if multi && !onefilepercell {
            let full_name: String = bamfile
                .iter()
                .map(|f| {
                    Path::new(f)
                        .file_stem()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .to_string()
                })
                .collect::<Vec<_>>()
                .join("_");
            if full_name.len() > 50 {
                format!("multi_input_{stem0}_XXXXX")
            } else {
                format!("multi_input_{full_name}_and_others_XXXXX")
            }
        } else if multi && onefilepercell {
            format!("onefilepercell_{stem0}_XXXXX")
        } else {
            format!("{stem0}_XXXXX")
        }
    }

    #[test]
    fn sampleid_single_bam_uses_stem() {
        let id = derive_sampleid(&["path/to/sample.bam"], false, false);
        assert!(id.starts_with("sample_"), "got: {id}");
    }

    #[test]
    fn sampleid_multi_bam_uses_multi_prefix() {
        let id = derive_sampleid(&["a.bam", "b.bam"], true, false);
        assert!(id.starts_with("multi_input_"), "got: {id}");
    }

    #[test]
    fn sampleid_onefilepercell_uses_onefilepercell_prefix() {
        let id = derive_sampleid(&["cell1.bam", "cell2.bam"], true, true);
        assert!(id.starts_with("onefilepercell_"), "got: {id}");
    }

    // ── load_mapping_tsv (dropest_bc_correct) ─────────────────────────────────

    #[test]
    fn load_mapping_tsv_parses_two_column_tsv() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "ORIGBC1\tCORRBC1").unwrap();
        writeln!(f, "ORIGBC2\tCORRBC2").unwrap();
        let map =
            crate::commands::dropest_bc_correct::load_mapping_tsv_pub(f.path().to_str().unwrap())
                .unwrap();
        assert_eq!(map.get("ORIGBC1").map(|s| s.as_str()), Some("CORRBC1"));
        assert_eq!(map.get("ORIGBC2").map(|s| s.as_str()), Some("CORRBC2"));
    }

    #[test]
    fn load_mapping_tsv_skips_blank_lines() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "A\tB").unwrap();
        writeln!(f, "").unwrap();
        writeln!(f, "C\tD").unwrap();
        let map =
            crate::commands::dropest_bc_correct::load_mapping_tsv_pub(f.path().to_str().unwrap())
                .unwrap();
        assert_eq!(map.len(), 2);
    }

    #[test]
    fn load_mapping_tsv_missing_file_returns_error() {
        let result =
            crate::commands::dropest_bc_correct::load_mapping_tsv_pub("/nonexistent/path.tsv");
        assert!(result.is_err());
    }

    // ── --cb-tag / --ub-tag CLI parsing ───────────────────────────────────────

    /// Thin wrapper so RunArgs (which derives Args, not Parser) can be parsed
    /// from a command-line slice in tests.
    #[derive(clap::Parser)]
    struct RunWrapper {
        #[command(flatten)]
        inner: RunArgs,
    }

    fn parse_run_args(argv: &[&str]) -> RunArgs {
        use clap::Parser;
        RunWrapper::try_parse_from(argv).unwrap().inner
    }

    #[test]
    fn cb_ub_tag_flags_parse_explicit_values() {
        let args = parse_run_args(&[
            "cmd",
            "sample.bam",
            "ref.gtf",
            "--cb-tag",
            "CR",
            "--ub-tag",
            "UR",
        ]);
        assert_eq!(args.cb_tag.as_deref(), Some("CR"));
        assert_eq!(args.ub_tag.as_deref(), Some("UR"));
    }

    #[test]
    fn cb_ub_tag_flags_default_to_none() {
        let args = parse_run_args(&["cmd", "sample.bam", "ref.gtf"]);
        assert!(args.cb_tag.is_none());
        assert!(args.ub_tag.is_none());
    }

    #[test]
    fn cb_tag_only_leaves_ub_tag_none() {
        let args = parse_run_args(&["cmd", "sample.bam", "ref.gtf", "--cb-tag", "GE"]);
        assert_eq!(args.cb_tag.as_deref(), Some("GE"));
        assert!(args.ub_tag.is_none());
    }

    // ── peek-skip condition ───────────────────────────────────────────────────

    fn should_skip_peek(cb_tag: Option<&str>, ub_tag: Option<&str>) -> bool {
        cb_tag.is_some() && ub_tag.is_some()
    }

    #[test]
    fn peek_skipped_when_both_tags_provided() {
        assert!(should_skip_peek(Some("CR"), Some("UR")));
    }

    #[test]
    fn peek_runs_when_only_cb_tag_provided() {
        assert!(!should_skip_peek(Some("CR"), None));
    }

    #[test]
    fn peek_runs_when_only_ub_tag_provided() {
        assert!(!should_skip_peek(None, Some("UR")));
    }

    #[test]
    fn peek_runs_when_no_tags_provided() {
        assert!(!should_skip_peek(None, None));
    }

    // ── loom dtype validation ─────────────────────────────────────────────────

    #[test]
    #[cfg(feature = "bam")]
    fn validate_loom_dtype_accepts_uint32() {
        assert!(validate_loom_dtype("uint32").is_ok());
    }

    #[test]
    #[cfg(feature = "bam")]
    fn validate_loom_dtype_accepts_uint16() {
        assert!(validate_loom_dtype("uint16").is_ok());
    }

    #[test]
    #[cfg(feature = "bam")]
    fn validate_loom_dtype_rejects_unknown() {
        // Typos and unsupported widths must be rejected, not silently coerced.
        for bad in ["uint8", "unit32", "float32", "u32", ""] {
            assert!(
                validate_loom_dtype(bad).is_err(),
                "expected '{bad}' to be rejected"
            );
        }
    }

    // ── sample-tag demultiplexing: composite cell-key → CellID / SampleID ──────
    //
    // Mirrors the inline logic in `run_inner` that splits each "{sample}|{bc}"
    // counting key into a CellID and a SampleID. Kept in sync by hand, in the
    // same test-only style as `parse_gem_grp` / `derive_sampleid` above.

    fn derive_cellids(
        all_bcs: &[&str],
        sampleid: &str,
        gem_grp: &str,
    ) -> (Vec<String>, Vec<String>, bool) {
        let mut cell_ids = Vec::new();
        let mut sample_ids = Vec::new();
        let mut any_sample = false;
        for key in all_bcs {
            let (sample, bc) = key.split_once('|').unwrap_or(("", key));
            if sample.is_empty() {
                cell_ids.push(format!("{sampleid}:{bc}{gem_grp}"));
            } else {
                any_sample = true;
                cell_ids.push(format!("{sampleid}_{sample}:{bc}{gem_grp}"));
            }
            sample_ids.push(sample.to_string());
        }
        (cell_ids, sample_ids, any_sample)
    }

    #[test]
    fn cellid_single_sample_keeps_legacy_format() {
        // sample empty → unchanged CellID, no SampleID attr written.
        let (cell_ids, sample_ids, any_sample) =
            derive_cellids(&["|ACGT", "|TTTT"], "mysample", "-1");
        assert_eq!(cell_ids, vec!["mysample:ACGT-1", "mysample:TTTT-1"]);
        assert_eq!(sample_ids, vec!["", ""]);
        assert!(!any_sample);
    }

    #[test]
    fn cellid_multi_sample_folds_sample_into_id() {
        let (cell_ids, sample_ids, any_sample) =
            derive_cellids(&["S1|ACGT", "S2|ACGT"], "mysample", "");
        // Same barcode in two samples must yield distinct CellIDs (no collision).
        assert_eq!(cell_ids, vec!["mysample_S1:ACGT", "mysample_S2:ACGT"]);
        assert_ne!(cell_ids[0], cell_ids[1]);
        assert_eq!(sample_ids, vec!["S1", "S2"]);
        assert!(any_sample);
    }

    #[test]
    fn cellid_multi_sample_preserves_gem_group_suffix() {
        let (cell_ids, _, _) = derive_cellids(&["S1|ACGT"], "mysample", "-1");
        assert_eq!(cell_ids, vec!["mysample_S1:ACGT-1"]);
    }
}
