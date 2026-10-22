# velocyto-rs

Faithful Rust translation of [velocyto.py](https://github.com/velocyto-team/velocyto.py) v0.17.16 — RNA velocity analysis for single-cell RNA-seq data.

**Reference snapshot**: velocyto.py v0.17.16 (no git hash available — copied source snapshot).
**Translation approach**: bottom-up, one Rust function per Python function, source files mirror Python module names.

## Why
BD Rhapsody single-cell RNA-seq data uses a non-standard, complex cell barcode scheme, which makes it difficult to match cell barcodes between STARsolo output and the BD Rhapsody™ Sequence Analysis Pipeline output. [velocyto.py](https://github.com/velocyto-team/velocyto.py) remains the most straightforward approach for generating spliced and unspliced count matrices required for RNA velocity analysis; however, the project has been unmaintained since June 2019. This project provides a modern, actively maintained, and high-performance Rust reimplementation capable of processing BD Rhapsody BAM files to produce the necessary count matrices.

> 🙏 **Heartfelt thanks to Prof. Johan Henriksson ([@mahogny](https://github.com/mahogny))** for the original software translation work that made this project possible. Your contribution and inspiration are deeply appreciated!

## New Feature

* output `anndata` or `loom`
* add `--cb-tag`, `--ub-tag`, `--sample-tag` to to suport BAM not from 10x (e.g. BD Rhaposody)



## Performance

Measured with `/usr/bin/time -v` on the same machine, same input files, same barcodes file. Both tools produce identical count output (see Output correctness below). Python: velocyto.py v0.17.16, CPython 3.9. Rust: `cargo build --release --features bam`.

> **Backend note:** the Rust numbers in the table below were measured with the original `rust-htslib` (C `htslib`) BAM backend. BAM/SAM reading has since moved to the pure-Rust [`noodles`](https://github.com/zaeleus/noodles) crates — see [BAM reader: noodles vs htslib](#bam-reader-noodles-vs-htslib), which re-baselines both backends on one machine. The Rust-vs-Python rows here are from a separate run and are indicative, not directly comparable to that same-machine table.

| Dataset | Tool | Wall time | Peak RSS | Speed ratio | RAM ratio |
|---|---|---|---|---|---|
| `mini_chr21.bam`<br>(39 MB BAM, 16 MB GTF) | Python | 22.1 s | 390 MB | — | — |
| | Rust | **2.2 s** | **36 MB** | **10× faster** | **11× less** |
| `s04.bam`<br>(1.7 GB BAM, 1.5 GB GTF) | Python | 14 min 47 s | 3.80 GB | — | — |
| | Rust | **2 min 17 s** | **3.26 GB** | **6.5× faster** | **1.2× less** |

Hardware: Linux x86-64, single machine.

The smaller RAM advantage on the large dataset is expected: the 1.5 GB GTF dominates RSS in both tools once read into memory. The speed advantage is consistent because Rust avoids Python's per-object allocation overhead and GIL contention during GTF parsing and BAM scanning.

### BAM reader: noodles vs htslib

BAM/SAM reading was migrated from `rust-htslib` (which links the C `htslib` library and needs `libclang` to build) to the pure-Rust [`noodles`](https://github.com/zaeleus/noodles) crates, removing the BAM reader's C build dependency. The new streaming output backend separately requires HDF5. Re-baselined on one machine, same inputs, `samtools sort` skipped equally, both built `--release --features bam`; output is bit-for-bit identical.

| Dataset | Metric | `rust-htslib` | `noodles` | Improvement |
|---|---|---|---|---|
| `mini_chr21`<br>(39 MB BAM, 16 MB GTF) | Wall time | 1.26 s | **0.76 s** | **1.7× faster** |
| | Peak RSS | 24.6 MB | **18.5 MB** | **1.3× less** |
| `s04`<br>(1.7 GB BAM, 1.5 GB GTF) | Wall time | 2 m 16 s | **1 m 54 s** | **1.2× faster** |
| | Peak RSS | 3.75 GB | **3.68 GB** | ~2% less |
| Release binary | Size | 93 MB | **59 MB** | **1.6× smaller** |

The small-input win is larger because noodles decodes records lazily (only the fields you touch) with no C FFI boundary; on the large input the unchanged 1.5 GB GTF load dominates wall time and RSS, diluting the reader speedup. The smaller binary reflects no statically-linked C `htslib`.

### Output correctness

Comparing Rust vs Python output on the same input (S08, 1974 cells × 39579 genes):

| Dataset dimension | Agreement |
|---|---|
| Gene set | Identical |
| Cell barcode set | Identical |
| Gene row order | Identical |
| Layer counts (spliced/unspliced/ambiguous) | **100% identical** — 0 of 78 M (gene, cell) pairs differ |
| Layer dtype | Rust uses `uint32` same as Python |

**Comparison note:** 69 gene names appear twice in the dataset (same name, different accession/gene ID). A naive gene-name-keyed comparison will silently cross-map these and produce false diffs. Always align by `row_attrs/Accession` (unique gene ID), not by `row_attrs/Gene` (gene name).

### Where the speedup comes from

The dominant cost for large reference genomes is loading the GTF annotation. Three optimisations in `read_transcriptmodels`:

1. **Schwartzian transform for GTF sorting** — Python's `list.sort(key=f)` computes the sort key once per element (built-in). The original Rust `sort_by(|a,b| f(a).cmp(f(b)))` recomputed the key — including a full line clone — on every comparison, O(n log n) times. Replaced with explicit key precomputation (O(n) extractions) then index sort.

2. **Filter to exon-only lines during collection** — The annotation builder only uses `feature_type == "exon"` lines. Gencode v42 is ~⅓ exon lines, so filtering at read time reduces both allocation count and sort size ~3×.

3. **1 MB read buffer** — reduces syscall overhead vs the default 8 KB `BufReader` for large files.

## Build

```bash
cargo build --release
```

BAM support is **on by default**, using pure-Rust [`noodles`](https://github.com/zaeleus/noodles). Streaming output requires the native HDF5 development library and `pkg-config` (tested with HDF5 1.14.6), plus Rust 1.89 or later. On Debian/Ubuntu install `libhdf5-dev pkg-config`; on Slurm load your site's HDF5 module before building and running. `--no-default-features` disables BAM support; checkpoint export still requires HDF5.

> **No samtools required at runtime.** Quantification uses the built-in Rust disk-spilling sort, sequentially across input BAMs, with a maximum merge fan-in of 32. Temporary sort runs are removed after a successful sort. On restart, a completed staged BAM is published only after its saved identity checks pass; incomplete sort state is removed and sorting starts again from the source BAM. Existing `--samtools-threads` and `--samtools-memory` names remain accepted and control Rust sort threads and the total chunk budget in MiB. A sorted BAM is reused only when its completion receipt matches the source, output and barcode tag.

## Usage

BAM support is on by default — no extra feature flag is needed. Run with `--help` for the full option list.

### `run10x` — 10X Chromium (CellRanger output)

The most common entry point. Pass the CellRanger sample folder and a GTF file; the BAM and barcodes list are found automatically.

```bash
velocyto-rs run10x /data/cellranger/sample1 /ref/gencode.v42.gtf
```

With a repeat-mask GTF and verbose logging:

```bash
velocyto-rs run10x /data/cellranger/sample1 /ref/gencode.v42.gtf \
    --mask /ref/repeats.gtf \
    --verbose
```

Output is written to `<samplefolder>/velocyto/<samplename>.h5ad` by default (AnnData).
Pass `--output-format loom` for the legacy loom file, or `--output-format both` for both.

### `run` — generic BAM

Use this for any platform not covered by a dedicated subcommand (BD Rhapsody, Parse Biosciences, etc.).

```bash
velocyto-rs run sample.bam /ref/gencode.v42.gtf \
    --bcfile barcodes.txt \
    --outputfolder ./results \
    --sampleid my_sample
```

Multiple BAMs merged into one loom:

```bash
velocyto-rs run lane1.bam lane2.bam lane3.bam /ref/gencode.v42.gtf \
    --bcfile barcodes.txt \
    --sampleid pooled_run
```

Key options:

| Flag | Default | Description |
|---|---|---|
| `-b / --bcfile` | — | Barcode whitelist (plain text or `.gz`, one per line) |
| `-o / --outputfolder` | `<bam-dir>/velocyto` | Output directory |
| `-e / --sampleid` | derived from BAM name | Loom filename stem |
| `-m / --mask` | — | GTF of genomic intervals to mask (e.g. repeats) |
| `-l / --logic` | `Default` | Molecule-filtering logic class |
| `-U / --without-umi` | false | Read-count mode (no UMI deduplication) |
| `-u / --umi-extension` | `no` | Extend UMI identity: `no`, `chr`, `Gene`, `Cluster`, `all` |
| `-M / --multimap` | false | Count non-uniquely mapped reads (not recommended) |
| `-t / --dtype` | `uint32` | Layer array dtype; use `uint16` for low-depth data |
| `--output-format` | `h5ad` | Output file format: `h5ad` (AnnData), `loom`, or `both` (see [Output formats](#output-formats)) |
| `--samtools-threads` | 16 | Built-in Rust sort threads, bounded by CPU affinity and Slurm allocation |
| `--samtools-memory` | 2048 | MiB per in-memory Rust sort chunk; bounded by available/cgroup/Slurm memory |
| `--cb-tag` | auto | BAM tag for cell barcode (e.g. `CB`, `XC`); skips auto-detection when both tags are set |
| `--ub-tag` | auto | BAM tag for UMI barcode (e.g. `UB`, `XM`); skips auto-detection when both tags are set |
| `--sample-tag` | — | BAM tag carrying sample identity (e.g. BD Rhapsody `ST`); demultiplexes a multi-sample BAM in place (see [Sample demultiplexing](#sample-demultiplexing)) |
| `--cell-batch-size` | 100 | Maximum raw barcode groups in a count batch; a barcode group is never split |
| `--count-memory-mb` | 512 | Budget for count read, evidence, and scratch buffers; not a total RSS limit |
| `--checkpoint-dir` | `<output>/<sample>.counts` | Directory for sparse, resumable count checkpoints |
| `--resume-from` | — | Resume a checkpoint after validating inputs and settings |
| `--telemetry-jsonl` | — | Append structured phase, progress, RSS, and throughput events |

## Library usage

Besides the CLI, `velocyto-rs` can be used as a Rust library. Every CLI subcommand is a
public function that takes an `Args` struct and runs the full pipeline — counting, then
writing the output file(s) to disk — returning `anyhow::Result<()>`.

> **Crate vs library name.** The package is `velocyto-rs`, but the library is imported as
> **`velocyto`** (`[lib] name = "velocyto"`). Depend on `velocyto-rs`, write `use velocyto::...`.

Add the dependency (not yet on crates.io — use a git dependency):

```toml
[dependencies]
velocyto-rs = { git = "https://github.com/Thokas99/velocyto-rs-lite", branch = "feature/sparse-counting" }
anyhow = "1"
# BAM support is enabled by default. For the dependency-light stub build
# (no BAM, commands error at runtime) use:
# velocyto-rs = { git = "https://github.com/Thokas99/velocyto-rs-lite", branch = "feature/sparse-counting", default-features = false }
```

### `run10x` — 10X Chromium

```rust
use velocyto::commands::run10x::{run10x, Run10xArgs};

fn main() -> anyhow::Result<()> {
    run10x(Run10xArgs {
        verbose: false,
        samplefolder: "/data/cellranger/sample1".into(),
        gtffile: "/ref/gencode.v42.gtf".into(),
        metadatatable: None,
        mask: None,                      // or Some("/ref/repeats.gtf".into())
        logic: "Default".into(),
        multimap: false,
        samtools_threads: 16,
        samtools_memory: 2048,
        dtype: "uint32".into(),
        output_format: "h5ad".into(),    // "h5ad" | "loom" | "both"
        dump: "0".into(),
        cb_tag: None,
        ub_tag: None,
    })?;
    Ok(())
}
```

This writes `/data/cellranger/sample1/velocyto/sample1.h5ad` as a side effect, exactly as the
CLI does.

### `run` — generic BAM (with sample demultiplexing)

The generic entry point takes one or more BAMs and the full option set, including `sample_tag`
for in-place demultiplexing (see [Sample demultiplexing](#sample-demultiplexing)):

```rust
use velocyto::commands::run::{run, CountingOptions, RunArgs};

fn main() -> anyhow::Result<()> {
    run(RunArgs {
        counting: CountingOptions::default(),
        verbose: false,
        bamfile: vec!["multiplexed.bam".into()],   // one or more BAMs, merged into one output
        gtffile: "/ref/gencode.v42.gtf".into(),
        bcfile: Some("barcodes.txt".into()),
        outputfolder: Some("./results".into()),
        sampleid: Some("pooled_run".into()),
        metadatatable: None,
        mask: None,
        onefilepercell: false,
        logic: "Default".into(),
        without_umi: false,
        umi_extension: "no".into(),                // "no" | "chr" | "Gene" | "Cluster" | "all"
        multimap: false,
        samtools_threads: 16,
        samtools_memory: 2048,
        dtype: "uint32".into(),
        output_format: "h5ad".into(),
        dump: "0".into(),
        cb_tag: None,                              // e.g. Some("XC".into()) for non-10x BAMs
        ub_tag: None,                              // e.g. Some("XM".into())
        sample_tag: Some("ST".into()),             // None for single-sample BAMs
    })?;
    Ok(())
}
```

The remaining subcommands follow the same shape:

| Function | Args struct | Module |
|---|---|---|
| `run10x` | `Run10xArgs` | `velocyto::commands::run10x` |
| `run` | `RunArgs` | `velocyto::commands::run` |
| `run_dropest` | `RunDropestArgs` | `velocyto::commands::run_dropest` |
| `run_smartseq2` | `RunSmartseq2Args` | `velocyto::commands::run_smartseq2` |
| `dropest_bc_correct` | `DropestBcCorrectArgs` | `velocyto::commands::dropest_bc_correct` |

> **Logging.** These functions do not configure a logger — the CLI sets one up in `main.rs`. To
> see progress output (the `verbose` field only controls log *level*), initialize a logger in
> your own code, e.g. `env_logger::builder().filter_level(log::LevelFilter::Info).init();`.

## Output formats

`--output-format` is available on all four run commands (`run10x`, `run`, `run_dropest`, `run_smartseq2`) and accepts:

| Value | Writes | Notes |
|---|---|---|
| `h5ad` *(default)* | `<sampleid>.h5ad` | AnnData, the modern scanpy/scVelo on-disk format |
| `loom` | `<sampleid>.loom` | Legacy loompy format (what upstream velocyto.py emits) |
| `both` | `<sampleid>.h5ad` and `<sampleid>.loom` | Both files from a single counting pass |

> **Note:** AnnData (h5ad) output is an extension beyond the v0.17.16 port — upstream velocyto only ever emitted loom. The two files carry identical counts; they differ only in on-disk layout.

CLI output streams from persistent sparse count checkpoints through the maintained `hdf5-metno` Rust binding to native HDF5. No Python is needed at runtime. Legacy library writers remain available for compatibility and comparison; they still collect dense arrays and are unsuitable for large exports.

### h5ad layout

The h5ad file targets the current anndata on-disk spec and round-trips through `anndata.read_h5ad`:

- **Cells × genes** (`obs` × `var`), the transpose of loom's genes × cells.
- **Sparse CSR** for `X` and every `layers/{spliced,unspliced,ambiguous}` — a single-cell matrix is mostly zeros, so CSR is far smaller than dense. `X` is the float32 sum of all layers; layer width follows `--dtype` (`uint32` lossless / `uint16` saturating).
- **`var` indexed by `Accession`** (unique Ensembl ID), with `Gene` as a column. Gene *names* collide (see the alignment note above), so the accession must be the index.
- **`obs` indexed by `CellID`**; `SampleID` appears as an `obs` column only when `--sample-tag` demultiplexing is active.

### Loom layout

Matches loompy's layout exactly (genes × cells): a float32 `matrix` (sum of all layers), `layers/{spliced,unspliced,ambiguous}` (`uint32` by default), and `row_attrs` / `col_attrs` string datasets including `Gene`, `Accession`, `Chromosome`, `Strand`, and `CellID`. This is byte-for-byte comparable to upstream velocyto.py output via `h5diff`.

## Sample demultiplexing

`--sample-tag` (generic `run` command only) lets a single multi-sample BAM be counted without `samtools split`. This is an extension beyond the v0.17.16 port — upstream velocyto expects one sample per BAM.

```bash
velocyto-rs run multiplexed.bam /ref/gencode.v42.gtf \
    --bcfile barcodes.txt \
    --sample-tag ST \
    --sampleid pooled_run
```

- Cell identity becomes `(sample, barcode)`, so the same bead barcode in two samples no longer collides. CellIDs are formatted `{sampleid}_{sample}:{bc}` and a `SampleID` column (loom `col_attrs/SampleID`, h5ad `obs` column) is added.
- The BAM is still sorted by `CB` only — no second sort, no split.
- Reads missing the sample tag are dropped (not pooled into a phantom sample). Reads that carry the tag with any value (e.g. BD `"Multiplet"`) are kept as their own sample — filter them downstream via `SampleID`.
- With `--sample-tag` off, output is byte-for-byte identical to before.


## Citation

If you use this software, cite the original velocyto paper:

> La Manno G, Soldatov R, Zeisel A, Braun E, Hochgerner H, Petukhov V, Lidschreiber K, Kastriti ME, Lönnerberg P, Furlan A, Fan J, Borm LE, Liu Z, van Bruggen D, Guo J, He X, Barker R, Sundström E, Castelo-Branco G, Cramer P, Adameyko I, Linnarsson S, Kharchenko PV.
> **RNA velocity of single cells.**
> *Nature* 560, 494–498 (2018). https://doi.org/10.1038/s41586-018-0414-6


## License

BSD 2-Clause License. See [LICENSE](LICENSE). Derived from [velocyto.py](https://github.com/velocyto-team/velocyto.py) by Peter Kharchenko, Sten Linnarsson and Gioele La Manno.

## Recoverable quantification

For BD Rhapsody, retain the explicit tag overrides:

```bash
velocyto-rs run /data/sample.bam /data/genes.gtf \
  -b /data/barcodes.txt -o /persistent/velocity -e sample \
  --cb-tag CB --ub-tag MA \
  --checkpoint-dir /persistent/velocity/sample.counts \
  --cell-batch-size 100 --count-memory-mb 512 \
  --telemetry-jsonl /persistent/velocity/sample.jsonl
```

Each completed raw-barcode group batch is compressed into sparse integer rows, synced to disk, and atomically added to a versioned manifest. Counts do not accumulate in a whole-dataset dense matrix. Small dense batch scratch preserves the existing molecule-assignment logic. Spliced, unspliced, ambiguous and any additional logic layers retain uint32 counts. The existing `uint16` export option remains explicitly lossy; checkpoints retain uint32 values.

`--telemetry-jsonl` appends phase completion, process RSS and peak RSS, sampled record throughput, active batch estimates, checkpoint commits, source completion and errors. Progress is sampled every 100,000 BAM records and after completed batches; each JSONL event is flushed. RSS describes the velocyto-rs process and does not include child processes or aggregate Slurm cgroup memory.

After a failed export, regenerate output **without BAM, GTF or barcode files**:

```bash
velocyto-rs write-output /persistent/velocity/sample.counts \
  -o /persistent/velocity --output-format h5ad
```

To resume interrupted counting, repeat the original `run` arguments and replace `--checkpoint-dir` with `--resume-from /persistent/velocity/sample.counts`. Completed batches are verified and counting seeks to the next unread BGZF record. Counting resume currently rebuilds annotation and reruns intron validation on the original BAMs; output-only restart does neither. Input identities use canonical paths, byte sizes and modification times, not cryptographic content hashes: inputs must remain immutable. Only generic `run` exposes count-budget/resume flags; all run wrappers use durable checkpoints with default budgets, and `write-output` works for their checkpoints too.

Use persistent project storage for count checkpoints, not job-local scratch that Slurm may remove. Successful sorts remove their spill files. An abrupt interruption can leave a sort work directory; the next run validates it, reuses a complete staged BAM, or deletes incomplete state and sorts again. Unmarked legacy `cellsorted_*` BAMs are refused rather than trusted; move them aside before a fresh sort. Checkpoint durability depends on filesystem locking, atomic rename, and fsync support.

`--count-memory-mb` limits estimated read, mapping-evidence and scratch components separately; it is not a total-RSS guarantee. Annotation/index structures, metadata and HDF5 caches add memory. Oversized individual barcode/UMI groups fail explicitly with checkpoints intact; disk-backed processing of a single oversized molecule is not implemented. Reduce the batch size or increase the component budget within the job allocation. Cell groups are never split to force a fit. H5AD uses int64 row offsets; `X` remains float32 and can round large integer totals, while integer layers remain exact.

See [audit/design](docs/memory-efficient-quantification-plan.md), [implementation and verification record](docs/quantification-implementation.md), and [benchmark results](docs/quantification-benchmarks.md).

### Measuring a Slurm run

The benchmark driver records wall time, peak process RSS, CPU time, logical `/proc` I/O, storage block I/O, exit status and sampled intermediate disk size. It writes a report even when the child is killed. Run it inside the allocation:

```bash
python3 scripts/benchmark_quantification.py \
  --report /persistent/velocity/run-metrics.json \
  --workdir /persistent/velocity -- \
  velocyto-rs run /data/sample.bam /data/genes.gtf \
    -b /data/barcodes.txt -o /persistent/velocity -e sample \
    --cb-tag CB --ub-tag MA --cell-batch-size 100 --count-memory-mb 512
```

Also retain Slurm `sacct` and cgroup peak-memory accounting: per-process RSS is not whole-job memory. Cached/tmpfs I/O can report zero storage block traffic; logical writes and final disk occupancy are separate measurements. If Slurm kills the entire job including the benchmark parent, its report cannot be finalized, but committed counting checkpoints survive.
