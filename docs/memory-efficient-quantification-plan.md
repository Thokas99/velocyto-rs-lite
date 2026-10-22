# Memory-efficient quantification: audit and proposed implementation

Date: 2026-10-20. Status: initial audit and prioritized design. Implementation and current verification are summarized in [quantification-implementation.md](quantification-implementation.md); measured synthetic storage results are in [quantification-benchmarks.md](quantification-benchmarks.md).

## Goal and evidence limits

Quantify approximately 187,000 BD Rhapsody WTA cells from approximately two billion alignment records without retaining dataset-sized dense count arrays. Preserve molecule classification, cell selection, CB/MA extraction, sample separation, and exact integer layer counts. Persist completed work and support output-only restart without reading BAMs or GTFs.

The reported Slurm termination near 150 GiB RSS is user-supplied evidence. No job log, BAM, GTF, cell whitelist, or phase-resolved memory trace was supplied. The allocations below are verified in source; attribution of the observed peak to any particular allocation remains a hypothesis. Later synthetic storage measurements are documented separately; no production performance improvement is claimed.

Repository: `velocyto-rs-lite`, branch `feature/sparse-counting`, baseline HEAD `17df61dae812f132c5b955b55b6e6e3efc12a1b7`. At the initial audit, GitNexus reported the indexed baseline commit, 804 symbols and 70 processes; that graph predates implementation. The current `detect_changes --scope all` report identifies 71 changed symbols, 34 affected execution flows, and critical aggregate risk. Sort, count-stream, and annotation impacts have regression and differential coverage. The new untracked telemetry module is absent from that graph; local re-indexing was blocked because the MCP GitNexus process held the database. Existing untracked `AGENTS.md` and `repomix-output.xml` were left untouched.

At audit time, baseline verification was blocked because `argmin` was not cached and `index.crates.io` was unreachable. The baseline was subsequently built and tested after dependencies became available; see the implementation record. Inspection of `hdf5-pure-rust 0.3.6` showed whole-chunk buffering, so the implementation uses a native HDF5 backend for slice writes.

## Observed architecture

| Location | Observed behavior | Consequence |
| --- | --- | --- |
| `src/commands/run.rs::run_inner`, counting section | Receives all dense batches from `count`; extends `all_layers` with `batches.iter().cloned()` | During transfer, original and cloned batches coexist. Moving ownership removes this duplicate, but does not eliminate accumulated dense storage. |
| `src/counter.rs::ExInCounter::count` (1142) | Retains every completed dense batch in `dict_list_arrays` until a BAM finishes | RAM grows with retained cells × annotated genes × layers. |
| `src/counter.rs::finalize_batch` (1642) | Allocates a zeroed genes × cells `Array2<u32>` for every logic layer | A small batch is acceptable as a compatibility bridge; storing all batches is not. |
| `count`, unfiltered mode | Computes spliced + unspliced totals, keeps cells with strictly >80 counts, and stacks kept columns into fresh arrays | Temporary dense copies; preserve this selection rule exactly, including exclusion of ambiguous counts. |
| `src/counter.rs::dump_anndata` (1898) | Concatenates all batches into full dense layers, builds a full dense float32 total, then constructs complete CSR vectors | Sparse output does not imply sparse or bounded-memory construction. `mapv` also creates temporary float arrays. |
| `src/counter.rs::dump_loom` (1700) | Concatenates all layers, constructs a dense float total, and flattens arrays into new vectors | Selecting `loom` or `both` retains the same scalability problem. |
| H5AD CSR conversion | Casts `data.len()` to `i32` without checking | Row offsets can wrap above 2,147,483,647 nonzeros. Alignment count alone does not determine nonzero count. |
| Both output writers | Delete an existing output before creating its replacement | A failed retry can destroy the previous valid output and leave a partial file at the final path. |
| `count_cell_batch_*_inner` (1474–1640) | Retains reads, sorts by genomic coordinates, and builds a `HashMap<String, Molitem>` | Peak is also driven by reads and mapping evidence within a barcode batch. |
| `src/molitem.rs::dictionary_intersect` | Concatenates `SegmentMatch` vectors for surviving transcript keys | Even a single UMI can retain arbitrarily many mapping records. Do not deduplicate evidence without semantic proof. |
| `count` batch boundary | Flushes only after 100 complete raw barcodes, across all interleaved samples | Preserves UMI evidence; does not bound read count, bytes, or composite-cell count. |
| `build_feature_indexes` (around 885) | Keeps annotation models, cloned `tms_flat`, and cloned indexed features | Annotation-sized duplication, lower priority than dataset-sized matrices; index identity and validation state make casual deduplication risky. |
| `src/bam_sort.rs::sort_bam_by_tag` | Spill directory is a `TempDir`; errors unwind through automatic cleanup; merge opens all runs | Recoverable spills can disappear on errors; many runs increase file handles and buffers without a fixed fan-in. SIGKILL bypasses Rust destructors, but that is not a retention guarantee for all failures. |
| `run_inner`, sort setup | Accepts existing `cellsorted_*` files solely by existence and writes sorts directly to those names | A partial sort can be reused after failure. |
| `read_mem_available_mb`, `num_cpus` | Reads host `/proc` resources; starts sorts for multiple BAMs concurrently | Host resources can exceed Slurm/cgroup allocation. External samtools processes also contribute to job memory. |
| `peek` (262) | Autodetects CB/UB or XC/XM; explicit tag overrides exist | BD CB/MA is already supported with `--cb-tag CB --ub-tag MA`; preserve this path. Autodetecting MA is separate scope. |

For C cells, G genes and L layers, dense count storage alone costs `4*C*G*L` bytes. At 187,000 cells, 30,000 genes and three layers, this is **62.7 GiB**. A second concatenated set costs another 62.7 GiB; a float32 total costs 20.9 GiB. These are illustrative dimensions, not a measured peak: allocation lifetimes, expression temporaries, CSR vectors, allocator behavior and HDF5 buffers affect RSS. Actual annotated gene count is needed for a dataset-specific estimate.

## Proposed architecture and alternatives

Recommended flow: existing read/molecule assignment → bounded batch counts → per-batch sparse CSR checkpoint → release batch → stream checkpoint data into H5AD and, when requested, tiled loom output.

1. **Recommended: retain dense scratch only within a bounded batch.** Preserve `Logic::count` and all classification branches initially. Convert immediately to cells × genes CSR, discard zeros, persist and drop. At 30,000 genes and 100 composite cells, three scratch layers cost about 34 MiB. Composite cells must be bounded separately from raw barcode count; this estimate excludes reads and mapping evidence.
2. **Alternative: change every logic implementation to write sparse counters.** Removes batch zero initialization, but changes the trait and seven implementations, expands the regression surface, and does not solve large read groups or output buffering. Defer until profiling shows scratch matrices matter.
3. **Insufficient alone: collect all sparse matrices in memory.** Improves typical memory, but grows with total nonzeros and remains unbounded with dataset size. Do not choose this as the final architecture.

Disk-backed sparse checkpoints are the durable source of truth. Neither a monolithic in-memory sparse matrix nor a memory-mapped file, by itself, proves bounded job memory.

## Prioritized implementation and acceptance checks

Each implementation step was tested before advancing. Before any future code edits, run GitNexus impact on the edited symbols. The current graph report is not a clean gate: aggregate risk is critical and the new telemetry file is not indexed. No commit was made. No unrelated analysis-module rewrite is needed.

### P0. Establish correctness and measurement baselines

Reach: existing `counter.rs` tests, new quantification fixtures/tests, a benchmark driver, and CI where necessary. Freeze the original HEAD as a separately built reference binary; record compiler, lockfile, logic, tags, flags, input identities and gene/cell ordering.

Use small synthetic BAM/GTF fixtures with CB/MA tags, duplicate UMIs, conflicting transcript evidence, validated and unvalidated introns, spliced/unspliced/ambiguous assignments, repeat masking, multimapping, absent tags, sample-tag interleaving, whitelist/no-whitelist selection, and multiple BAMs. Cover all logic variants, including the optional `spanning` layer. Include >65,535 counts, empty layers, empty input, and multiple output batches.

Compare every integer layer element after alignment by exact CellID and gene accession, not hash iteration order or totals alone. Preserve duplicate CellIDs across BAMs if that is the reference behavior; use source/occurrence identity for comparison and never silently merge them. Store the baseline outputs as independent fixtures so tests do not simply mirror the optimized implementation. A same-library HDF5 round-trip is useful but insufficient: independently open H5AD with Python AnnData in integration verification.

Baseline benchmark cases: increasing cell count at fixed genes/read depth, increasing read depth for a single barcode, increasing duplicate depth for one UMI, increasing genes/sparsity, and a representative real BD subset. Repeated runs use identical data and resources; record cache state. Do not extrapolate a two-billion-record success from a tiny fixture.

### P1. Remove counting clones

Reach: `src/commands/run.rs::run_inner`; its callers `run`, `run10x`, `run_dropest`, `run_smartseq2` keep their existing signatures.

Consume returned maps and move owned arrays into `all_layers` instead of cloning. Preserve layer and barcode order. Test a multi-BAM, multi-layer fixture against the baseline; confirm transferred arrays retain their allocation where practical. Expected effect: removes the extra per-BAM transfer copy. This is a useful isolated fix, not a complete memory solution.

### P2. Add sparse, atomic batch checkpoints and a batch sink

Reach: `counter.rs::count`, its batch-flush block and callers; a new small `count_checkpoint` module; `lib.rs`; metadata extraction currently embedded in the output writers. Retain the public collecting `count` API as a compatibility adapter if needed, while the CLI uses the sink path.

The sink receives one finalized batch after existing filtering. Convert each layer to cells × genes CSR with sorted gene indices, `u32` values, and checked 64-bit row offsets/nonzero totals. Record all cell/sample identities and row order. Emit even all-zero cell rows. Do not change assignment branches in `logic.rs` or evidence aggregation in `molitem.rs`.

Use a versioned, explicitly encoded binary format with bounded length validation, dimensions, a checksum, and no platform-dependent `usize` fields. Standard-library little-endian I/O is sufficient; reuse installed compression/checksum facilities if appropriate. Cap reader allocations against declared schema and configured limits, validate monotonic pointers, index bounds, lengths and trailing data. Read payloads in bounded chunks rather than loading complete checkpoint files.

Commit protocol: create a unique persistent run directory; write a batch `.partial` using exclusive creation; finish compression/checksum; flush and `sync_all`; rename to an immutable complete batch on the same filesystem; sync the directory. Atomically publish the manifest entry only after the batch is durable. A crash between these steps may leave an orphan, which is retained and ignored until explicitly reconciled. A committed manifest must never reference a partial shard.

The manifest includes schema and assignment version, program build identity, ordered input identities, logic/settings, CB/MA/sample tags, whitelist/mask identities, sample naming/GEM suffix, ordered gene metadata, intron-validation state identity, batch ranges, per-layer totals/nonzeros/checksums, and counting-complete status. Persist validated annotation state separately for counting resume. Canonicalize or preserve serialized transcript/feature indices; rebuilding randomized indices and restoring flags by position is unsafe.

Never remove checkpoint files, partial files or manifests automatically after an error or success; cleanup remains explicit. Refuse accidental checkpoint-directory reuse and simultaneous writers. Report recovery paths in logs.

Tests: dense-to-CSR exact equivalence; zero/empty rows; shuffled gene ordering; overflow and malformed lengths; truncated/corrupt shards; failures before/after rename and before manifest publication; no double application of committed batches; all checkpoint artifacts retained on error. Expected effect: count storage RAM scales with one active batch rather than all retained cells; disk scales with nonzeros.

### P3. Stream output and make output-only restart independent

Reach: `counter.rs::dump_anndata`, `dump_loom`, shared metadata extraction; new checkpoint-backed output functions; `src/main.rs` and command exports; README. Preserve existing run flags and defaults. Add an additive `write-output <checkpoint-dir> --output-format h5ad|loom|both` command; it requires a counting-complete manifest and needs no BAM/GTF access. Add `--checkpoint-dir` and `--resume-from` only where needed, preserving existing wrappers.

First gate: inspect and test locked `hdf5-pure-rust 0.3.6` for incremental writes, metadata finalization, per-dataset caching and internal byte copies. Deflate/chunk options on `write(&whole_vector)` do not prove streaming. If supported, precompute dataset lengths from shard metadata, create datasets once, and write bounded slices. Otherwise use an established Rust HDF5 binding with hyperslab writes, accepting/documenting its native HDF5 build requirement. Do not implement a custom HDF5 serializer or silently retain all output buffers. Resolve this gate before claiming P3 feasible with the current dependency.

For each layer, stream data/indices and adjust local row pointers with checked global offsets; use signed 64-bit `indptr` in H5AD and checked gene-index widths. Write metadata without constructing full count arrays. Derive sparse `X` by merging the layer rows within one cell/batch; preserve its documented float32 representation and test rounding separately. Exact preservation is promised for integer layers; float32 cannot represent every integer above 2^24 exactly. Use deterministic layer order and document any difference in high-count float rounding from the old hash-dependent addition.

For loom, write bounded genes × cells tiles and tile totals. `both` must execute both bounded writers sequentially against the same checkpoint source. If the backend cannot stream loom, P3 is incomplete; do not fall back silently to full dense concatenation.

Keep `uint32` as the default lossless layer export. Existing `uint16` is explicitly lossy and saturates; preserve this opt-in compatibility behavior and warning, while checkpoints always retain full counts. Add checked accumulation/overflow handling at the counting boundary before claiming exactness outside the existing u32 range; do not describe u32 as unbounded.

Write to a unique sibling partial output. Flush, close, validate and sync it before atomically publishing the final path. Preserve an existing valid output on failure. Output restart regenerates from the committed checkpoints and retains any previous partial output. Disk-full errors leave the checkpoint source intact.

Tests: compare old/new writer output by identity and layer values across uneven batches; 64-bit offsets without allocating billions of entries; zero cells/genes/nonzeros; all requested output formats; independent AnnData read; fault injection mid-output; output-only success with original BAM/GTF paths inaccessible; rerun produces identical integer values. Expected effect: output RAM bounded by buffer/tile sizes plus metadata and backend caches, rather than C×G or total nonzeros.

### P4. Counting resume and safe sort reuse

Reach: counting reader and flush protocol in `counter.rs`; startup/sort handling in `run_inner`; `bam_sort.rs::sort_bam_by_tag` and merge setup; CLI/checkpoint manifest.

Checkpoint only after a whole raw CB group across all samples is finalized. Persist the next unread record's BGZF virtual position and sorted-input identity. Capture the cursor before reading any lookahead record; do not advance past the first record of the next batch. Validate seeking with the locked noodles versions. If seeking is unavailable, explicitly document replay-and-skip as a slower fallback; it rereads the prefix and must not reassign completed molecules. Do not conflate output-only restart with counting resume.

A last committed batch is the recovery unit. Do not attempt to recover arbitrary in-memory reads lost to SIGKILL. Restore validated annotation state and all filtering settings before continuing. Refuse mismatched inputs/configuration, incomplete validation state, unsupported schema, or corrupt committed shards. Full content hashes are strongest but incur an extra pass at this scale; compute during existing reads where possible, and explicitly document weaker path/size/mtime checks if offered. Never silently assume a replacement file is identical.

Write sorted BAMs to partial names; finalize BGZF, verify integrity/order as appropriate, then publish atomically with a completion manifest bound to inputs, sort tag and sorter settings. Never reuse an unmarked file solely because it exists. Make spill directories persistent before errors can unwind; publish a bounded spill-run manifest. On restart, reuse only validated complete runs; retain incomplete runs for inspection. Limit merge fan-in with multi-pass merging, retaining source runs until output is committed.

Tests: interrupts at group boundaries, lookahead and EOF; CB values shared by multiple samples; filtering of all reads in a group; resume in the second of multiple BAMs; corrupted cursors and mismatched input; failed sort and truncated existing BAM; injected merge failure retaining runs; bounded file handles with many runs. Expected effect: completed groups survive process failure, and failed sorts are neither discarded nor mistaken for successful inputs.

### P5. Bound pathological barcode and UMI evidence

Reach: `count` read batching, three `count_cell_batch_*_inner` paths, `Molitem` evidence handling, and bounded disk spill support. This is required before claiming generally bounded counting memory.

First shrink batches at complete raw-barcode boundaries using measured read-byte and composite-cell limits. For an oversized individual barcode, spill decoded reads to persistent runs and externally reproduce the existing stable genomic comparator (include original ordinal for ties). Preserve all samples of that barcode. Use bounded merge fan-in.

The current Molitem vectors remain an unbounded second source even after read spill. Stage mapping evidence to disk keyed by complete molecule identity and genomic processing ordinal, then process molecule evidence through a bounded representation. Either prove a compact sufficient-state representation equivalent for every logic branch (including evidence multiplicity), or provide a disk-backed evidence iterator with a carefully tested classification adapter. Do not replace concatenated vectors with sets on intuition: singleton/non-singleton decisions can change. This is the highest semantic-risk portion and needs dedicated impact analysis and differential tests before implementation.

Define an explicit behavior if an individual record/evidence object exceeds the configured budget: retain work and fail with an actionable error rather than exceed the limit or split a molecule. A memory cap with a diagnostic failure is preferable to a false guarantee that every adversarial input completes. Acceptance on target BD inputs additionally requires completion, not merely safe failure.

Tests: one deeply sequenced barcode; one UMI with repeated evidence; same UMI spanning chromosomes/strands; empty transcript intersections; repeated equal-coordinate reads; tiny forced budgets; multiple interleaved samples; all logic variants; restart with spill files. Compare exact output to the original for all representable fixtures. Expected effect: memory becomes independent of barcode read depth, subject to verified backend/allocator limits.

### P6. Slurm budgets and performance acceptance

Reach: `read_mem_available_mb`, `num_cpus`, sort orchestration, README/Slurm example, benchmark driver. Honor process affinity/cgroup CPU and memory limits plus explicit CLI budgets; treat Slurm environment values as additional bounds with node/task semantics checked. Reserve headroom for annotations, active count batches, output buffers and compression. Avoid launching a full resource budget per input concurrently. Pure-Rust record-size estimates are not hard RSS limits.

Benchmark original and optimized release builds in fresh processes, with fixed thread counts and dataset/config identities. Record wall time, CPU time, exit status, peak process RSS, phase RSS, logical bytes read/written, physical storage I/O where available, checkpoint/output sizes and peak intermediate disk use. GNU time process RSS does not measure aggregate concurrent job memory: also collect cgroup memory peak and Slurm accounting, including sorter subprocesses. Distinguish logical I/O, physical I/O and page-cache effects.

Use `/usr/bin/time -v` plus phase sampling of `/proc/<pid>/status` and `io`, and Slurm/cgroup accounting on the target cluster. Keep raw measurements and commands. Include counting-only, output-only, full pipeline, forced-crash recovery, and full 187k-cell output tests. Small scaling cases should show flat output RSS as cells/nonzeros grow; counting must also pass the heavy-barcode and heavy-UMI tests. Establish the numerical RSS target from the actual Slurm allocation after measuring annotations/backend overhead; no unsupported fixed GiB promise.

Test after each change, then run full default-feature tests and CLI compatibility checks. Static parsing, synthetic correctness, local benchmarking and a successful full Slurm job are separate evidence levels. Stop optimization when acceptance passes; do not refactor unrelated velocity analysis code.

## Graph impact and review scope

GitNexus upstream impact reports `run_inner`: LOW, direct callers `run`, `run10x`, `run_dropest`, `run_smartseq2`, then `main`. `finalize_batch`: LOW, three batch-counter paths feeding `count`. `dump_anndata`: LOW, its round-trip test appears as the only caller in the graph. Direct source inspection also finds its production call in `run_inner`; the graph is incomplete for that method access despite matching HEAD, so its LOW label is not proof of isolated impact.

`sort_bam_by_tag`: **HIGH risk**, 11 upstream symbols, six direct callers, including `run_inner` and five sort tests at the pre-edit audit. The sort path was subsequently changed and tested with fast-path, forced-spill and bounded-merge cases. Recheck impact on the current graph before further sort edits or any commit.

## Feasibility assessment

High confidence: ownership moves, bounded scratch arrays, sparse checkpoint storage and separation of counting from output address concrete source-level allocations without changing molecule assignment. Checkpoints make the completed counting result recoverable even if export OOMs.

Conditional feasibility: true bounded HDF5 writes depend on backend behavior not yet verified. Exact counting under extreme barcode/UMI depth requires the separate evidence-spill work; sparse matrices alone do not solve it. Persistent storage must be on a filesystem that outlives the Slurm job, with sufficient quota; job-local scratch cannot guarantee recovery after node loss. `sync_all` and atomic rename require filesystem semantics to hold and do not supersede storage-system guarantees.

Validation unit: exact per-cell/per-gene layer element and complete molecule evidence, with performance measured per process/job on fixed inputs. Cells from one benchmark run are not independent performance replicates. Global count totals, correlations and successful file creation cannot establish count equivalence or scalability.

Scientific-critical-thinking was used as a procedural appraisal tool to distinguish observations, hypotheses and acceptance evidence; it supplies no empirical evidence about this repository. Ponytail scope: preserve assignment rules, reuse existing crates where verified, keep batch scratch initially, avoid a wholesale application rewrite, and do not add a database or distributed scheduler.
