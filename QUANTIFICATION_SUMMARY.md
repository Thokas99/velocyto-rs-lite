# RNA velocity quantification: change and verification summary

## Changes

| Area | Implementation |
|---|---|
| Counting | Bounded barcode batches; compressed sparse checkpoints; resumable counting and output-only export. |
| Output | H5AD sparse writes and tiled loom export; integer count layers preserved as `u32`. |
| Memory | Release oversized batch buffers, GTF line capacity, and duplicated annotation models before counting. |
| BAM path | Stream CIGAR parsing; account for CIGAR, auxiliary payloads, and cached tag keys in sort chunk estimates. |
| Recovery | Validate count shards and sort provenance; publish completed sorts atomically; restart incomplete sorts safely. |
| Telemetry | Optional JSONL phase, RSS/HWM, progress, batch estimates, checkpoint, throughput, completion, and error events. Write failures disable telemetry. |
| Compatibility | Keep molecule assignment and UMI logic; preserve BAM/GTF/whitelist workflow and BD Rhapsody `CB`/`MA` tag overrides. |

No new Cargo dependencies were added. Output uses the existing native HDF5 runtime dependency.

## Verification

| Check | Status | Evidence |
|---|---|---|
| Default test suite | **Pass** | 90 passed; 2 ignored. |
| No-default-feature suite | **Pass** | 66 passed. |
| Differential counts | **Pass on fixtures** | 24 runs against original `17df61d`; 8 logic variants × 3 BAM fixtures; exact layer counts in H5AD and loom, aligned by CellID and accession. |
| Formatting and compile | **Pass** | `cargo fmt --all -- --check`; all-target checks with and without BAM support. |
| Clippy | **Pass with warnings** | Legacy lint warnings remain; no new warning was identified in the audit fixes. |
| Storage benchmark | **Measured, synthetic only** | Sparse checkpoint + streamed H5AD used 29,776 KiB vs 609,844 KiB dense in the 3,000 × 10,000 fixture. |
| Production BAM / Dragon3 | **Pending** | No full BAM or cluster allocation was available. |

The synthetic benchmark measures storage and output, not BAM parsing or whole-pipeline RSS. See the [benchmark record](docs/quantification-benchmarks.md) and [implementation details](docs/quantification-implementation.md).

## Remaining limits

- `--count-memory-mb` limits estimated count buffers, not whole-process or Slurm cgroup RSS.
- Flat transcript models and feature indexes still retain separate feature copies; broader de-duplication would alter core representations.
- An unusually large single barcode/UMI group can exceed the budget; it fails while retaining completed checkpoints. Evidence is not disk-backed.
- Resume rebuilds annotation and repeats intron validation. Sort state reuses a verified completed stage or restarts incomplete sorting.
- Checkpoint input identity uses canonical path, size, and modification time, not a content hash; source files must remain unchanged.
- H5AD `X` is float32; integer layers are exact. The artifact dynamically links HDF5 and system libraries; Dragon3 compatibility is unverified.

## Packaging

Local executables and build sidecars are ignored under `artifacts/`; they are not installable release packages. Packaging and target-cluster compatibility remain future work.
