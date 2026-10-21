# RNA velocity quantification memory work

## What changed

Counting now writes completed barcode batches as durable sparse checkpoints instead of retaining every dense batch in memory. The existing molecule-assignment paths and BD Rhapsody `CB`/`MA` tag support remain in use. Checkpoints support counting resume, while `write-output` regenerates H5AD or loom without reopening the BAM, GTF or whitelist. H5AD layers are written incrementally; loom output uses bounded tiles. Sort spills and failed output partials remain on disk for recovery or inspection.

## Verification

- `cargo test --locked --offline`: 83 passed; 2 external-data tests ignored.
- `cargo test --locked --offline --no-default-features`: 65 passed.
- Differential fixture against the original executable: passed for seven logic variants, comparing integer layer values by cell identity.
- AnnData independently read the synthetic 187,000-cell × 30,000-gene H5AD and verified all generated integer-layer values.
- Rust formatting, all-target cargo check, CLI help and benchmark-driver self-check passed.

## Synthetic storage benchmark

Fresh release processes, 3,000 cells × 10,000 genes, one nonzero per cell in each of three layers:

| Writer | Peak process RSS | Wall time | Final directory size |
| --- | ---: | ---: | ---: |
| Dense writer | 609,844 KiB | 0.810 s | 1,839,452 B |
| Sparse checkpoint + streamed H5AD | 29,776 KiB | 0.456 s | 2,675,702 B |

The checkpoint path retains recovery files, so its disk footprint is larger. A separate 187,000 × 30,000 synthetic output run used 67,680 KiB peak process RSS, took 94.280 s, and occupied 22,406,932 B. These tests measure storage and output only; they do not parse BAMs or assign molecules. Raw reports and measurement limits are documented in [the benchmark record](docs/quantification-benchmarks.md).

## Limits

No production BD BAM, two-billion-record run or Slurm job was available. The count budget estimates reads, mapping evidence and batch scratch separately; it does not guarantee a whole-process RSS ceiling. An oversized single barcode/UMI group fails while retaining prior checkpoints. Counting resume rebuilds annotation and reruns intron validation, and partial sort runs are retained but not resumed. H5AD `X` is float32; integer layers preserve exact counts.

See the [implementation and verification record](docs/quantification-implementation.md) and [technical audit and plan](docs/memory-efficient-quantification-plan.md).
