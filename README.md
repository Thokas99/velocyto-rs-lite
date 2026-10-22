# velocyto-rs-lite

A Rust implementation of [velocyto.py](https://github.com/velocyto-team/velocyto.py) for spliced, unspliced, and ambiguous RNA counts. The project keeps the existing BAM, GTF, whitelist, and molecule-counting workflow, including BD Rhapsody `CB` and `MA` tags.

## Quick start

Build with BAM support (enabled by default):

```bash
cargo build --release
```

Building requires Rust 1.89+ and the native HDF5 development library. See the [full README](README.full.md) for platform details.

Run one BAM with its annotation and barcode whitelist:

```bash
./target/release/velocyto-rs run sample.bam genes.gtf \
  --bcfile whitelist.txt \
  --cb-tag CB --ub-tag MA \
  --sampleid sample \
  --output-format h5ad
```

The default output is H5AD. Use `--output-format loom` or `both` to choose another output.

## Large BAMs and recovery

Counting is processed in barcode batches and saved to a compressed sparse checkpoint. The default checkpoint is `<output-folder>/<sample-id>.counts`.

| Option | Purpose |
|---|---|
| `--cell-batch-size N` | Limit the number of barcodes in each batch; barcode groups stay intact. |
| `--count-memory-mb N` | Set the working-memory budget for count buffers. |
| `--checkpoint-dir PATH` | Choose where to save count progress. |
| `--resume-from PATH` | Resume matching inputs and settings from a checkpoint. |
| `--telemetry-jsonl PATH` | Write phase, memory, progress, and throughput events as JSONL. |

Memory estimates cover count buffers, not total process RSS. See `--help` for the full CLI and existing options.

Rebuild output from a completed checkpoint without rereading the BAM:

```bash
./target/release/velocyto-rs write-output results/sample.counts \
  --outputfolder results --output-format h5ad
```

## Documentation and status

- [Full README and historical project documentation](README.full.md)
- [Change summary and verification status](CHANGES.md)
- [Quantification implementation notes](docs/quantification-implementation.md)

Fixture-level differential checks are documented in `CHANGES.md`. A full production-BAM comparison and Dragon3 run are still pending.
