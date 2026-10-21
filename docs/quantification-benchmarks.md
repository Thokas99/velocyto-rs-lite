# Quantification storage benchmarks

Measurements collected 2026-10-21 in fresh release processes with `scripts/benchmark_quantification.py`. The example generates synthetic counts with one nonzero per cell in each of three layers; it does not parse BAMs or run molecule assignment. The 3,000-cell dense/checkpoint runs use 10,000 genes. The larger checkpoint run uses 187,000 cells and 30,000 genes.

| Mode | Shape | Wall time | Peak process RSS | Logical I/O (`/proc`) | Final directory size |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original dense writer | 3,000 × 10,000 | 0.810 s | 609,844 KiB | 55,619 B read; 42,394 B written | 1,839,452 B |
| Sparse checkpoint + streamed H5AD | 3,000 × 10,000 | 0.456 s | 29,776 KiB | 742,997 B read; 2,188,806 B written | 2,675,702 B |
| Sparse checkpoint + streamed H5AD | 187,000 × 30,000 | 94.280 s | 67,680 KiB | 8,828,485 B read; 202,253,315 B written | 22,406,932 B |

For the 3,000-cell synthetic case, the new path used about 95% less peak process RSS and took about 44% less wall time. The checkpoint path uses more disk because it retains recoverable sparse batches alongside the H5AD. This comparison measures storage architecture on the same synthetic generator; it does not establish full-pipeline speed on biological inputs. The low density is deliberately controlled and may not match a production library.

For the 187,000-cell case, the benchmark driver was run without its repeated directory-size sampler; final size was measured afterward with `du -sb`. Consequently the raw report records zero for its disk-size fields; the table contains the separate final directory measurement. `/proc` logical writes were about 202 MB, while kernel storage block I/O was zero on this run, consistent with cached or memory-backed storage. These values do not predict writes or runtime on the Slurm filesystem.

The independent AnnData read reported shape `(187000, 30000)`, 187,000 nonzeros in each named integer layer and 561,000 nonzeros in `X`. The synthetic generator sets those exact counts by construction. This validates file readability and expected sparse shape, not biological count correctness.

The benchmark driver records peak child-process RSS and child CPU/wall time. It does not measure the complete Slurm job or concurrent process tree. A Slurm acceptance run still needs cgroup or `sacct` peak memory, production storage measurements, representative BAM counting, output-only restart, and crash recovery. Do not extrapolate the synthetic 94-second runtime to two billion alignments.

Reproduce the large synthetic writer run from the repository root:

```bash
cargo build --locked --release --example count_storage_bench
python3 scripts/benchmark_quantification.py \
  --report /persistent/bench/checkpoint-187000x30000.json \
  -- target/release/examples/count_storage_bench checkpoint \
    187000 30000 /persistent/bench/checkpoint-187000x30000
du -sb /persistent/bench/checkpoint-187000x30000
```

The three raw driver reports are preserved under [`benchmark-results`](benchmark-results/). The `dense` and `checkpoint` modes share the synthetic count generator and differ in accumulation/output path; run each in a separate fresh directory. Benchmark outputs and partial checkpoints are intentionally retained.
