# Change summary

**Scope:** original Rust project at `17df61d` to `feature/sparse-counting` at `8f0888f`, plus the current uncommitted audit fixes. The molecule-classification and UMI-counting algorithms remain unchanged.

## Changes

| Area | Change | Quality evidence |
|---|---|---|
| Counting memory | Count bounded barcode batches and save compressed sparse checkpoints; export H5AD and loom from checkpoint data. | Fixture count comparisons; synthetic storage benchmark only. |
| Recovery | Validate checkpoint shards, input metadata, and sort receipts; resume committed count batches and recover or restart interrupted sorting. | Corruption, interruption, and resume tests. |
| Memory use | Drop annotation-map copies after indexes are built; free emptied GTF line buffers; release oversized read buffers after batch flush. | Code-path audit; exact output regression tests. |
| BAM hot path | Parse CIGAR operations as a stream; include CIGAR, auxiliary payloads, and cached keys in sort chunk estimates. | CIGAR regression and sort-budget tests. |
| Telemetry | Optional JSONL phases, RSS/HWM, records, barcodes, batch estimates, commits, throughput, elapsed time, and errors. | JSONL field/flush tests; write failures disable telemetry without stopping counts. |
| Compatibility | Keep BAM/GTF/whitelist inputs, existing CLI options and counting logic; `--cb-tag CB --ub-tag MA` supports BD Rhapsody. | 24 reference CLI comparisons across 8 logic variants and 3 fixtures. |
| Documentation | Add this change record and retain the prior full README as [`README.full.md`](README.full.md); replace the root README with a short guide. | README Rust example corrected; old sort/recovery claims updated. |

## Verification status

| Check | Status | Evidence and limits |
|---|---|---|
| Unit and regression tests | **Pass** | 90 passed; 2 ignored (external reference and large real-BAM fixture). |
| Differential counts | **Pass on fixtures** | 24 runs vs the original executable; exact integer layer counts aligned by CellID and unique gene accession in H5AD and loom. Not a full production BAM. |
| Graph impact | **Critical aggregate risk** | `detect_changes --scope all` reports 71 changed symbols and 34 affected execution flows; sort, count-stream, and annotation flows have broad callers. Relevant fixture and recovery tests pass. |
| Sparse export | **Pass** | H5AD and loom counts agree on fixtures; integer layers remain exact. |
| Storage benchmark | **Measured, synthetic only** | 3,000 × 10,000 storage/export fixture used about 95% less peak RSS with sparse streaming; this does not predict full-pipeline RSS. |
| Full production run | **Pending** | No representative production BAM or Slurm allocation was tested. |
| Build and lint | **Pass locally** | Formatting, default/no-default all-target checks, release build, and packaging list pass. Clippy completes with repository warnings. |
| Dragon3 compatibility | **Pending** | Build metadata records local dynamic-library requirements; Dragon3 runtime libraries have not been checked. |

## Packaging

Local Linux executables and their build sidecars are ignored under `artifacts/`; they are not release packages. An installable artifact and target-cluster compatibility check remain future work.

**Writing reference:** Kassis et al. (2026), *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents* ([arXiv:2609.00065](https://arxiv.org/abs/2609.00065)).
