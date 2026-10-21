# Quantification implementation and verification

The approved design is in [memory-efficient-quantification-plan.md](memory-efficient-quantification-plan.md). The work is implemented on the existing feature branch. No merge, publication or commit was performed.

## Changes

- The CLI now consumes completed count batches instead of retaining all dense batches or cloning them during transfer. The public collecting counter API remains available.
- Count batches are serialized as compressed sparse integer rows in a durable checkpoint. A small dense batch remains in use so the molecule assignment code and CB/MA behavior stay on the existing paths. The checkpoint keeps full `u32` values, a progress cursor per BAM and ordered cell/gene metadata. Version 1 inline metadata is migrated to version 2 on resume; version 2 stores static gene metadata once, outside the frequently updated progress manifest.
- H5AD and loom writers read checkpoint batches incrementally. H5AD uses checked `i64` row pointers and per-batch sparse writes; loom writes bounded tiles. Both publish through a sibling partial file and leave previous outputs and failed partials intact.
- `write-output` regenerates H5AD, loom or both from a completed checkpoint without BAM, GTF or whitelist inputs. `run --resume-from` validates its serialized settings and input path/size/mtime identities, checks committed shards and seeks to the saved BGZF cursor.
- Sorted BAMs use a completion receipt before reuse. Rust sort spill runs live in persistent directories and multi-pass merging bounds open runs. Temporary inputs are retained on success and failure.
- Count budgets bound estimates for reads, mapping evidence and batch scratch separately. They do not cap whole-process RSS. An oversized single barcode/UMI group fails with a diagnostic and leaves its checkpoint intact.

## Verification performed

- Baseline on original source: 66 passed, 1 ignored. A reference executable was built from the original revision.
- Reference differential test against that executable passed for the seven available logic variants, comparing integer layer values by cell identity. The same fixture covers BD `CB`/`MA`, UMI evidence, multiple batches, interrupted resume, and output-only restart. The external-reference test remains ignored unless `VELOCYTO_REFERENCE_BIN` is supplied.
- Current `cargo test --locked --offline`: 83 passed, 2 ignored. The ignored tests need the original external BAM fixtures and reference executable.
- Current `cargo test --locked --offline --no-default-features`: 65 passed, 0 ignored.
- The latest 187,000-cell by 30,000-gene synthetic output was opened independently with AnnData. Shape was `(187000, 30000)`; each named integer layer contained exactly 187,000 nonzeros and `X` contained 561,000 nonzeros, as designed by the synthetic fixture.
- A checkpoint migration regression verifies that v1 inline metadata resumes and atomically becomes v2. Checkpoint corruption, exact sparse round-trip, export failure retention, and multi-BAM cursor resume tests pass.
- The storage benchmark results and commands are recorded in [quantification-benchmarks.md](quantification-benchmarks.md).

## Remaining limits

No real BD BAM, two-billion-record run, Slurm allocation or production filesystem was available for validation. The synthetic benchmark measures checkpoint and H5AD storage, not BAM decoding, sorting, molecule assignment, concurrent job memory or cluster storage behavior. Per-process RSS is not Slurm job RSS.

Counting resume rebuilds the annotation index and reruns intron validation against the original BAMs. Existing partial sort runs are retained but are not resumed; the sort is repeated unless a validated completed sort receipt exists. Input identity uses canonical path, file size and modification time rather than a full content hash. A single pathological barcode or UMI is rejected if its estimate exceeds the configured component budget; disk-backed evidence processing for such groups is not implemented. Use persistent storage that survives the Slurm job and supports locking, `fsync` and atomic rename.

The critical-thinking workflow applied here was a procedural aid for separating source observations, synthetic checks and unverified production claims; it is not evidence for performance or biological validity: Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065
