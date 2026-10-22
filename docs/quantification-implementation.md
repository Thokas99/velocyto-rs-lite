# Quantification implementation and verification

The initial source audit and approved design are in [memory-efficient-quantification-plan.md](memory-efficient-quantification-plan.md). Implementation and verification details for the `feature/sparse-counting` branch are recorded here.

## Changes

- The CLI consumes count batches instead of retaining or cloning all dense batches. The collecting library API remains available.
- Batch flushes release unusually large empty `Read` buffers. GTF parsing releases the emptied raw-line vector after moving strings to sorted storage.
- Production `run` releases the original annotation model map after indexes and validation state are ready. The complete models used by counting remain available in the flat/indexed representations.
- Count batches are compressed sparse integer rows in durable checkpoints. Checkpoints retain `u32` counts, ordered cell/gene metadata, and a per-BAM cursor. Version 1 metadata migrates to version 2, where static gene metadata is written once.
- H5AD streams sparse batches; loom writes bounded tiles. Both publish through sibling partial files and preserve prior outputs on failure. `write-output` can regenerate formats from a completed checkpoint.
- Counting resume verifies settings, input path/size/mtime, and committed shards before seeking to the saved BGZF cursor. It rebuilds annotation and reruns intron validation.
- BAM CIGAR parsing no longer allocates a temporary tuple vector per record. Sort chunk estimates include sequence, CIGAR, auxiliary payload, and cached tag key sizes.
- Sorting uses bounded merge fan-in and a staged work directory. A completed staged BAM is published after provenance checks; an incomplete sort is removed and repeated. Completed sort receipts bind the output to source identity and tag.
- Optional telemetry records phases, process RSS/HWM, record/barcode progress, batch estimates, checkpoint commits, throughput, completion, and errors. Events flush when written. A telemetry write failure disables telemetry and does not abort quantification.
- `--version` reports the package version. No Cargo dependencies or lockfile entries were added. The executable remains dynamically linked to the HDF5 library available at build time.

## Verification

- Reference executable was built from source commit `17df61dae812f132c5b955b55b6e6e3efc12a1b7`.
- Differential tests passed for 8 logic variants × 3 generated BAM fixtures (24 runs). Counts match exactly by CellID, unique gene accession, and integer layer for both H5AD and loom. Fixtures cover transcript overlaps, repeated UMIs, spliced/unspliced/ambiguous evidence, high-depth barcode and UMI groups, malformed/missing tags, whitelist filtering, and independent sort paths.
- Default tests: 90 passed, 2 ignored. No-default-feature tests: 66 passed. The ignored tests require an external reference executable or large real BAM fixture.
- `cargo fmt --all -- --check`, all-target checks with and without BAM support, and release build passed. Clippy completes with legacy warnings; strict warning-free linting is not established.
- Synthetic H5AD storage benchmark: 3,000 × 10,000 dense writer used 609,844 KiB peak RSS; sparse checkpoint + streamed H5AD used 29,776 KiB. This tests output storage, not BAM counting or production memory.
- Local executables and build sidecars are excluded under `artifacts/`; no installable package is published.

## Remaining limits

No full production BAM, Dragon3 allocation, or production filesystem was available. Fixture comparisons do not establish equivalence on the complete production dataset. The synthetic benchmark does not measure BAM decoding, sorting, molecule assignment, total process RSS, or Slurm cgroup memory.

`--count-memory-mb` limits estimated read, evidence, and scratch components, not total RSS. A single oversized barcode/UMI group can fail explicitly; evidence spilling to disk is not implemented. Resume reruns annotation parsing and intron validation. Input identity uses canonical path, file size, and modification time rather than a full content hash. Atomic rename, fsync, and locking depend on filesystem support. H5AD `X` is float32; integer layers retain exact values. HDF5/system-library compatibility on Dragon3 remains unverified.

The critical-thinking workflow was used to separate observed source behavior, fixture evidence, synthetic measurement, and unverified production claims. It is not empirical evidence of performance or biological validity: Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. [arXiv:2609.00065](https://arxiv.org/abs/2609.00065).
