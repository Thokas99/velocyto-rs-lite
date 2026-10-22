use crate::{counter::ExInCounter, logic::logic_from_name};
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    io::Cursor,
    path::Path,
};

pub(crate) fn fixture(dir: &Path) -> (String, String) {
    let gtf = dir.join("genes.gtf");
    std::fs::write(&gtf, concat!(
        "chr1\ttest\texon\t100\t199\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\"; exon_number \"1\";\n",
        "chr1\ttest\texon\t300\t399\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"Gene1\"; exon_number \"2\";\n",
        "chr1\ttest\texon\t100\t249\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\"; gene_name \"Gene1\"; exon_number \"1\";\n",
        "chr1\ttest\texon\t300\t399\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\"; gene_name \"Gene1\"; exon_number \"2\";\n",
    )).unwrap();
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n");
    for cell in ["a", "b", "c"] {
        for (i, (pos, cigar, umi)) in [
            (120, "30M", "exon"),
            (120, "30M", "exon"),
            (210, "30M", "mixed"),
            (260, "30M", "intron"),
            (180, "20M100N20M", "junction"),
        ]
        .iter()
        .enumerate()
        {
            sam.push_str(&format!("r{cell}{i}\t0\tchr1\t{pos}\t255\t{cigar}\t*\t0\t0\t*\t*\tCB:Z:{cell}\tMA:Z:{umi}\tNH:i:1\tST:Z:S1\n"));
        }
    }
    let bam = dir.join("cellsorted.bam");
    write_bam(&bam, &sam);
    (bam.to_string_lossy().into(), gtf.to_string_lossy().into())
}

fn write_bam(path: &Path, sam: &str) {
    use noodles_sam::alignment::io::Write;
    let mut reader = noodles_sam::io::Reader::new(Cursor::new(sam));
    let header = reader.read_header().unwrap();
    let mut writer = noodles_bam::io::Writer::new(std::fs::File::create(path).unwrap());
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer
            .write_alignment_record(&header, &record.unwrap())
            .unwrap();
    }
    writer.try_finish().unwrap();
}

fn high_depth_fixture(dir: &Path) -> (String, String) {
    let (bam, gtf) = fixture(dir);
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n");
    for umi in 0..100 {
        for read in 0..100 {
            sam.push_str(&format!(
                "r{umi}_{read}\t0\tchr1\t120\t255\t30M\t*\t0\t0\t*\t*\tCB:Z:a\tMA:Z:u{umi}\tNH:i:1\tST:Z:S1\n"
            ));
        }
    }
    for read in 0..90_000 {
        sam.push_str(&format!(
            "excluded{read}\t0\tchr1\t120\t255\t30M\t*\t0\t0\t*\t*\tCB:Z:z\tMA:Z:outside\tNH:i:1\tST:Z:S1\n"
        ));
    }
    write_bam(Path::new(&bam), &sam);
    (bam, gtf)
}

fn tag_edge_fixture(dir: &Path) -> (String, String) {
    let (bam, gtf) = fixture(dir);
    let sam = concat!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n",
        "valid\t0\tchr1\t120\t255\t30M\t*\t0\t0\t*\t*\tCB:Z:a\tMA:Z:valid\tNH:i:1\tST:Z:S1\n",
        "missing_cb\t0\tchr1\t120\t255\t30M\t*\t0\t0\t*\t*\tMA:Z:no_cb\tNH:i:1\tST:Z:S1\n",
        "missing_umi\t0\tchr1\t120\t255\t30M\t*\t0\t0\t*\t*\tCB:Z:a\tNH:i:1\tST:Z:S1\n",
        "typed_cb\t0\tchr1\t120\t255\t30M\t*\t0\t0\t*\t*\tCB:i:1\tMA:Z:typed_cb\tNH:i:1\tST:Z:S1\n",
        "typed_umi\t0\tchr1\t120\t255\t30M\t*\t0\t0\t*\t*\tCB:Z:a\tMA:i:3\tNH:i:1\tST:Z:S1\n",
    );
    write_bam(Path::new(&bam), sam);
    (bam, gtf)
}

fn counter(gtf: &str, bam: &str, logic: &str) -> ExInCounter {
    let whitelist: HashSet<_> = ["a", "b", "c"].into_iter().map(str::to_string).collect();
    let mut c = ExInCounter::new(
        "test".into(),
        logic_from_name(logic),
        Some(whitelist),
        "no",
        false,
        "0",
        "/tmp".into(),
        "uint32".into(),
    )
    .unwrap();
    c.umibarcode_str = "MA".into();
    c.read_transcriptmodels(gtf).unwrap();
    c.mark_up_introns(bam, false).unwrap();
    c
}

fn unique_ids(file: &hdf5::File, path: &str) -> Vec<String> {
    let values = file
        .dataset(path)
        .unwrap()
        .read_raw::<hdf5::types::VarLenUnicode>()
        .unwrap();
    let values: Vec<_> = values.iter().map(ToString::to_string).collect();
    assert_eq!(
        values.iter().collect::<HashSet<_>>().len(),
        values.len(),
        "{path} must be unique"
    );
    values
}

fn read_h5ad_counts(path: &Path, layer: &str) -> BTreeMap<String, BTreeMap<String, u32>> {
    let file = hdf5::File::open(path).unwrap();
    let ids = unique_ids(&file, "obs/_index");
    let genes = unique_ids(&file, "var/Accession");
    let ptr = file
        .dataset(&format!("layers/{layer}/indptr"))
        .unwrap()
        .read_raw::<i64>()
        .unwrap();
    let indices = file
        .dataset(&format!("layers/{layer}/indices"))
        .unwrap()
        .read_raw::<i32>()
        .unwrap();
    let data = file
        .dataset(&format!("layers/{layer}/data"))
        .unwrap()
        .read_raw::<u32>()
        .unwrap();
    assert_eq!(ptr.len(), ids.len() + 1);
    assert_eq!(indices.len(), data.len());
    let mut result = BTreeMap::new();
    for (cell, id) in ids.into_iter().enumerate() {
        let mut row = BTreeMap::new();
        for i in usize::try_from(ptr[cell]).unwrap()..usize::try_from(ptr[cell + 1]).unwrap() {
            let gene = genes.get(usize::try_from(indices[i]).unwrap()).unwrap();
            assert!(
                row.insert(gene.clone(), data[i]).is_none(),
                "duplicate sparse gene entry"
            );
        }
        assert!(result.insert(id, row).is_none(), "duplicate CellID");
    }
    result
}

fn read_loom_counts(path: &Path, layer: &str) -> BTreeMap<String, BTreeMap<String, u32>> {
    let file = hdf5::File::open(path).unwrap();
    let ids = unique_ids(&file, "col_attrs/CellID");
    let genes = unique_ids(&file, "row_attrs/Accession");
    let dataset = file.dataset(&format!("layers/{layer}")).unwrap();
    assert_eq!(dataset.shape(), [genes.len(), ids.len()]);
    let values = dataset.read_raw::<u32>().unwrap();
    let n_cells = ids.len();
    let mut result = BTreeMap::new();
    for (cell, id) in ids.into_iter().enumerate() {
        let mut row = BTreeMap::new();
        for (gene, accession) in genes.iter().enumerate() {
            let value = values[gene * n_cells + cell];
            if value != 0 {
                row.insert(accession.clone(), value);
            }
        }
        assert!(result.insert(id, row).is_none(), "duplicate CellID");
    }
    result
}

#[test]
fn streaming_matches_collecting_and_resumes_without_losing_lookahead() {
    let dir = tempfile::tempdir().unwrap();
    let (bam, gtf) = fixture(dir.path());
    for logic in [
        "Default",
        "Permissive10X",
        "Intermediate10X",
        "ValidatedIntrons10X",
        "Stricter10X",
        "ObservedSpanning10X",
        "Discordant10X",
        "SmartSeq2",
    ] {
        let mut original = counter(&gtf, &bam, logic);
        let (layers, keys) = original.count(&bam, false, 100).unwrap();
        let mut expected = HashMap::new();
        for (layer, arrays) in layers {
            for (col, key) in keys.iter().enumerate() {
                expected.insert((layer.clone(), key.clone()), arrays[0].column(col).to_vec());
            }
        }
        let mut optimized = counter(&gtf, &bam, logic);
        let mut observed = HashMap::new();
        let mut cursor = 0;
        let interrupted =
            optimized.count_stream(&bam, false, 1, 0, 1024 * 1024, |batch, keys, next| {
                for (layer, arr) in batch {
                    for (col, key) in keys.iter().enumerate() {
                        observed.insert((layer.clone(), key.clone()), arr.column(col).to_vec());
                    }
                }
                cursor = next;
                anyhow::bail!("simulated crash after checkpoint")
            });
        assert!(interrupted.is_err());
        assert!(cursor > 0);
        optimized
            .count_stream(&bam, false, 1, cursor, 1024 * 1024, |batch, keys, _| {
                for (layer, arr) in batch {
                    for (col, key) in keys.iter().enumerate() {
                        assert!(observed
                            .insert((layer.clone(), key.clone()), arr.column(col).to_vec())
                            .is_none());
                    }
                }
                Ok(())
            })
            .unwrap();
        assert_eq!(observed, expected, "logic {logic}");
    }
}

#[test]
fn count_limits_fail_explicitly_instead_of_splitting_molecules() {
    let dir = tempfile::tempdir().unwrap();
    let (bam, gtf) = fixture(dir.path());
    let mut c = counter(&gtf, &bam, "Default");
    let mut emitted = false;
    let result = c.count_stream(&bam, false, 100, 0, 1, |_, _, _| {
        emitted = true;
        Ok(())
    });
    assert!(result.unwrap_err().to_string().contains("budget"));
    assert!(!emitted);
}

#[test]
fn unfiltered_stream_does_not_cache_every_observed_barcode() {
    let dir = tempfile::tempdir().unwrap();
    let (bam, gtf) = fixture(dir.path());
    let mut c = counter(&gtf, &bam, "Default");
    c.filter_mode = false;
    c.valid_bcset.clear();
    c.count_stream(&bam, false, 1, 0, 1_048_576, |_, keys, _| {
        assert!(
            keys.is_empty(),
            "fixture cells have fewer than 81 molecules"
        );
        Ok(())
    })
    .unwrap();
    assert!(
        c.valid_bcset.is_empty(),
        "unfiltered runs must not retain all raw barcodes"
    );
}

#[test]
fn cli_checkpoint_resume_and_output_only_without_inputs() {
    use clap::Parser;
    #[derive(Parser)]
    struct Cli {
        #[command(flatten)]
        args: crate::commands::run::RunArgs,
    }
    let dir = tempfile::tempdir().unwrap();
    let (bam, gtf) = fixture(dir.path());
    let whitelist = dir.path().join("barcodes.txt");
    std::fs::write(&whitelist, "a\nb\nc\n").unwrap();
    let folder = dir.path().join("out");
    let telemetry = dir.path().join("telemetry.jsonl");
    let args = vec![
        "cmd",
        &bam,
        &gtf,
        "-b",
        whitelist.to_str().unwrap(),
        "-o",
        folder.to_str().unwrap(),
        "-e",
        "test",
        "--cb-tag",
        "CB",
        "--ub-tag",
        "MA",
        "--sample-tag",
        "ST",
        "--samtools-threads",
        "1",
        "--cell-batch-size",
        "1",
        "--telemetry-jsonl",
        telemetry.to_str().unwrap(),
    ];
    crate::commands::run::run(Cli::parse_from(&args).args).unwrap();
    let events: Vec<serde_json::Value> = std::fs::read_to_string(&telemetry)
        .unwrap()
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert!(events.iter().any(|e| e["phase"] == "bam_sort"));
    assert!(events.iter().any(|e| e["phase"] == "counting"));
    assert!(events.iter().any(|e| e["event"] == "checkpoint_commit"));
    assert!(events.iter().any(|e| e["event"] == "run_complete"));
    let commit = events
        .iter()
        .find(|e| e["event"] == "checkpoint_commit")
        .unwrap();
    assert!(commit["batch_size"]["records"].is_number());
    assert!(commit["estimated_batch_memory_bytes"]["read_batch"].is_number());
    assert!(commit["estimated_batch_memory_bytes"]["dense_scratch"].is_number());
    assert!(events
        .iter()
        .all(|e| e["rss_kib"].is_number() && e["peak_rss_kib"].is_number()));
    let root = folder.join("test.counts");
    let original = crate::count_checkpoint::Checkpoint::open(&root).unwrap();
    assert!(original.manifest.complete());
    assert_eq!(original.manifest.rows().unwrap(), 3);
    let n_batches = original.manifest.batches.len();
    let mut resume = args;
    resume.extend(["--resume-from", root.to_str().unwrap()]);
    crate::commands::run::run(Cli::parse_from(resume).args).unwrap();
    assert_eq!(
        crate::count_checkpoint::Checkpoint::open(&root)
            .unwrap()
            .manifest
            .batches
            .len(),
        n_batches
    );
    std::fs::rename(&bam, format!("{bam}.offline")).unwrap();
    std::fs::rename(&gtf, format!("{gtf}.offline")).unwrap();
    std::fs::remove_file(&whitelist).unwrap();
    crate::commands::write_output::write_output(crate::commands::write_output::WriteOutputArgs {
        checkpoint_dir: root,
        outputfolder: Some(dir.path().join("recovered")),
        sampleid: None,
        output_format: "both".into(),
        verbose: false,
    })
    .unwrap();
    let file = hdf5::File::open(dir.path().join("recovered/test.h5ad")).unwrap();
    assert_eq!(file.dataset("obs/_index").unwrap().size(), 3);
}

#[test]
fn fresh_counter_resumes_second_bam_with_interleaved_samples() {
    use crate::count_checkpoint::{Checkpoint, InputIdentity, Metadata};
    use noodles_sam::alignment::io::Write;
    let dir = tempfile::tempdir().unwrap();
    let (bam, gtf) = fixture(dir.path());
    let interleaved = dir.path().join("interleaved.bam");
    let mut reader = noodles_bam::io::Reader::new(std::fs::File::open(&bam).unwrap());
    let header = reader.read_header().unwrap();
    let mut writer = noodles_bam::io::Writer::new(std::fs::File::create(&interleaved).unwrap());
    writer.write_header(&header).unwrap();
    let mut record = noodles_sam::alignment::RecordBuf::default();
    let mut index = 0;
    while reader.read_record_buf(&header, &mut record).unwrap() != 0 {
        record.data_mut().insert(
            noodles_sam::alignment::record::data::field::Tag::from(*b"ST"),
            noodles_sam::alignment::record_buf::data::field::Value::String(
                if index % 2 == 0 { "S1" } else { "S2" }.into(),
            ),
        );
        writer.write_alignment_record(&header, &record).unwrap();
        index += 1;
    }
    writer.try_finish().unwrap();
    let mut c = counter(&gtf, interleaved.to_str().unwrap(), "Default");
    c.sample_tag = Some("ST".into());
    let metadata = Metadata::from_counter(
        &c,
        "".into(),
        vec![
            InputIdentity::read(&interleaved).unwrap(),
            InputIdentity::read(&interleaved).unwrap(),
        ],
        "fixture".into(),
    )
    .unwrap();
    let root = dir.path().join("counts");
    let mut checkpoint = Checkpoint::create(&root, metadata.clone()).unwrap();
    c.count_stream(
        interleaved.to_str().unwrap(),
        false,
        1,
        0,
        1_048_576,
        |layers, keys, next| checkpoint.commit_batch(0, next, &layers, &keys),
    )
    .unwrap();
    checkpoint.finish_source(0).unwrap();
    let result = c.count_stream(
        interleaved.to_str().unwrap(),
        false,
        1,
        0,
        1_048_576,
        |layers, keys, next| {
            checkpoint.commit_batch(1, next, &layers, &keys)?;
            anyhow::bail!("simulated crash in second source")
        },
    );
    assert!(result.is_err());
    assert_eq!(checkpoint.manifest.rows().unwrap(), 8); // six cells in source 1 plus first two in source 2
    drop(checkpoint);
    drop(c);
    let mut c = counter(&gtf, interleaved.to_str().unwrap(), "Default");
    c.sample_tag = Some("ST".into());
    let mut checkpoint = Checkpoint::resume(&root, &metadata).unwrap();
    c.count_stream(
        interleaved.to_str().unwrap(),
        false,
        1,
        checkpoint.manifest.sources[1].offset,
        1_048_576,
        |layers, keys, next| checkpoint.commit_batch(1, next, &layers, &keys),
    )
    .unwrap();
    checkpoint.finish_source(1).unwrap();
    assert!(checkpoint.manifest.complete());
    assert_eq!(checkpoint.manifest.rows().unwrap(), 12);
    let mut observed = Vec::new();
    for i in 0..checkpoint.manifest.batches.len() {
        let batch = checkpoint.read_batch(i).unwrap();
        for (l, name) in metadata.layers.iter().enumerate() {
            for (col, key) in batch.keys.iter().enumerate() {
                observed.push((name.clone(), key.clone(), batch.layers[l][col].clone()));
            }
        }
    }
    let (layers, keys) = c.count(interleaved.to_str().unwrap(), false, 100).unwrap();
    let mut expected = Vec::new();
    for _ in 0..2 {
        for (name, arrays) in &layers {
            for (col, key) in keys.iter().enumerate() {
                let row: Vec<_> = arrays[0]
                    .column(col)
                    .iter()
                    .enumerate()
                    .filter_map(|(g, &v)| (v > 0).then_some((g as u32, v)))
                    .collect();
                expected.push((name.clone(), key.clone(), row));
            }
        }
    }
    observed.sort();
    expected.sort();
    assert_eq!(observed, expected);
}

#[test]
#[ignore = "requires VELOCYTO_REFERENCE_BIN pointing to the original 17df61d executable"]
fn differential_against_original_cli() {
    use clap::Parser;
    let executable = std::env::var("VELOCYTO_REFERENCE_BIN").expect("set VELOCYTO_REFERENCE_BIN");
    let root = tempfile::tempdir().unwrap();
    type FixtureFn = fn(&Path) -> (String, String);
    let cases: [(&str, FixtureFn); 3] = [
        ("core", fixture),
        ("high-depth", high_depth_fixture),
        ("tag-edges", tag_edge_fixture),
    ];
    let logics = [
        "Default",
        "Permissive10X",
        "Intermediate10X",
        "ValidatedIntrons10X",
        "Stricter10X",
        "ObservedSpanning10X",
        "Discordant10X",
        "SmartSeq2",
    ];
    for (case, make_fixture) in cases {
        let case_dir = root.path().join(case);
        std::fs::create_dir_all(&case_dir).unwrap();
        let (bam, gtf) = make_fixture(&case_dir);
        let whitelist = case_dir.join("barcodes.txt");
        std::fs::write(&whitelist, "a\nb\nc\n").unwrap();
        for logic in logics {
            let reference_out = case_dir.join(format!("reference-{logic}"));
            let optimized_out = case_dir.join(format!("optimized-{logic}"));
            let result = std::process::Command::new(&executable)
                .args([
                    "run",
                    &bam,
                    &gtf,
                    "-b",
                    whitelist.to_str().unwrap(),
                    "-o",
                    reference_out.to_str().unwrap(),
                    "-e",
                    "reference",
                    "-l",
                    logic,
                    "--cb-tag",
                    "CB",
                    "--ub-tag",
                    "MA",
                    "--sample-tag",
                    "ST",
                    "--samtools-threads",
                    "1",
                    "--output-format",
                    "both",
                ])
                .output()
                .unwrap();
            assert!(
                result.status.success(),
                "{case}/{logic}: {}",
                String::from_utf8_lossy(&result.stderr)
            );
            // The two implementations must sort independently from the same raw BAM.
            let sorted = case_dir.join("cellsorted_cellsorted.bam");
            if sorted.exists() {
                std::fs::remove_file(&sorted).unwrap();
            }
            let receipt = sorted.with_extension("bam.sort.json");
            if receipt.exists() {
                std::fs::remove_file(receipt).unwrap();
            }
            #[derive(clap::Parser)]
            struct Cli {
                #[command(flatten)]
                args: crate::commands::run::RunArgs,
            }
            let telemetry = case_dir.join(format!("telemetry-{logic}.jsonl"));
            let args = Cli::parse_from([
                "cmd",
                &bam,
                &gtf,
                "-b",
                whitelist.to_str().unwrap(),
                "-o",
                optimized_out.to_str().unwrap(),
                "-e",
                "reference",
                "-l",
                logic,
                "--cb-tag",
                "CB",
                "--ub-tag",
                "MA",
                "--sample-tag",
                "ST",
                "--samtools-threads",
                "1",
                "--samtools-memory",
                "256",
                "--cell-batch-size",
                "1",
                "--count-memory-mb",
                "256",
                "--output-format",
                "both",
                "--telemetry-jsonl",
                telemetry.to_str().unwrap(),
            ])
            .args;
            crate::commands::run::run(args).unwrap();
            let checkpoint =
                crate::count_checkpoint::Checkpoint::open(&optimized_out.join("reference.counts"))
                    .unwrap();
            let events: Vec<serde_json::Value> = std::fs::read_to_string(&telemetry)
                .unwrap()
                .lines()
                .map(|line| serde_json::from_str(line).unwrap())
                .collect();
            assert!(events.iter().any(|e| e["event"] == "checkpoint_commit"));
            assert!(events.iter().any(|e| e["event"] == "run_complete"));
            if case == "high-depth" {
                assert!(events.iter().any(|e| {
                    e["event"] == "count_progress"
                        && e["records_processed"].as_u64().unwrap_or(0) >= 100_000
                        && e["records_per_second"].is_number()
                }));
            }
            for layer in &checkpoint.manifest.metadata.layers {
                let reference_h5ad = read_h5ad_counts(&reference_out.join("reference.h5ad"), layer);
                let optimized_h5ad = read_h5ad_counts(&optimized_out.join("reference.h5ad"), layer);
                let reference_loom = read_loom_counts(&reference_out.join("reference.loom"), layer);
                let optimized_loom = read_loom_counts(&optimized_out.join("reference.loom"), layer);
                assert_eq!(
                    reference_h5ad, optimized_h5ad,
                    "H5AD {case}/{logic}/{layer}"
                );
                assert_eq!(
                    reference_loom, optimized_loom,
                    "loom {case}/{logic}/{layer}"
                );
                assert_eq!(
                    reference_h5ad, reference_loom,
                    "reference format counts {case}/{logic}/{layer}"
                );
                assert_eq!(
                    optimized_h5ad, optimized_loom,
                    "optimized format counts {case}/{logic}/{layer}"
                );
                if case == "core" && logic == "Default" {
                    let total: u64 = reference_h5ad
                        .values()
                        .flat_map(BTreeMap::values)
                        .map(|&v| u64::from(v))
                        .sum();
                    assert!(total > 0, "core fixture did not exercise {layer}");
                }
                if case == "high-depth" && logic == "Default" {
                    assert_eq!(
                        reference_h5ad.len(),
                        1,
                        "whitelist admitted an outside barcode"
                    );
                    assert_eq!(reference_h5ad.keys().next().unwrap(), "reference_S1:ax");
                }
            }
        }
    }
}
