//! Synthetic count-storage/export benchmark, not a BAM or biological validation.
//! Run in a fresh process for each mode: dense|checkpoint CELLS GENES OUTPUT_DIR.
use anyhow::{ensure, Context, Result};
use std::{collections::HashMap, path::Path};
use velocyto::{
    count_checkpoint::{Checkpoint, Gene, InputIdentity, Metadata},
    counter::ExInCounter,
    logic::logic_from_name,
};

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    ensure!(
        args.len() == 5,
        "usage: count_storage_bench dense|checkpoint CELLS GENES OUTPUT_DIR"
    );
    let cells: usize = args[2].parse()?;
    let genes: usize = args[3].parse()?;
    ensure!(
        cells > 0 && genes > 0 && genes <= 1_000_000,
        "positive dimensions required; genes <= 1,000,000"
    );
    ensure!(
        matches!(args[1].as_str(), "dense" | "checkpoint"),
        "unknown mode"
    );
    let folder = Path::new(&args[4]);
    std::fs::create_dir(folder).context("choose a new output directory for each benchmark")?;
    let metadata = Metadata {
        sample: "synthetic".into(),
        gem: "".into(),
        dtype: "uint32".into(),
        layers: vec!["ambiguous".into(), "spliced".into(), "unspliced".into()],
        genes: (0..genes)
            .map(|g| Gene {
                id: format!("G{g}"),
                name: format!("gene{g}"),
                chrom: "1".into(),
                strand: "+".into(),
                start: 1,
                end: 100,
            })
            .collect(),
        inputs: vec![InputIdentity {
            path: "synthetic".into(),
            bytes: 0,
            modified_ns: 0,
        }],
        settings: format!("synthetic cells={cells} genes={genes} nonzeros_per_layer_per_cell=1"),
    };
    let mut counter = ExInCounter::new(
        "synthetic".into(),
        logic_from_name("Default"),
        None,
        "no",
        false,
        "0",
        folder.to_string_lossy().into(),
        "uint32".into(),
    )?;
    for (i, gene) in metadata.genes.iter().enumerate() {
        counter.geneid2ix.insert(gene.id.clone(), i);
        counter.genes.insert(
            gene.id.clone(),
            velocyto::gene_info::GeneInfo::new(gene.name.clone(), gene.id.clone(), "1+", 1, 100),
        );
    }
    let mut checkpoint = if args[1] == "checkpoint" {
        Some(Checkpoint::create(&folder.join("counts"), metadata)?)
    } else {
        None
    };
    let mut dense: HashMap<String, Vec<ndarray::Array2<u32>>> = HashMap::new();
    let mut all_keys = vec![];
    for start in (0..cells).step_by(100) {
        let width = (cells - start).min(100);
        let mut layers = HashMap::new();
        for (l, name) in ["ambiguous", "spliced", "unspliced"].iter().enumerate() {
            let mut array = ndarray::Array2::zeros((genes, width));
            for c in 0..width {
                array[((start + c + l) % genes, c)] = (l + 1) as u32;
            }
            layers.insert(name.to_string(), array);
        }
        let keys: Vec<_> = (start..start + width)
            .map(|c| format!("|cell{c}"))
            .collect();
        if let Some(c) = &mut checkpoint {
            c.commit_batch(0, 0, &layers, &keys)?;
        } else {
            for (name, arr) in layers {
                dense.entry(name).or_default().push(arr);
            }
            all_keys.extend(keys.into_iter().map(|k| format!("synthetic:{}", &k[1..])));
        }
    }
    if let Some(mut checkpoint) = checkpoint {
        checkpoint.finish_source(0)?;
        drop(counter);
        velocyto::count_output::write_h5ad(&checkpoint, &folder.join("synthetic.h5ad"))?;
    } else {
        counter.dump_anndata(
            folder.join("synthetic.h5ad").to_str().unwrap(),
            &dense,
            &all_keys,
            None,
        )?;
    }
    Ok(())
}
