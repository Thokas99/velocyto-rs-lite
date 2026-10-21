use crate::count_checkpoint::Checkpoint;
use anyhow::{ensure, Result};
use clap::Args;
use std::path::{Path, PathBuf};

#[derive(Args)]
pub struct WriteOutputArgs {
    pub checkpoint_dir: PathBuf,
    #[arg(short, long)]
    pub outputfolder: Option<PathBuf>,
    #[arg(short = 'e', long)]
    pub sampleid: Option<String>,
    #[arg(long, default_value = "h5ad", value_parser = ["h5ad", "loom", "both"])]
    pub output_format: String,
    #[arg(short, long)]
    pub verbose: bool,
}
pub fn write_output(args: WriteOutputArgs) -> Result<()> {
    let checkpoint = Checkpoint::open(&args.checkpoint_dir)?;
    let folder = args
        .outputfolder
        .as_deref()
        .unwrap_or_else(|| args.checkpoint_dir.parent().unwrap_or(Path::new(".")));
    let sample = args
        .sampleid
        .as_deref()
        .unwrap_or(&checkpoint.manifest.metadata.sample);
    export(&checkpoint, folder, sample, &args.output_format)
}
pub fn export(checkpoint: &Checkpoint, folder: &Path, sample: &str, format: &str) -> Result<()> {
    ensure!(
        matches!(format, "h5ad" | "loom" | "both"),
        "unsupported output format"
    );
    ensure!(
        checkpoint.manifest.complete(),
        "counting checkpoint is incomplete; resume counting first"
    );
    ensure!(
        !sample.is_empty()
            && Path::new(sample).components().count() == 1
            && sample != "."
            && sample != "..",
        "sample ID must be a filename stem"
    );
    std::fs::create_dir_all(folder)?;
    if matches!(format, "h5ad" | "both") {
        let path = folder.join(format!("{sample}.h5ad"));
        log::info!("Streaming H5AD output to {}", path.display());
        crate::count_output::write_h5ad(checkpoint, &path)?;
    }
    if matches!(format, "loom" | "both") {
        let path = folder.join(format!("{sample}.loom"));
        log::info!("Streaming loom output to {}", path.display());
        crate::count_output::write_loom(checkpoint, &path)?;
    }
    log::info!(
        "Count checkpoints retained at {}",
        checkpoint.root.display()
    );
    Ok(())
}
