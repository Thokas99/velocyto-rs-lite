//! Durable sparse count batches. Complete shards are immutable; failures never clean them up.
use anyhow::{bail, ensure, Context, Result};
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    fs::{self, File, OpenOptions},
    io::{BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};

pub type DenseLayers = HashMap<String, Array2<u32>>;
pub type SparseRow = Vec<(u32, u32)>;
const VERSION: u32 = 2;
const MAX_SHARD_BYTES: u64 = 512 * 1024 * 1024;
const MAX_MANIFEST_BYTES: u64 = 64 * 1024 * 1024;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Gene {
    pub id: String,
    pub name: String,
    pub chrom: String,
    pub strand: String,
    pub start: i64,
    pub end: i64,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct InputIdentity {
    pub path: PathBuf,
    pub bytes: u64,
    pub modified_ns: u128,
}
impl InputIdentity {
    pub fn read(path: &Path) -> Result<Self> {
        let m = fs::metadata(path)?;
        Ok(Self {
            path: fs::canonicalize(path)?,
            bytes: m.len(),
            modified_ns: m
                .modified()?
                .duration_since(std::time::UNIX_EPOCH)?
                .as_nanos(),
        })
    }
}
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Metadata {
    pub sample: String,
    pub gem: String,
    pub dtype: String,
    pub layers: Vec<String>,
    pub genes: Vec<Gene>,
    pub inputs: Vec<InputIdentity>,
    pub settings: String,
}
impl Metadata {
    pub fn from_counter(
        counter: &crate::counter::ExInCounter,
        gem: String,
        inputs: Vec<InputIdentity>,
        settings: String,
    ) -> Result<Self> {
        let mut genes = vec![None; counter.geneid2ix.len()];
        for (id, &i) in &counter.geneid2ix {
            ensure!(i < genes.len(), "gene index out of bounds");
            let g = counter.genes.get(id).context("missing gene metadata")?;
            ensure!(genes[i].is_none(), "duplicate gene index");
            genes[i] = Some(Gene {
                id: id.clone(),
                name: g.genename.clone(),
                chrom: g.chrom.clone(),
                strand: g.strand.to_string(),
                start: g.start,
                end: g.end,
            });
        }
        let mut layers: Vec<_> = counter
            .logic
            .layers()
            .iter()
            .map(|s| s.to_string())
            .collect();
        layers.sort();
        Ok(Self {
            sample: counter.sampleid.clone(),
            gem,
            dtype: counter.loom_numeric_dtype.clone(),
            layers,
            genes: genes
                .into_iter()
                .map(|g| g.context("missing gene index"))
                .collect::<Result<_>>()?,
            inputs,
            settings,
        })
    }
    pub fn cell_id(&self, key: &str) -> (String, String) {
        let (sample, bc) = key.split_once('|').unwrap_or(("", key));
        let id = if sample.is_empty() {
            format!("{}:{bc}{}", self.sample, self.gem)
        } else {
            format!("{}_{sample}:{bc}{}", self.sample, self.gem)
        };
        (id, sample.into())
    }
    #[cfg(test)]
    pub(crate) fn test_metadata() -> Self {
        Self {
            sample: "test".into(),
            gem: "".into(),
            dtype: "uint32".into(),
            layers: vec!["ambiguous".into(), "spliced".into(), "unspliced".into()],
            genes: (0..2)
                .map(|i| Gene {
                    id: format!("G{i}"),
                    name: format!("gene{i}"),
                    chrom: "1".into(),
                    strand: "+".into(),
                    start: 1,
                    end: 100,
                })
                .collect(),
            inputs: vec![InputIdentity {
                path: "/test.bam".into(),
                bytes: 0,
                modified_ns: 0,
            }],
            settings: "test".into(),
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Progress {
    pub offset: u64,
    pub complete: bool,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BatchInfo {
    pub file: String,
    pub rows: usize,
    pub nnz: Vec<u64>,
    pub totals: Vec<u64>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Manifest {
    pub version: u32,
    #[serde(skip_serializing, default)]
    pub metadata: Metadata,
    pub sources: Vec<Progress>,
    pub batches: Vec<BatchInfo>,
}
impl Manifest {
    pub fn complete(&self) -> bool {
        self.sources.iter().all(|s| s.complete)
    }
    pub fn rows(&self) -> Result<usize> {
        self.batches.iter().try_fold(0usize, |a, b| {
            a.checked_add(b.rows).context("cell count overflow")
        })
    }
}
pub struct SparseBatch {
    pub keys: Vec<String>,
    pub layers: Vec<Vec<SparseRow>>,
}
pub struct Checkpoint {
    pub root: PathBuf,
    pub manifest: Manifest,
    _lock: Option<File>,
}
impl Checkpoint {
    pub fn create(root: &Path, metadata: Metadata) -> Result<Self> {
        ensure!(!metadata.inputs.is_empty(), "no counting inputs");
        fs::create_dir(root).with_context(|| {
            format!(
                "Creating checkpoint {}; choose a new path or use --resume-from",
                root.display()
            )
        })?;
        let lock = Self::lock(root)?;
        let sources = metadata
            .inputs
            .iter()
            .map(|_| Progress {
                offset: 0,
                complete: false,
            })
            .collect();
        let out = Self {
            root: root.into(),
            manifest: Manifest {
                version: VERSION,
                metadata,
                sources,
                batches: vec![],
            },
            _lock: Some(lock),
        };
        out.save_metadata()?;
        out.save()?;
        Ok(out)
    }
    fn lock(root: &Path) -> Result<File> {
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(root.join("writer.lock"))?;
        file.try_lock()
            .context("checkpoint is already in use by a counting writer")?;
        Ok(file)
    }
    pub fn open(root: &Path) -> Result<Self> {
        let f = File::open(root.join("manifest.json"))?;
        ensure!(
            f.metadata()?.len() <= MAX_MANIFEST_BYTES,
            "manifest exceeds size limit"
        );
        let mut manifest: Manifest = serde_json::from_reader(BufReader::new(f))?;
        ensure!(
            matches!(manifest.version, 1 | VERSION),
            "unsupported checkpoint version"
        );
        if manifest.version == VERSION {
            let metadata = File::open(root.join("metadata.json"))?;
            ensure!(
                metadata.metadata()?.len() <= MAX_MANIFEST_BYTES,
                "metadata exceeds size limit"
            );
            manifest.metadata = serde_json::from_reader(BufReader::new(metadata))?;
        }
        ensure!(
            !manifest.sources.is_empty()
                && manifest.sources.len() == manifest.metadata.inputs.len(),
            "invalid source manifest"
        );
        ensure!(
            manifest.metadata.genes.len() <= 1_000_000 && manifest.metadata.layers.len() <= 16,
            "checkpoint dimensions exceed limits"
        );
        ensure!(
            matches!(manifest.metadata.dtype.as_str(), "uint16" | "uint32"),
            "invalid checkpoint dtype"
        );
        for b in &manifest.batches {
            ensure!(
                b.file.starts_with("batch-")
                    && b.file.ends_with(".gz")
                    && Path::new(&b.file).components().count() == 1,
                "invalid shard filename"
            );
            ensure!(
                b.rows <= 100_000
                    && b.nnz.len() == manifest.metadata.layers.len()
                    && b.totals.len() == b.nnz.len(),
                "invalid batch metadata"
            );
        }
        Ok(Self {
            root: root.into(),
            manifest,
            _lock: None,
        })
    }
    pub fn resume(root: &Path, metadata: &Metadata) -> Result<Self> {
        let lock = Self::lock(root)?;
        let mut out = Self::open(root)?;
        ensure!(
            &out.manifest.metadata == metadata,
            "checkpoint inputs, genes or settings differ; refusing resume"
        );
        out._lock = Some(lock);
        for i in 0..out.manifest.batches.len() {
            out.read_batch(i)?;
        }
        if out.manifest.version == 1 {
            out.save_metadata()?;
            out.manifest.version = VERSION;
            out.save()?;
        }
        Ok(out)
    }
    fn save_metadata(&self) -> Result<()> {
        let partial = self
            .root
            .join(format!("metadata-{:016x}.partial", rand::random::<u64>()));
        let mut f = BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&partial)?,
        );
        serde_json::to_writer(&mut f, &self.manifest.metadata)?;
        f.flush()?;
        ensure!(
            f.get_ref().metadata()?.len() <= MAX_MANIFEST_BYTES,
            "metadata exceeds size limit"
        );
        f.get_ref().sync_all()?;
        fs::rename(partial, self.root.join("metadata.json"))?;
        File::open(&self.root)?.sync_all()?;
        Ok(())
    }
    fn save(&self) -> Result<()> {
        let partial = self
            .root
            .join(format!("manifest-{:016x}.partial", rand::random::<u64>()));
        let mut f = BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&partial)?,
        );
        serde_json::to_writer(&mut f, &self.manifest)?;
        f.flush()?;
        ensure!(
            f.get_ref().metadata()?.len() <= MAX_MANIFEST_BYTES,
            "manifest exceeds size limit"
        );
        f.get_ref().sync_all()?;
        fs::rename(partial, self.root.join("manifest.json"))?;
        File::open(&self.root)?.sync_all()?;
        Ok(())
    }
    pub fn commit_batch(
        &mut self,
        source: usize,
        offset: u64,
        layers: &DenseLayers,
        keys: &[String],
    ) -> Result<()> {
        ensure!(self._lock.is_some(), "checkpoint opened read-only");
        ensure!(
            source < self.manifest.sources.len() && !self.manifest.sources[source].complete,
            "invalid counting source"
        );
        let meta = &self.manifest.metadata;
        ensure!(
            layers.len() == meta.layers.len() && keys.len() <= 100_000,
            "invalid batch shape"
        );
        let mut size = 16u64;
        for key in keys {
            ensure!(key.len() <= 1_048_576, "cell identity too long");
            size += 8 + key.len() as u64;
        }
        let mut nnz = vec![];
        let mut totals = vec![];
        for name in &meta.layers {
            let arr = layers.get(name).context("missing counting layer")?;
            ensure!(
                arr.dim() == (meta.genes.len(), keys.len()),
                "count dimensions differ from metadata"
            );
            let n = arr.iter().filter(|&&v| v != 0).count() as u64;
            nnz.push(n);
            totals.push(arr.iter().try_fold(0u64, |s, &v| {
                s.checked_add(v as u64).context("layer total overflow")
            })?);
            size = size
                .checked_add(keys.len() as u64 * 8 + n * 8)
                .context("shard size overflow")?;
        }
        ensure!(
            size <= MAX_SHARD_BYTES,
            "count batch exceeds 512 MiB checkpoint limit; reduce --cell-batch-size"
        );
        let file = format!(
            "batch-{:08}-{:016x}.gz",
            self.manifest.batches.len(),
            rand::random::<u64>()
        );
        let partial = self.root.join(format!("{file}.partial"));
        let mut w = GzEncoder::new(
            BufWriter::new(
                OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(&partial)?,
            ),
            Compression::new(2),
        );
        w.write_all(b"VRSCOUNT")?;
        put_u64(&mut w, keys.len() as u64)?;
        for key in keys {
            put_u64(&mut w, key.len() as u64)?;
            w.write_all(key.as_bytes())?;
        }
        for name in &meta.layers {
            let arr = &layers[name];
            for c in 0..keys.len() {
                put_u64(
                    &mut w,
                    arr.column(c).iter().filter(|&&v| v != 0).count() as u64,
                )?;
                for g in 0..arr.nrows() {
                    let v = arr[(g, c)];
                    if v != 0 {
                        w.write_all(&(g as u32).to_le_bytes())?;
                        w.write_all(&v.to_le_bytes())?;
                    }
                }
            }
        }
        let mut w = w.finish()?;
        w.flush()?;
        w.get_ref().sync_all()?;
        fs::rename(&partial, self.root.join(&file))?;
        File::open(&self.root)?.sync_all()?;
        self.manifest.batches.push(BatchInfo {
            file,
            rows: keys.len(),
            nnz,
            totals,
        });
        self.manifest.sources[source].offset = offset;
        self.save()
    }
    pub fn finish_source(&mut self, source: usize) -> Result<()> {
        ensure!(
            self._lock.is_some() && source < self.manifest.sources.len(),
            "invalid counting source"
        );
        self.manifest.sources[source].complete = true;
        self.save()
    }
    pub fn read_batch(&self, index: usize) -> Result<SparseBatch> {
        let info = self.manifest.batches.get(index).context("missing batch")?;
        let mut r = GzDecoder::new(BufReader::new(File::open(self.root.join(&info.file))?))
            .take(MAX_SHARD_BYTES + 1);
        let mut magic = [0; 8];
        r.read_exact(&mut magic)?;
        ensure!(&magic == b"VRSCOUNT", "invalid checkpoint shard");
        ensure!(
            get_u64(&mut r)? == info.rows as u64,
            "checkpoint row count mismatch"
        );
        let mut keys = Vec::with_capacity(info.rows);
        for _ in 0..info.rows {
            let n = get_u64(&mut r)?;
            ensure!(n <= 1_048_576, "cell identity exceeds size limit");
            let mut key = vec![0; n as usize];
            r.read_exact(&mut key)?;
            keys.push(String::from_utf8(key)?);
        }
        let mut layers = vec![];
        for l in 0..self.manifest.metadata.layers.len() {
            let mut rows = Vec::with_capacity(info.rows);
            let (mut nnz, mut total) = (0u64, 0u64);
            for _ in 0..info.rows {
                let n = get_u64(&mut r)?;
                ensure!(
                    n <= self.manifest.metadata.genes.len() as u64,
                    "invalid row nonzero count"
                );
                let mut row = Vec::with_capacity((n as usize).min(4096));
                let mut last = None;
                for _ in 0..n {
                    let g = get_u32(&mut r)?;
                    let v = get_u32(&mut r)?;
                    ensure!(
                        (g as usize) < self.manifest.metadata.genes.len()
                            && v > 0
                            && last.is_none_or(|p| g > p),
                        "invalid sparse row"
                    );
                    last = Some(g);
                    row.push((g, v));
                    total = total
                        .checked_add(v as u64)
                        .context("layer total overflow")?;
                }
                nnz += n;
                rows.push(row);
            }
            ensure!(
                nnz == info.nnz[l] && total == info.totals[l],
                "checkpoint counts differ from manifest"
            );
            layers.push(rows);
        }
        let mut end = [0];
        ensure!(
            r.read(&mut end)? == 0 && r.limit() > 0,
            "checkpoint trailing data or size limit exceeded"
        );
        Ok(SparseBatch { keys, layers })
    }
}
fn put_u64(w: &mut impl Write, n: u64) -> Result<()> {
    w.write_all(&n.to_le_bytes())?;
    Ok(())
}
fn get_u64(r: &mut impl Read) -> Result<u64> {
    let mut b = [0; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}
fn get_u32(r: &mut impl Read) -> Result<u32> {
    let mut b = [0; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn immutable_gene_metadata_is_not_rewritten_per_batch() {
        let dir = tempfile::tempdir().unwrap();
        let mut c =
            Checkpoint::create(&dir.path().join("counts"), Metadata::test_metadata()).unwrap();
        let metadata = std::fs::read(c.root.join("metadata.json")).unwrap();
        let layers = ["ambiguous", "spliced", "unspliced"]
            .into_iter()
            .map(|name| (name.to_string(), ndarray::Array2::zeros((2, 0))))
            .collect();
        c.commit_batch(0, 1, &layers, &[]).unwrap();
        assert_eq!(
            std::fs::read(c.root.join("metadata.json")).unwrap(),
            metadata
        );
        let manifest: serde_json::Value =
            serde_json::from_reader(std::fs::File::open(c.root.join("manifest.json")).unwrap())
                .unwrap();
        assert!(manifest.get("metadata").is_none());
    }

    #[test]
    fn resumes_inline_v1_metadata_and_migrates_to_v2() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path().join("counts");
        let metadata = Metadata::test_metadata();
        let store = Checkpoint::create(&root, metadata.clone()).unwrap();
        drop(store);

        let mut manifest: serde_json::Value =
            serde_json::from_reader(File::open(root.join("manifest.json")).unwrap()).unwrap();
        manifest["version"] = 1.into();
        manifest["metadata"] =
            serde_json::from_reader(File::open(root.join("metadata.json")).unwrap()).unwrap();
        serde_json::to_writer(File::create(root.join("manifest.json")).unwrap(), &manifest)
            .unwrap();
        fs::remove_file(root.join("metadata.json")).unwrap();

        let migrated = Checkpoint::resume(&root, &metadata).unwrap();
        assert_eq!(migrated.manifest.version, VERSION);
        assert_eq!(migrated.manifest.metadata, metadata);
        assert!(root.join("metadata.json").exists());
        let manifest: serde_json::Value =
            serde_json::from_reader(File::open(root.join("manifest.json")).unwrap()).unwrap();
        assert_eq!(manifest["version"], VERSION);
        assert!(manifest.get("metadata").is_none());
    }

    #[test]
    fn sparse_checkpoint_roundtrip_and_completion() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path().join("counts");
        let metadata = Metadata::test_metadata();
        let mut store = Checkpoint::create(&root, metadata.clone()).unwrap();
        let layers = [
            ("spliced".into(), array![[0, 70_000], [3, 0]]),
            ("unspliced".into(), array![[0, 2], [0, 0]]),
            ("ambiguous".into(), array![[1, 0], [0, 0]]),
        ]
        .into();
        store
            .commit_batch(0, 1234, &layers, &["|a".into(), "S|b".into()])
            .unwrap();
        assert!(!store.manifest.complete());
        store.finish_source(0).unwrap();
        assert!(store.manifest.complete());
        let batch = store.read_batch(0).unwrap();
        assert_eq!(batch.keys, ["|a", "S|b"]);
        assert_eq!(batch.layers[1][1], [(0, 70_000)]);
        assert_eq!(batch.layers[0][0], [(0, 1)]);
        assert_eq!(store.manifest.sources[0].offset, 1234);
        drop(store);
        let reopened = Checkpoint::resume(&root, &metadata).unwrap();
        assert_eq!(reopened.read_batch(0).unwrap().layers, batch.layers);
    }

    #[test]
    fn corrupt_shards_fail_without_removal() {
        let dir = tempfile::tempdir().unwrap();
        let mut store =
            Checkpoint::create(&dir.path().join("counts"), Metadata::test_metadata()).unwrap();
        let layers = [
            ("spliced".into(), ndarray::Array2::zeros((2, 0))),
            ("unspliced".into(), ndarray::Array2::zeros((2, 0))),
            ("ambiguous".into(), ndarray::Array2::zeros((2, 0))),
        ]
        .into();
        store.commit_batch(0, 0, &layers, &[]).unwrap();
        let path = store.root.join(&store.manifest.batches[0].file);
        std::fs::OpenOptions::new()
            .write(true)
            .open(&path)
            .unwrap()
            .set_len(12)
            .unwrap();
        assert!(store.read_batch(0).is_err());
        assert!(path.exists());
    }

    #[test]
    fn rejects_reuse_concurrent_writers_and_incompatible_resume() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path().join("counts");
        let metadata = Metadata::test_metadata();
        let store = Checkpoint::create(&root, metadata.clone()).unwrap();
        assert!(Checkpoint::create(&root, metadata.clone()).is_err());
        assert!(Checkpoint::resume(&root, &metadata).is_err());
        drop(store);
        let mut incompatible = metadata;
        incompatible.settings = "different tags".into();
        assert!(Checkpoint::resume(&root, &incompatible).is_err());
        assert!(root.join("manifest.json").exists());
    }
}
