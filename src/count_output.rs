//! Bounded HDF5 slice writes from immutable sparse count checkpoints.
use crate::count_checkpoint::{Checkpoint, SparseBatch, SparseRow};
use anyhow::{ensure, Context, Result};
use hdf5::{types::VarLenUnicode, Dataset, Group, H5Type, Location};
use std::{collections::BTreeMap, fs, path::Path, str::FromStr};

fn string_attr(location: &Location, name: &str, value: &str) -> Result<()> {
    location
        .new_attr::<VarLenUnicode>()
        .shape(())
        .create(name)?
        .write_scalar(&VarLenUnicode::from_str(value)?)?;
    Ok(())
}
fn encoding(location: &Location, kind: &str, version: &str) -> Result<()> {
    string_attr(location, "encoding-type", kind)?;
    string_attr(location, "encoding-version", version)
}
fn strings(group: &Group, name: &str, values: &[String], encoded: bool) -> Result<()> {
    let data: Vec<_> = values
        .iter()
        .map(|s| VarLenUnicode::from_str(s))
        .collect::<std::result::Result<_, _>>()?;
    let ds = group
        .new_dataset::<VarLenUnicode>()
        .shape(data.len())
        .create(name)?;
    ds.write_raw(&data)?;
    if encoded {
        encoding(&ds, "string-array", "0.2.0")?;
    }
    Ok(())
}
fn vector<T: H5Type>(group: &Group, name: &str, len: usize) -> Result<Dataset> {
    let builder = group.new_dataset::<T>().shape(len);
    Ok(if len == 0 {
        builder.create(name)?
    } else {
        builder.chunk(len.min(65_536)).deflate(2).create(name)?
    })
}
fn dataframe(group: &Group, index: &str, columns: &[&str]) -> Result<()> {
    encoding(group, "dataframe", "0.2.0")?;
    string_attr(group, "_index", index)?;
    let cols: Vec<_> = columns
        .iter()
        .map(|s| VarLenUnicode::from_str(s))
        .collect::<std::result::Result<_, _>>()?;
    group
        .new_attr::<VarLenUnicode>()
        .shape(cols.len())
        .create("column-order")?
        .write_raw(&cols)?;
    Ok(())
}
fn metadata(
    file: &hdf5::File,
    checkpoint: &Checkpoint,
    ids: &[String],
    samples: &[String],
    anndata: bool,
) -> Result<()> {
    let genes = &checkpoint.manifest.metadata.genes;
    let obs = file.create_group(if anndata { "obs" } else { "col_attrs" })?;
    let has_samples = samples.iter().any(|s| !s.is_empty());
    if anndata {
        dataframe(
            &obs,
            "_index",
            if has_samples { &["SampleID"] } else { &[] },
        )?;
    }
    strings(
        &obs,
        if anndata { "_index" } else { "CellID" },
        ids,
        anndata,
    )?;
    if has_samples {
        strings(&obs, "SampleID", samples, anndata)?;
    }
    let var = file.create_group(if anndata { "var" } else { "row_attrs" })?;
    if anndata {
        dataframe(
            &var,
            "Accession",
            &["Gene", "Chromosome", "Strand", "Start", "End"],
        )?;
    }
    for (name, values) in [
        (
            "Accession",
            genes.iter().map(|g| g.id.clone()).collect::<Vec<_>>(),
        ),
        ("Gene", genes.iter().map(|g| g.name.clone()).collect()),
        (
            "Chromosome",
            genes.iter().map(|g| g.chrom.clone()).collect(),
        ),
        ("Strand", genes.iter().map(|g| g.strand.clone()).collect()),
    ] {
        strings(&var, name, &values, anndata)?;
    }
    for (name, values) in [
        ("Start", genes.iter().map(|g| g.start).collect::<Vec<_>>()),
        ("End", genes.iter().map(|g| g.end).collect()),
    ] {
        let ds = vector::<i64>(&var, name, values.len())?;
        ds.write_raw(&values)?;
        if anndata {
            encoding(&ds, "array", "0.2.0")?;
        }
    }
    for name in if anndata {
        &["uns", "obsm", "varm", "obsp", "varp"][..]
    } else {
        &["row_graphs", "col_graphs"][..]
    } {
        let group = file.create_group(name)?;
        if anndata {
            encoding(&group, "dict", "0.1.0")?;
        }
    }
    if !anndata {
        let attrs = file.create_group("attrs")?;
        strings(&attrs, "LOOM_SPEC_VERSION", &["3.0.0".into()], false)?;
    }
    Ok(())
}
fn total_rows(batch: &SparseBatch) -> Result<Vec<Vec<(u32, f32)>>> {
    (0..batch.keys.len())
        .map(|c| {
            let mut row = BTreeMap::<u32, u64>::new();
            for layer in &batch.layers {
                for &(g, v) in &layer[c] {
                    let entry = row.entry(g).or_default();
                    *entry = entry
                        .checked_add(v as u64)
                        .context("total count overflow")?;
                }
            }
            Ok(row.into_iter().map(|(g, v)| (g, v as f32)).collect())
        })
        .collect()
}
struct CsrTarget {
    data: Dataset,
    indices: Dataset,
    indptr: Dataset,
    offset: usize,
}
impl CsrTarget {
    fn create<T: H5Type>(group: &Group, cells: usize, genes: usize, nnz: u64) -> Result<Self> {
        encoding(group, "csr_matrix", "0.1.0")?;
        group
            .new_attr::<i64>()
            .shape(2)
            .create("shape")?
            .write_raw(&[i64::try_from(cells)?, i64::try_from(genes)?])?;
        let len = usize::try_from(nnz)?;
        ensure!(nnz <= i64::MAX as u64, "CSR offset exceeds int64");
        Ok(Self {
            data: vector::<T>(group, "data", len)?,
            indices: vector::<i32>(group, "indices", len)?,
            indptr: vector::<i64>(
                group,
                "indptr",
                cells.checked_add(1).context("row count overflow")?,
            )?,
            offset: 0,
        })
    }
    fn write<T: H5Type + Copy>(
        &mut self,
        rows: &[Vec<(u32, T)>],
        cell_offset: usize,
    ) -> Result<()> {
        let mut ptr = Vec::with_capacity(rows.len() + 1);
        let mut indices = vec![];
        let mut data = vec![];
        ptr.push(i64::try_from(self.offset)?);
        for row in rows {
            for &(g, v) in row {
                indices.push(i32::try_from(g)?);
                data.push(v);
            }
            ptr.push(i64::try_from(
                self.offset
                    .checked_add(data.len())
                    .context("CSR overflow")?,
            )?);
        }
        if !data.is_empty() {
            let end = self
                .offset
                .checked_add(data.len())
                .context("CSR overflow")?;
            self.data.write_slice(&data[..], self.offset..end)?;
            self.indices.write_slice(&indices[..], self.offset..end)?;
            self.offset = end;
        }
        self.indptr
            .write_slice(&ptr[..], cell_offset..cell_offset + ptr.len())?;
        Ok(())
    }
}
fn atomic_output(path: &Path, write: impl FnOnce(&Path) -> Result<()>) -> Result<()> {
    let name = path
        .file_name()
        .context("output needs a filename")?
        .to_string_lossy();
    let partial = path.with_file_name(format!("{name}.partial-{:016x}", rand::random::<u64>()));
    write(&partial).with_context(|| {
        format!(
            "output failed; checkpoints and partial file {} retained",
            partial.display()
        )
    })?;
    let f = fs::File::open(&partial)?;
    f.sync_all()?;
    // Verify the backend can reopen the completed container before publication.
    hdf5::File::open(&partial)?.close()?;
    fs::rename(&partial, path)?;
    fs::File::open(
        path.parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new(".")),
    )?
    .sync_all()?;
    Ok(())
}
pub fn write_h5ad(checkpoint: &Checkpoint, path: &Path) -> Result<()> {
    ensure!(
        checkpoint.manifest.complete(),
        "counting is incomplete; resume counting before export"
    );
    let cells = checkpoint.manifest.rows()?;
    let meta = &checkpoint.manifest.metadata;
    let mut ids = Vec::with_capacity(cells);
    let mut samples = Vec::with_capacity(cells);
    let mut x_nnz = 0u64;
    let mut layer_nnz = vec![0u64; meta.layers.len()];
    for i in 0..checkpoint.manifest.batches.len() {
        let batch = checkpoint.read_batch(i)?;
        for key in &batch.keys {
            let (id, sample) = meta.cell_id(key);
            ids.push(id);
            samples.push(sample);
        }
        for row in total_rows(&batch)? {
            x_nnz = x_nnz
                .checked_add(row.len() as u64)
                .context("X size overflow")?;
        }
        for (l, total) in layer_nnz.iter_mut().enumerate() {
            *total = total
                .checked_add(checkpoint.manifest.batches[i].nnz[l])
                .context("layer size overflow")?;
        }
    }
    atomic_output(path, |partial| {
        let file = hdf5::File::create_excl(partial)?;
        encoding(&file, "anndata", "0.1.0")?;
        metadata(&file, checkpoint, &ids, &samples, true)?;
        let mut x =
            CsrTarget::create::<f32>(&file.create_group("X")?, cells, meta.genes.len(), x_nnz)?;
        let group = file.create_group("layers")?;
        encoding(&group, "dict", "0.1.0")?;
        let narrow = meta.dtype == "uint16";
        let mut targets = vec![];
        for (l, name) in meta.layers.iter().enumerate() {
            let sg = group.create_group(name)?;
            targets.push(if narrow {
                CsrTarget::create::<u16>(&sg, cells, meta.genes.len(), layer_nnz[l])?
            } else {
                CsrTarget::create::<u32>(&sg, cells, meta.genes.len(), layer_nnz[l])?
            });
        }
        let mut cell_offset = 0;
        let mut overflow = 0u64;
        for i in 0..checkpoint.manifest.batches.len() {
            let batch = checkpoint.read_batch(i)?;
            x.write(&total_rows(&batch)?, cell_offset)?;
            for (l, target) in targets.iter_mut().enumerate() {
                if narrow {
                    let rows: Vec<Vec<_>> = batch.layers[l]
                        .iter()
                        .map(|row| {
                            row.iter()
                                .map(|&(g, v)| {
                                    if v > u16::MAX as u32 {
                                        overflow += 1;
                                    }
                                    (g, v.min(u16::MAX as u32) as u16)
                                })
                                .collect()
                        })
                        .collect();
                    target.write(&rows, cell_offset)?;
                } else {
                    target.write(&batch.layers[l], cell_offset)?;
                }
            }
            cell_offset += batch.keys.len();
        }
        if overflow > 0 {
            log::warn!("{overflow} layer values saturated by --dtype uint16; checkpoints retain exact uint32 counts");
        }
        ensure!(x.offset == usize::try_from(x_nnz)?, "X length mismatch");
        for (l, target) in targets.iter().enumerate() {
            ensure!(
                target.offset == usize::try_from(layer_nnz[l])?,
                "layer length mismatch"
            );
        }
        file.flush()?;
        drop(targets);
        drop(group);
        drop(x);
        file.close()?;
        Ok(())
    })
}
pub fn write_loom(checkpoint: &Checkpoint, path: &Path) -> Result<()> {
    ensure!(
        checkpoint.manifest.complete(),
        "counting is incomplete; resume counting before export"
    );
    let cells = checkpoint.manifest.rows()?;
    let meta = &checkpoint.manifest.metadata;
    let genes = meta.genes.len();
    let mut ids = vec![];
    let mut samples = vec![];
    for i in 0..checkpoint.manifest.batches.len() {
        for key in checkpoint.read_batch(i)?.keys {
            let (id, sample) = meta.cell_id(&key);
            ids.push(id);
            samples.push(sample);
        }
    }
    atomic_output(path, |partial| {
        let file = hdf5::File::create_excl(partial)?;
        metadata(&file, checkpoint, &ids, &samples, false)?;
        let total_dataset = matrix::<f32>(&file, "matrix", genes, cells)?;
        let group = file.create_group("layers")?;
        let narrow = meta.dtype == "uint16";
        let mut targets = vec![];
        for name in &meta.layers {
            targets.push(if narrow {
                matrix::<u16>(&group, name, genes, cells)?
            } else {
                matrix::<u32>(&group, name, genes, cells)?
            });
        }
        let mut cell_offset = 0;
        let mut overflow = 0u64;
        for i in 0..checkpoint.manifest.batches.len() {
            let batch = checkpoint.read_batch(i)?;
            let width = batch.keys.len();
            if width == 0 {
                continue;
            }
            for start in (0..genes).step_by(256) {
                let end = (start + 256).min(genes);
                let mut total = ndarray::Array2::<f32>::zeros((end - start, width));
                for (l, target) in targets.iter().enumerate() {
                    let mut tile = ndarray::Array2::<u32>::zeros((end - start, width));
                    for (c, row) in batch.layers[l].iter().enumerate() {
                        let first = row.partition_point(|&(g, _)| (g as usize) < start);
                        let last = row.partition_point(|&(g, _)| (g as usize) < end);
                        for &(g, v) in &row[first..last] {
                            tile[(g as usize - start, c)] = v;
                        }
                    }
                    for ((r, c), &v) in tile.indexed_iter() {
                        total[(r, c)] += v as f32;
                    }
                    if narrow {
                        let narrowed = tile.mapv(|v| {
                            if v > u16::MAX as u32 {
                                overflow += 1;
                            }
                            v.min(u16::MAX as u32) as u16
                        });
                        target.write_slice(
                            narrowed.view(),
                            (start..end, cell_offset..cell_offset + width),
                        )?;
                    } else {
                        target.write_slice(
                            tile.view(),
                            (start..end, cell_offset..cell_offset + width),
                        )?;
                    }
                }
                total_dataset
                    .write_slice(total.view(), (start..end, cell_offset..cell_offset + width))?;
            }
            cell_offset += width;
        }
        if overflow > 0 {
            log::warn!("{overflow} layer values saturated by --dtype uint16; checkpoints retain exact uint32 counts");
        }
        file.flush()?;
        drop(targets);
        drop(group);
        drop(total_dataset);
        file.close()?;
        Ok(())
    })
}
fn matrix<T: H5Type>(group: &Group, name: &str, genes: usize, cells: usize) -> Result<Dataset> {
    let builder = group.new_dataset::<T>().shape((genes, cells));
    Ok(if genes == 0 || cells == 0 {
        builder.create(name)?
    } else {
        builder
            .chunk((genes.min(256), cells.min(100)))
            .deflate(2)
            .create(name)?
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::count_checkpoint::{Checkpoint, Metadata};
    use ndarray::array;

    #[test]
    fn failed_output_retains_partial_and_previous_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("existing.h5ad");
        std::fs::write(&path, "old").unwrap();
        let result = atomic_output(&path, |partial| {
            std::fs::write(partial, "incomplete")?;
            anyhow::bail!("simulated OOM or disk failure");
        });
        assert!(result.is_err());
        assert_eq!(std::fs::read_to_string(path).unwrap(), "old");
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 2);
    }

    #[test]
    fn streams_all_formats_and_preserves_exact_layers() {
        let dir = tempfile::tempdir().unwrap();
        let mut checkpoint =
            Checkpoint::create(&dir.path().join("counts"), Metadata::test_metadata()).unwrap();
        for (key, n) in [("|a", 70_000), ("S|b", 3)] {
            let layers = [
                ("spliced".into(), array![[n], [0]]),
                ("unspliced".into(), array![[0], [2]]),
                ("ambiguous".into(), array![[0], [0]]),
            ]
            .into();
            checkpoint
                .commit_batch(0, 0, &layers, &[key.into()])
                .unwrap();
        }
        checkpoint.finish_source(0).unwrap();
        let h5ad = dir.path().join("test.h5ad");
        let loom = dir.path().join("test.loom");
        write_h5ad(&checkpoint, &h5ad).unwrap();
        write_loom(&checkpoint, &loom).unwrap();
        let file = hdf5::File::open(&h5ad).unwrap();
        assert_eq!(
            file.dataset("layers/spliced/data")
                .unwrap()
                .read_raw::<u32>()
                .unwrap(),
            [70_000, 3]
        );
        assert_eq!(
            file.dataset("X/indptr").unwrap().read_raw::<i64>().unwrap(),
            [0, 2, 4]
        );
        assert_eq!(
            file.dataset("X/data").unwrap().read_raw::<f32>().unwrap(),
            [70_000., 2., 3., 2.]
        );
        let ids = file
            .dataset("obs/_index")
            .unwrap()
            .read_raw::<hdf5::types::VarLenUnicode>()
            .unwrap();
        assert_eq!(
            ids.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            ["test:a", "test_S:b"]
        );
        let file = hdf5::File::open(&loom).unwrap();
        assert_eq!(
            file.dataset("layers/spliced")
                .unwrap()
                .read_raw::<u32>()
                .unwrap(),
            [70_000, 3, 0, 0]
        );
        assert_eq!(
            file.dataset("matrix").unwrap().read_raw::<f32>().unwrap(),
            [70_000., 3., 2., 2.]
        );
    }

    #[test]
    fn incomplete_or_corrupt_checkpoint_cannot_replace_existing_output() {
        let dir = tempfile::tempdir().unwrap();
        let c = Checkpoint::create(&dir.path().join("counts"), Metadata::test_metadata()).unwrap();
        let out = dir.path().join("old.h5ad");
        std::fs::write(&out, "previous valid output").unwrap();
        assert!(write_h5ad(&c, &out).is_err());
        assert_eq!(
            std::fs::read_to_string(&out).unwrap(),
            "previous valid output"
        );
    }

    #[test]
    fn empty_checkpoint_writes_valid_zero_cell_shapes() {
        let dir = tempfile::tempdir().unwrap();
        let mut c =
            Checkpoint::create(&dir.path().join("counts"), Metadata::test_metadata()).unwrap();
        c.finish_source(0).unwrap();
        let out = dir.path().join("empty.h5ad");
        write_h5ad(&c, &out).unwrap();
        let file = hdf5::File::open(&out).unwrap();
        assert_eq!(
            file.dataset("X/indptr").unwrap().read_raw::<i64>().unwrap(),
            [0]
        );
        assert_eq!(file.dataset("X/data").unwrap().size(), 0);
    }
}
