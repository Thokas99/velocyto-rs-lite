#[cfg(feature = "bam")]
use crate::counter::CountProgress;
use anyhow::{Context, Result};
use serde_json::{json, Value};
#[cfg(feature = "bam")]
use std::time::Duration;
use std::{
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
    path::Path,
    time::{Instant, SystemTime, UNIX_EPOCH},
};

pub(crate) struct Telemetry {
    writer: BufWriter<File>,
    started: Instant,
    started_unix_ms: u128,
}

impl Telemetry {
    pub(crate) fn open(path: &Path) -> Result<Self> {
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .with_context(|| format!("Opening telemetry log {}", path.display()))?;
        Ok(Self {
            writer: BufWriter::new(file),
            started: Instant::now(),
            started_unix_ms: SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_millis(),
        })
    }

    #[cfg(feature = "bam")]
    pub(crate) fn phase(&mut self, name: &str, elapsed: Duration) -> Result<()> {
        self.write_event(json!({
            "event": "phase_complete",
            "phase": name,
            "phase_elapsed_seconds": elapsed.as_secs_f64(),
        }))
    }

    #[cfg(feature = "bam")]
    pub(crate) fn count_progress(
        &mut self,
        source: usize,
        checkpoint_commits: u64,
        progress: CountProgress,
    ) -> Result<()> {
        let records_per_second = if progress.elapsed_seconds > 0.0 {
            progress.records_processed as f64 / progress.elapsed_seconds
        } else {
            0.0
        };
        let event = if progress.batch_committed {
            json!({
                "event": "checkpoint_commit",
                "source_index": source,
                "checkpoint_commits": checkpoint_commits + 1,
                "records_processed": progress.records_processed,
                "completed_barcodes": progress.completed_barcodes,
                "completed_cells": progress.completed_cells,
                "cursor": progress.cursor,
                "batch_size": {
                    "raw_barcodes": progress.active_barcodes,
                    "cells": progress.active_cells,
                    "records": progress.active_records,
                },
                "estimated_batch_memory_bytes": {
                    "read_batch": progress.estimated_read_bytes,
                    "dense_scratch": progress.estimated_dense_bytes,
                },
                "records_per_second": records_per_second,
            })
        } else {
            json!({
                "event": "count_progress",
                "source_index": source,
                "records_processed": progress.records_processed,
                "completed_barcodes": progress.completed_barcodes,
                "active_batch_size": {
                    "raw_barcodes": progress.active_barcodes,
                    "cells": progress.active_cells,
                    "records": progress.active_records,
                },
                "estimated_memory_bytes": {
                    "read_batch": progress.estimated_read_bytes,
                    "dense_scratch": progress.estimated_dense_bytes,
                },
                "records_per_second": records_per_second,
            })
        };
        self.write_event(event)
    }

    #[cfg(feature = "bam")]
    pub(crate) fn source_complete(&mut self, source: usize, batches: usize) -> Result<()> {
        self.write_event(json!({
            "event": "source_complete",
            "source_index": source,
            "checkpoint_batches_total": batches,
        }))
    }

    #[cfg(feature = "bam")]
    pub(crate) fn complete(&mut self) -> Result<()> {
        self.write_event(json!({ "event": "run_complete" }))
    }

    pub(crate) fn error(&mut self, error: &anyhow::Error) -> Result<()> {
        self.write_event(json!({
            "event": "error",
            "error": format!("{error:#}"),
        }))
    }

    fn write_event(&mut self, mut event: Value) -> Result<()> {
        let event = event
            .as_object_mut()
            .context("telemetry event must be a JSON object")?;
        event.insert("run_started_unix_ms".into(), json!(self.started_unix_ms));
        event.insert(
            "elapsed_seconds".into(),
            json!(self.started.elapsed().as_secs_f64()),
        );
        let (rss_kib, peak_rss_kib) = process_memory_kib();
        event.insert("rss_kib".into(), json!(rss_kib));
        event.insert("peak_rss_kib".into(), json!(peak_rss_kib));
        serde_json::to_writer(&mut self.writer, &event)?;
        self.writer.write_all(b"\n")?;
        self.writer.flush()?;
        Ok(())
    }
}

fn process_memory_kib() -> (Option<u64>, Option<u64>) {
    let Ok(status) = std::fs::read_to_string("/proc/self/status") else {
        return (None, None);
    };
    let value = |name: &str| {
        status
            .lines()
            .find_map(|line| line.strip_prefix(name))?
            .split_whitespace()
            .next()?
            .parse()
            .ok()
    };
    (value("VmRSS:"), value("VmHWM:"))
}

#[cfg(all(test, feature = "bam"))]
mod tests {
    use super::*;

    #[test]
    fn appends_flushed_jsonl_phase_and_error_events_with_memory_fields() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("telemetry.jsonl");
        let mut telemetry = Telemetry::open(&path).unwrap();
        telemetry
            .phase("annotation", Duration::from_millis(20))
            .unwrap();
        telemetry
            .error(&anyhow::anyhow!("fixture failure"))
            .unwrap();
        drop(telemetry);

        let lines: Vec<Value> = std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .map(|line| serde_json::from_str(line).unwrap())
            .collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0]["phase"], "annotation");
        assert_eq!(lines[1]["event"], "error");
        for event in lines {
            assert!(event["elapsed_seconds"].is_number());
            assert!(event["rss_kib"].is_number() || event["rss_kib"].is_null());
            assert!(event["peak_rss_kib"].is_number() || event["peak_rss_kib"].is_null());
        }
    }
}
