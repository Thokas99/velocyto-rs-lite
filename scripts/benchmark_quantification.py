#!/usr/bin/env python3
"""Measure one Linux process; retain reports and all user-supplied work files."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def disk_bytes(root):
    total = 0
    if root is not None and root.exists():
        for parent, _, files in os.walk(root):
            for name in files:
                try:
                    total += (Path(parent) / name).stat().st_size
                except FileNotFoundError:
                    pass
    return total


def measure(command, workdir=None, interval=0.05):
    start = time.perf_counter()
    process = subprocess.Popen(command)
    peak_disk = disk_bytes(workdir)
    last_io = {}
    next_disk = start + 1
    while True:
        exited = os.waitid(os.P_PID, process.pid, os.WEXITED | os.WNOHANG | os.WNOWAIT)
        try:
            fields = Path(f"/proc/{process.pid}/io").read_text().splitlines()
            last_io = {key: int(value) for key, value in (line.split(":", 1) for line in fields)}
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
        if exited is not None:
            _, status, usage = os.wait4(process.pid, 0)
            process.returncode = os.waitstatus_to_exitcode(status)
            break
        if time.perf_counter() >= next_disk:
            peak_disk = max(peak_disk, disk_bytes(workdir))
            next_disk += 1
        time.sleep(interval)
    final_disk = disk_bytes(workdir)
    return {
        "command": command,
        "exit_code": process.returncode,
        "wall_seconds": time.perf_counter() - start,
        "cpu_user_seconds": usage.ru_utime,
        "cpu_system_seconds": usage.ru_stime,
        "peak_rss_kib": usage.ru_maxrss,
        "storage_read_bytes": usage.ru_inblock * 512,
        "storage_write_bytes": usage.ru_oublock * 512,
        "process_io": last_io,
        "peak_sampled_disk_bytes": max(peak_disk, final_disk),
        "final_disk_bytes": final_disk,
        "scope": "Linux child process rusage; use Slurm/cgroup accounting for whole-job peak. /proc I/O includes the final zombie-state sample. Disk sampled every second; cache state uncontrolled.",
    }


def self_test():
    import tempfile
    with tempfile.TemporaryDirectory() as root:
        path = Path(root) / "retained.bin"
        report = measure([sys.executable, "-c",
                          "import pathlib,sys; pathlib.Path(sys.argv[1]).write_bytes(b'x'*65536); sys.exit(7)",
                          str(path)], Path(root), 0.001)
        assert report["exit_code"] == 7
        assert report["peak_rss_kib"] > 0
        assert report["wall_seconds"] > 0
        assert report["final_disk_bytes"] == 65536
        assert path.exists()
    print("benchmark self-check passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--workdir", type=Path, help="Directory to sample for intermediate/output disk size")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or args.report is None:
        parser.error("--report and a command after -- are required")
    args.report.parent.mkdir(parents=True, exist_ok=True)
    result = measure(command, args.workdir)
    # Even failed/OOM-killed children leave a report. Never delete their intermediates.
    args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in ("exit_code", "wall_seconds", "peak_rss_kib", "storage_read_bytes", "storage_write_bytes")}))
    return result["exit_code"] if result["exit_code"] >= 0 else 128 - result["exit_code"]


if __name__ == "__main__":
    sys.exit(main())
