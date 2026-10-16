#!/usr/bin/env python3
"""Summarize actual serialized Show packet measurements from the C demo."""
import json
import re
import statistics
import sys
from pathlib import Path

results = Path(sys.argv[1])
log = (results / "show.log").read_text(errors="replace")
rows = [json.loads(line.split("[q88 full packet measurement] ", 1)[1])
        for line in log.splitlines() if "[q88 full packet measurement] " in line]
if not rows:
    raise SystemExit("Missing full packet measurements")
worker = (results / "worker.log").read_text(errors="replace")
if (results / "last_show_packet.bin").stat().st_size != rows[-1]["packet_bytes"]:
    raise SystemExit("Serialized Show artifact size disagrees with log")
if (results / "last_accumulator_envelope.bin").stat().st_size != rows[-1]["accumulator_bytes"]:
    raise SystemExit("Accumulator artifact size disagrees with log")
summary = {
    "profile": "production complete q88 CRT compact v3",
    "server": "autodl-pqa11",
    "generation_includes_serialization": True,
    "excludes_initialization_issue_and_revocation": True,
    "measurements": rows,
    "accumulator_worker_generation_s": [float(v) for v in
        re.findall(r"prove request total: ([0-9.]+) s", worker)],
    "accumulator_worker_verification_s": [float(v) for v in
        re.findall(r"verify request total: ([0-9.]+) s", worker)],
}
summary["workflow_measurements"] = [
    json.loads(line.split("[q88 workflow measurement] ", 1)[1])
    for line in log.splitlines() if "[q88 workflow measurement] " in line]
summary["post_update_show"] = [
    json.loads(line.split("[q88 post-update packet measurement] ", 1)[1])
    for line in log.splitlines() if "[q88 post-update packet measurement] " in line]
summary["mean"] = {key: statistics.mean(row[key] for row in rows)
                   for key in rows[0] if key != "run"}
summary["range"] = {key: [min(row[key] for row in rows), max(row[key] for row in rows)]
                    for key in rows[0] if key != "run"}
summary["target_met_all_runs"] = all(row["generation_seconds"] < 3 and
                                         row["verification_seconds"] < 1.5 for row in rows)
(results / "measurement.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
