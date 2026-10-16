#!/usr/bin/env python3
"""Aggregate independent fresh q88 lifecycle runs; retain operation scopes."""
import json
import statistics
import sys
from pathlib import Path

root = Path(sys.argv[1])
files = sorted(root.glob("run*/measurement.json"))
if not files:
    raise SystemExit("No completed lifecycle measurements")
runs = []
for path in files:
    data = json.loads(path.read_text())
    phases = {row["phase"]: row["seconds"] for row in data["workflow_measurements"]}
    if "full_test_wall_including_checks_and_cleanup" not in phases or not data["post_update_show"]:
        raise SystemExit(f"Incomplete lifecycle: {path}")
    phases["issue_three_credentials_total"] = sum(phases[k] for k in
        ("issue_uav1", "issue_uav2", "issue_uav_local"))
    phases["initial_show_generation_including_serialization"] = data["mean"]["generation_seconds"]
    phases["initial_show_wire_verification"] = data["mean"]["verification_seconds"]
    runs.append({"run": path.parent.name, "seconds": phases,
                 "initial_packet_bytes": data["measurements"][0]["packet_bytes"],
                 "updated_packet_bytes": data["post_update_show"][0]["packet_bytes"]})
keys = list(runs[0]["seconds"])
if any(set(run["seconds"]) != set(keys) for run in runs):
    raise SystemExit("Lifecycle phase sets differ")
stats = {key: {"mean_seconds": statistics.mean(run["seconds"][key] for run in runs),
               "min_seconds": min(run["seconds"][key] for run in runs),
               "max_seconds": max(run["seconds"][key] for run in runs)} for key in keys}
summary = {"profile": "production q88 compact CRT v3", "independent_runs": len(runs),
           "includes_fresh_setup_issue_revoke_update": True,
           "full_test_wall_includes_negative_checks_and_extra_sync_tests": True,
           "remote_sync_uses_local_authenticated_callback_not_network_latency": True,
           "runs": runs, "statistics": stats}
(root / "lifecycle_measurement.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(stats, indent=2))
