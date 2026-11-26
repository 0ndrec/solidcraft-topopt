"""Writes a solidcraft-cli `run` script that installs the TopOpt zip into a scratch plugins folder
and runs a whole study on the test bracket (box 80×40×20 with a Ø10 hole):

    python -I scripts/topopt_headless.py OUT.json ZIP APPROVAL RESULTS_DIR [max_elements] [iterations]

Then: solidcraft-cli --plugin-limits limits.json --plugins-dir DIR --keep-plugin-approvals \
        run OUT.json --in part.sc --results --timings --json-out out.json
"""
import json
import os
import sys

out, zip_path, approval, results = sys.argv[1:5]
max_elements = int(sys.argv[5]) if len(sys.argv) > 5 else 8000
iterations = int(sys.argv[6]) if len(sys.argv) > 6 else 40


def entry(name, params):
    return {"command": f"topopt.{name}", "params": {"approval": approval, "params": params}}


def face(i, point):
    return {"type": "face", "body": "Bracket", "index": i, "point": point}


script = [
    {"command": "plugin.install", "params": {"path": zip_path, "approval": approval, "replace": True, "approve": {
        "command_groups": ["file"], "fuel": 10_000_000_000_000, "memory_mib": 1024,
        "folders": [{"name": "results", "path": os.path.normpath(results), "read": True, "write": True}]}}},
    entry("design-body", {"body": [{"type": "body", "name": "Bracket"}], "material": "aluminium-6061", "volume_percent": 30,
                    "element_size": 0, "max_elements": max_elements, "min_member": 0}),
    entry("add-support", {"faces": [face(2, [0.0, 32.36, 7.64])], "fix_x": True, "fix_y": True, "fix_z": True}),
    entry("add-load", {"faces": [face(6, [61.3, 24.83, 16.45])], "mode": "vector", "fx": 0, "fy": 0, "fz": -2000, "magnitude": 0, "case": 1}),
    # A panel field change, as the host sends it.
    entry("panel", {"panel": "topopt", "event": "change", "path": "/mesh/triangles", "row": None, "value": {"number": 30000}}),
    entry("panel", {"panel": "topopt", "event": "change", "path": "/settings/max_iterations", "row": None, "value": {"number": iterations}}),
    entry("refresh", {}),
    entry("optimize", {}),
    entry("create-body", {}),
]
with open(out, "w", encoding="utf-8") as f:
    json.dump(script, f, indent=1)
