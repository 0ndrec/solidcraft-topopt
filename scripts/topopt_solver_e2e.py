"""End-to-end test of TopOpt with TopOpt Solver through solidcraft-cli (headless).

    python -I scripts/topopt_solver_e2e.py WORKDIR ZIP [max_elements] [iterations] [use_solver 0|1]

WORKDIR holds part.sc (the test bracket), limits.json, plugins/ and results/; TopOpt Solver must
already watch WORKDIR/results when use_solver is 1. Runs: install, design body, support, load,
options, Optimize; waits for the solver's result; then Refresh and Create body. Prints timings
and the inserted body's measure.
"""
import json
import os
import subprocess
import sys
import time

work, zip_path = sys.argv[1], sys.argv[2]
max_elements = int(sys.argv[3]) if len(sys.argv) > 3 else 30000
iterations = int(sys.argv[4]) if len(sys.argv) > 4 else 60
use_solver = (sys.argv[5] if len(sys.argv) > 5 else "1") == "1"
results = os.path.join(work, "results")
plugins = os.path.join(work, "plugins")
cli = ["solidcraft-cli", "--plugin-limits", os.path.join(work, "limits.json"), "--plugins-dir", plugins, "--keep-plugin-approvals"]

dry = json.loads(subprocess.run(["solidcraft-cli", "--plugins-dir", plugins, "exec", "plugin.install", json.dumps({"path": zip_path, "dry_run": True})],
                                capture_output=True, text=True, check=True).stdout)
approval = dry["approval"]


def entry(name, params):
    return {"command": f"topopt.{name}", "params": {"approval": approval, "params": params}}


def change(path, value):
    key = "flag" if isinstance(value, bool) else "text" if isinstance(value, str) else "number"
    return entry("panel", {"panel": "topopt", "event": "change", "path": path, "row": None, "value": {key: value}})


def face(i, point):
    return {"type": "face", "body": "Bracket", "index": i, "point": point}


def run(script, name, design_in, design_out=None):
    path = os.path.join(work, f"{name}.json")
    out = os.path.join(work, f"{name}-out.json")
    json.dump(script, open(path, "w"))
    args = cli + ["run", path, "--in", design_in, "--results", "--timings", "--json-out", out]
    if design_out:
        args += ["--save", design_out]
    p = subprocess.run(args, capture_output=True, text=True, encoding="utf-8")
    if p.returncode != 0:
        print(p.stderr[-2000:])
        sys.exit(1)
    return json.load(open(out, encoding="utf-8"))


part1 = [
    {"command": "plugin.install", "params": {"path": zip_path, "approval": approval, "replace": True, "approve": {
        "command_groups": ["file"], "fuel": 10_000_000_000_000, "memory_mib": 1024,
        "folders": [{"name": "results", "path": os.path.normpath(results), "read": True, "write": True}]}}},
    entry("design-body", {"body": [{"type": "body", "name": "Bracket"}], "material": "aluminium-6061", "volume_percent": 30,
                          "element_size": 0, "max_elements": max_elements, "min_member": 0}),
    entry("add-support", {"faces": [face(2, [0.0, 32.36, 7.64])], "fix_x": True, "fix_y": True, "fix_z": True}),
    entry("add-load", {"faces": [face(6, [61.3, 24.83, 16.45])], "mode": "vector", "fx": 0, "fy": 0, "fz": -2000, "magnitude": 0, "case": 1}),
    change("/settings/max_iterations", iterations),
    change("/settings/use_solver", use_solver),
    # Extra panel options as JSON, e.g. TOPOPT_E2E_OPTIONS='{"/mesh/upsample": 1}'.
    *[change(k, v) for k, v in json.loads(os.environ.get("TOPOPT_E2E_OPTIONS", "{}")).items()],
    # Kept surfaces as JSON [[face index, [x, y, z]], ...], e.g. TOPOPT_E2E_KEEP='[[1, [71.2, 40, 8.8]]]'.
    *[entry("keep-face", {"faces": [face(i, pt)], "thickness": 0}) for i, pt in json.loads(os.environ.get("TOPOPT_E2E_KEEP", "[]"))],
    entry("optimize", {}),
]
t0 = time.time()
d1 = run(part1, "e2e-1", os.path.join(work, "part.sc"), os.path.join(work, "e2e-1.sc"))
sent = d1["results"][-1]["result"]
print("optimize:", json.dumps(sent)[:300], f"({d1['timings'][-1]['ms'] / 1000:.1f} s)")

if use_solver:
    if "sent_to_solver" not in sent:
        print("FAIL: the job was not handed to the solver")
        sys.exit(1)
    rf = os.path.join(results, "topopt", "solver", "results", sent["sent_to_solver"] + ".json")
    while not os.path.exists(rf):
        if time.time() - t0 > 900:
            print("FAIL: no solver result after 15 min")
            sys.exit(1)
        time.sleep(1)
    print(f"solver finished {time.time() - t0:.1f} s after the job was sent")

part2 = [part1[0], entry("refresh", {}), entry("create-body", {}), {"command": "storage.get", "params": {"plugin": "topopt"}}]
d2 = run(part2, "e2e-2", os.path.join(work, "e2e-1.sc"), os.path.join(work, "e2e-2.sc"))
body = d2["results"][2]["result"]
m = body["measure"]
doc = d2["results"][3]["document"]
r = doc["result"]
print("create-body:", body["body"], "closed", m["closed"], "volume", round(m["volume_mm3"], 1), "faces", m["faces"])
print("result:", {k: r[k] for k in ["outcome", "iterations", "seconds", "mass_g", "volume_fraction", "triangles", "inserted"]}, r["warnings"])
print("pending after:", doc.get("pending"), "error:", doc.get("last_error"))
