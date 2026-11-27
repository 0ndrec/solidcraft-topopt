"""Runs the FPV frame study headless with the TopOpt plugin, as you would click it in SolidCraft:
design body, supports, loads, kept surfaces, Optimize, Create body.

    python -I run_study.py PLUGIN_ZIP WORKDIR [--solver]

WORKDIR gets plugins/, results/ (the plugin's results folder), limits.json and the designs:
study.sc (set up and optimized), optimized.sc (with the result as a solid body). With --solver,
TopOpt Solver (or TopOpt CUDA Solver) must already be watching WORKDIR/results: the plugin then
hands the job to it. frame.sc (from build_frame.py) must be next to this script.

Loads (see README.md for where they come from):
- case 1, full throttle: 4 x 4.12 N up on the motor seats; the frame held at the stack holes,
  where the battery and electronics are (their inertia is what holds the frame back);
- cases 2-5, a side hit on one arm: 10 N on the arm's end, in the frame's plane, square to the arm.
"""
import json
import math
import os
import subprocess
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
zip_path, work = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
use_solver = "--solver" in sys.argv
plugins, results = os.path.join(work, "plugins"), os.path.join(work, "results")
os.makedirs(results, exist_ok=True)
json.dump({"max_fuel": 10_000_000_000_000, "max_memory_mib": 2048}, open(os.path.join(work, "limits.json"), "w"))
cli = ["solidcraft-cli", "--plugin-limits", os.path.join(work, "limits.json"), "--plugins-dir", plugins, "--keep-plugin-approvals"]
dry = json.loads(subprocess.run(["solidcraft-cli", "--plugins-dir", plugins, "exec", "plugin.install", json.dumps({"path": zip_path, "dry_run": True})],
                                capture_output=True, text=True, check=True).stdout)
approval = dry["approval"]

THRUST = 4.12  # N per motor: 420 g, a 1404 motor's static thrust on a 3" prop at 4S
HIT = 10.0     # N, a side hit on one arm
R = 60 / math.sqrt(2)


def entry(name, params):
    return {"command": f"topopt.{name}", "params": {"approval": approval, "params": params}}


def change(path, value):
    key = "flag" if isinstance(value, bool) else "number"
    return entry("panel", {"panel": "topopt", "event": "change", "path": path, "row": None, "value": {key: value}})


def face(i, point):
    return {"type": "face", "body": "Frame", "index": i, "point": point}


seats = [face(14, [R, R, 6.5]), face(16, [-R, R, 6.5]), face(18, [-R, -R, 6.5]), face(20, [R, -R, 6.5])]
# The round ends of the arms (FR, FL, RL, RR), picked at their outermost point.
ends = [(13, (1, 1)), (15, (-1, 1)), (17, (-1, -1)), (21, (1, -1))]
stack_holes = [face(i, None) for i in range(42, 46)]

install = {"command": "plugin.install", "params": {"path": zip_path, "approval": approval, "replace": True, "approve": {
    "command_groups": ["file"], "fuel": 10_000_000_000_000, "memory_mib": 1024,
    "folders": [{"name": "results", "path": results, "read": True, "write": True}]}}}
steps = [
    install,
    entry("design-body", {"body": [{"type": "body", "name": "Frame"}], "material": "petg", "volume_percent": 30,
                          "element_size": 0.8, "max_elements": 200000, "min_member": 2.0}),
    entry("add-support", {"faces": stack_holes, "fix_x": True, "fix_y": True, "fix_z": True}),
    entry("add-load", {"faces": seats, "mode": "vector", "fx": 0, "fy": 0, "fz": 4 * THRUST, "magnitude": 0, "case": 1}),
]
for case, (i, (sx, sy)) in enumerate(ends, start=2):
    d = math.hypot(sx, sy)
    point = [sx * (R + 11 / d), sy * (R + 11 / d), 2.5]
    # Square to the arm, in the plane: the arm runs along (sx, sy); turn it a quarter turn.
    f = [-sy / d * HIT, sx / d * HIT, 0.0]
    steps.append(entry("add-load", {"faces": [face(i, point)], "mode": "vector", "fx": f[0], "fy": f[1], "fz": 0, "magnitude": 0, "case": case}))
steps += [
    # The motor seats are mating faces: kept exactly as drawn (the holes stay exact through the
    # pinning of the body's faces).
    entry("keep-face", {"faces": seats, "thickness": 0}),
    change("/settings/use_solver", use_solver),
    change("/settings/max_iterations", 200),
    change("/settings/time_limit", 3600),
    entry("optimize", {}),
    {"command": "storage.get", "params": {"plugin": "topopt"}},
]
# Faces picked by index only get their point from model.faces.
faces_script, faces_out = os.path.join(work, "faces.json"), os.path.join(work, "faces-out.json")
json.dump([{"command": "model.faces", "params": {"body": "Frame"}}], open(faces_script, "w"))
subprocess.run(["solidcraft-cli", "run", faces_script, "--in", os.path.join(here, "frame.sc"), "--results", "--json-out", faces_out], capture_output=True, check=True)
faces = {f["index"]: f["point"] for f in json.load(open(faces_out, encoding="utf-8"))[0]["faces"]}
for s in steps:
    for f in s.get("params", {}).get("params", {}).get("faces", []):
        if f["point"] is None:
            f["point"] = faces[f["index"]]


def run(script, name, design_in, design_out):
    path, out = os.path.join(work, f"{name}.json"), os.path.join(work, f"{name}-out.json")
    json.dump(script, open(path, "w"), indent=1)
    p = subprocess.run(cli + ["run", path, "--in", design_in, "--results", "--timings", "--json-out", out, "--save", design_out],
                       capture_output=True, text=True, encoding="utf-8")
    if p.returncode != 0:
        raise SystemExit(p.stderr[-3000:])
    return json.load(open(out, encoding="utf-8"))


t0 = time.time()
d1 = run(steps, "study", os.path.join(here, "frame.sc"), os.path.join(work, "study.sc"))
print("optimize:", json.dumps(d1["results"][-2]["result"])[:300])
if use_solver:
    job = d1["results"][-2]["result"]["sent_to_solver"]
    rf = os.path.join(results, "topopt", "solver", "results", job + ".json")
    while not os.path.exists(rf):
        if time.time() - t0 > 3600:
            raise SystemExit("no solver result after an hour")
        time.sleep(2)
    print(f"solver finished after {time.time() - t0:.0f} s")
d2 = run([install, entry("refresh", {}), entry("create-body", {}), {"command": "storage.get", "params": {"plugin": "topopt"}}],
         "create", os.path.join(work, "study.sc"), os.path.join(work, "optimized.sc"))
body = d2["results"][2]["result"]
result = d2["results"][3]["document"]["result"]
json.dump({"body": body, "result": result}, open(os.path.join(work, "result.json"), "w"), indent=1)
print("create-body:", body["body"], {k: body["measure"][k] for k in ["closed", "volume_mm3", "faces"]})
print("result:", json.dumps({k: v for k, v in result.items() if k != "cases"}, ensure_ascii=False)[:1500])
for c in result.get("cases", []):
    print("  case", c)
