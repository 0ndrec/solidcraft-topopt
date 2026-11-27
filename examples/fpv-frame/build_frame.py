"""Writes frame.json: a solidcraft-cli script that builds the solid blank of a 3-inch FPV frame
(the design body for TopOpt), and runs it.

    python build_frame.py            # writes frame.json, builds frame.sc, prints the faces

Sizes (mm): 120 mm wheelbase (motor to motor, diagonally); a 5 mm plate, 46 mm square core and
22 mm wide arms with round ends; raised motor seats 22 mm across and 1.5 mm high; 12 x 12 mm
M2 motor holes and a 5.5 mm shaft relief; a 20 x 20 mm M2 stack pattern.
"""
import json
import math
import os
import subprocess

WB, T, ARM, CORE, PAD_D, PAD_H = 120.0, 5.0, 22.0, 46.0, 22.0, 1.5
M2, SHAFT, MOTOR_SQ, STACK_SQ = 2.2, 5.5, 12.0, 20.0
R = WB / 2 / math.sqrt(2)  # motor centre offset along X and Y
MOTORS = [(sx * R, sy * R) for sx, sy in [(1, 1), (-1, 1), (-1, -1), (1, -1)]]
NAMES = ["FR", "FL", "RL", "RR"]  # front is +X

steps = [
    {"command": "sketch.create", "params": {"plane": "XY", "name": "Outline"}},
    {"command": "sketch.rectangle.center", "params": {"center": [0, 0], "corner": [CORE / 2, CORE / 2]}},
    {"command": "sketch.slot.center_to_center", "params": {"p0": [-R, -R], "p1": [R, R], "width": ARM}},
    {"command": "sketch.slot.center_to_center", "params": {"p0": [-R, R], "p1": [R, -R], "width": ARM}},
    {"command": "sketch.finish", "params": {}},
    {"command": "solid.extrude", "params": {"distance": T, "name": "Plate"}},
    {"command": "body.rename", "params": {"body": "Body1", "name": "Frame"}},
]
for (x, y), n in zip(MOTORS, NAMES):
    steps.append({"command": "solid.cylinder", "params": {"diameter": PAD_D, "height": PAD_H, "base": [x, y, T], "operation": "join", "name": f"Seat {n}"}})
through = T + PAD_H + 2
for (x, y), n in zip(MOTORS, NAMES):
    steps.append({"command": "solid.cylinder", "params": {"diameter": SHAFT, "height": through, "base": [x, y, -1], "operation": "cut", "name": f"Shaft {n}"}})
    for k, (dx, dy) in enumerate([(1, 1), (-1, 1), (-1, -1), (1, -1)]):
        steps.append({"command": "solid.cylinder", "params": {"diameter": M2, "height": through, "base": [x + dx * MOTOR_SQ / 2, y + dy * MOTOR_SQ / 2, -1],
                                                                 "operation": "cut", "name": f"Motor hole {n}{k + 1}"}})
for k, (dx, dy) in enumerate([(1, 1), (-1, 1), (-1, -1), (1, -1)]):
    steps.append({"command": "solid.cylinder", "params": {"diameter": M2, "height": through, "base": [dx * STACK_SQ / 2, dy * STACK_SQ / 2, -1],
                                                             "operation": "cut", "name": f"Stack hole {k + 1}"}})
steps += [{"command": "inspect.measure", "params": {"bodies": ["Frame"]}}, {"command": "model.faces", "params": {"body": "Frame"}}]
json.dump(steps, open("frame.json", "w"), indent=1)

p = subprocess.run(["solidcraft-cli", "run", "frame.json", "--save", "frame.sc", "--results", "--json-out", "frame-out.json"], capture_output=True, text=True)
if p.returncode:
    raise SystemExit(p.stderr[-3000:])
out = json.load(open("frame-out.json", encoding="utf-8"))
m = out[-2]["bodies"][0]
print({k: m[k] for k in ["volume_mm3", "area_mm2", "bbox", "faces", "closed"]})
faces = out[-1]["faces"]
os.remove("frame-out.json")
for f in faces:
    print(f["index"], f["name"], round(f["area"], 2), [round(c, 2) for c in f["point"]], f.get("plane_normal"), f.get("cylinder", {}).get("radius") if f.get("cylinder") else "")
