"""Independent STL check with trimesh: welded corners, watertight, consistent winding, positive
volume, single body, Euler characteristic.

    tools/reference/.venv/Scripts/python tools/reference/check_stl.py FILE.stl
"""
import json
import sys

import trimesh

m = trimesh.load(sys.argv[1], force="mesh", process=True)  # process=True welds equal corners
bodies = m.split(only_watertight=False)
print(json.dumps({
    "vertices": len(m.vertices),
    "faces": len(m.faces),
    "is_watertight": bool(m.is_watertight),
    "is_winding_consistent": bool(m.is_winding_consistent),
    "is_volume": bool(m.is_volume),
    "volume_mm3": round(float(m.volume), 1),
    "bodies": len(bodies),
    "euler_number": int(m.euler_number),
    "bbox": [list(map(float, m.bounds[0])), list(map(float, m.bounds[1]))],
}, indent=1))
