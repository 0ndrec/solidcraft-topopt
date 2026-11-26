# TopOpt Solver

`tools/topopt-solver` runs topology optimization on every CPU core. It is three things in one package:

- **a command line for pipelines**: a study described in JSON goes in; a closed mesh (PLY/STL), the density field and a summary come out. Progress is printed as JSON lines on stderr, the summary as JSON on stdout, and the exit code tells how it ended;
- **a Rust library**: `run()` with a progress closure and a cancel flag, which writes files only when asked;
- **the background solver of the SolidCraft plugin** (`watch`): it takes the plugin's jobs from its `results` folder and solves them natively instead of in the single-threaded plugin sandbox.

Every mode runs the same code ([topopt-core](topopt-core.md) → [topopt-mesh](topopt-mesh.md), with the problem assembled by [topopt-study](topopt-study.md)), so a study solved here gives the same result as in the plugin.

Contents:
1. [Build](#1-build)
2. [Quick start](#2-quick-start)
3. [Commands](#3-commands)
4. [The spec format](#4-the-spec-format)
5. [Outputs](#5-outputs)
6. [Progress, exit codes, cancelling](#6-progress-exit-codes-cancelling)
7. [Pipeline examples](#7-pipeline-examples)
8. [Rust API](#8-rust-api)
9. [Watch mode: the plugin protocol](#9-watch-mode-the-plugin-protocol)
10. [Performance and grid size](#10-performance-and-grid-size)
11. [Troubleshooting](#11-troubleshooting)

---

## 1. Build

```sh
scripts/build-solver.sh
```

The script builds the release binary and copies it with its launcher `start-topopt-solver.cmd` into `build/solver/`. By hand:

```sh
cargo +stable-x86_64-pc-windows-gnu build --release -p topopt-solver   # -> target/release/topopt-solver.exe
```

The result is one self-contained executable (~3.7 MB). Nothing in the solver is Windows-specific: on Linux or macOS, build it with the native toolchain.

## 2. Quick start

```sh
# 1. An example study: a 120×40×40 mm cantilever, clamped at x = 0, 1000 N down on the far bottom edge
topopt-solver example > cantilever.json

# 2. Check it: the grid, the nodes the supports and loads caught, a time estimate
topopt-solver check cantilever.json

# 3. Solve it: out/cantilever.ply, out/cantilever.density.json, out/cantilever.result.json
topopt-solver solve cantilever.json --out out
```

`solve` prints its progress to stderr (`--progress text`, the default):

```
TopOpt Solver 0.1.0: cantilever on 6 threads
  0%      0.0 s  voxelizing (36×12×12, 3.34 mm)
  5%      0.1 s  iteration 1 — compliance …, β 1
 ...
 90%      9.6 s  making the mesh
100%     10.0 s  done
done: iteration limit in 10.0 s, 150 iterations
  wrote out/cantilever.ply
  ...
```

and the summary to stdout ([§5.3](#53-summary-resultjson)). More specs are in [`examples/`](../examples).

## 3. Commands

```
topopt-solver solve INPUT.json [options]   solve one study (a spec, or a TopOpt plugin job)
topopt-solver check INPUT.json             build the problem only; print a JSON report
topopt-solver remesh INPUT.json --density FILE --mesh FILE
                                           mesh a run's densities again with INPUT's mesh settings
topopt-solver example                      print an example spec
topopt-solver watch [FOLDER] [--once]      serve the TopOpt plugin (its results folder)
topopt-solver version
topopt-solver --help
```

`topopt-solver FOLDER` and no arguments at all mean `watch`.

`INPUT` is a JSON file, or `-` for standard input. Relative paths inside a spec (an STL, say) resolve from **the spec file's folder**, or from the current folder when the spec comes from standard input.

The input kind is detected from its content: JSON with the keys `problem` and `study` is a job written by the plugin (already voxelized, [§9](#9-watch-mode-the-plugin-protocol)); anything else is a spec ([§4](#4-the-spec-format)).

### 3.1. `solve`

| Option | Effect |
|---|---|
| `--out DIR` | Write all three files into `DIR`: `NAME.ply`, `NAME.density.json`, `NAME.result.json`. This is the default, with the current folder, when no output option is given |
| `--mesh FILE` | Write the mesh; `.stl` gives binary STL, anything else binary PLY |
| `--density FILE` | Write the density field (JSON, [§5.2](#52-densities-densityjson)) |
| `--result FILE` | Write the summary (JSON); it is printed to stdout anyway |
| `--name NAME` | The `NAME` of `--out` files. Default: the spec's `name`, the body's name for a plugin job, else `topopt` |
| `--no-files` | Write no files; only print the summary |
| `--progress MODE` | `json` (one JSON event per line), `text` (readable lines, the default) or `none`. Always on **stderr** |
| `--time-limit S` | Stop optimizing after S seconds. The mesh is still made from the design reached; the summary says `outcome: "time limit"` and warns |
| `--max-iterations N` | At most N optimization iterations |
| `--cancel-file FILE` | Stop as soon as `FILE` exists (exit code 3) |
| `--threads N` | Worker threads. Default: every logical core |

Explicit `--mesh`, `--density` and `--result` combine with `--out`: the explicit paths win and the other files go into `--out`. `--mesh x.stl` alone writes only the mesh.

`--time-limit` and `--max-iterations` override `settings.time_limit` and `settings.max_iterations` of the input, so a pipeline can tune them without editing the spec.

### 3.2. `check`

`check` builds the problem exactly as `solve` would: it reads the STL files, voxelizes the domain, picks the support and load nodes and checks that the problem can be solved. It does not optimize. The report goes to stdout:

```json
{
  "name": "cantilever",
  "grid": [24, 8, 8],
  "h": 5.03,
  "elements": 1536,
  "design_elements": 1464,
  "solid_elements": 72,
  "fixed_nodes": 81,
  "load_cases": 1,
  "loaded_nodes": 9,
  "volfrac": 0.2656,
  "filter_radius": 7.55,
  "estimate_seconds": 6.6,
  "threads": 6,
  "warnings": []
}
```

| Field | Meaning |
|---|---|
| `grid`, `h` | Grid size in elements and element edge, mm |
| `elements` | Elements inside the domain. Of these, `design_elements` are optimized and `solid_elements` always stay solid (pads under supports and loads, `keep_in`) |
| `fixed_nodes`, `loaded_nodes`, `load_cases` | Nodes caught by the supports and loads, and the number of load cases |
| `volfrac` | The material fraction **of the design elements**: `settings.volfrac` (a fraction of the whole domain) corrected for the pads and clamped to 0.05…0.95 |
| `filter_radius` | Density filter radius, mm: it sets the thinnest member |
| `estimate_seconds` | Rough run time on this machine. Calibrated before the multigrid thresholds were lowered, so it usually overestimates now |

A broken study (a region that catches no node, no supports, an open domain…) gives exit code 2 and the reason on stderr, which makes `check` a cheap validation step in CI before a long run.

### 3.3. `remesh`

```
topopt-solver remesh INPUT.json --density RUN.density.json --mesh OUT.ply|stl [--trim PERCENT] [--threads N]
```

Makes the mesh again from the densities of an earlier `solve` of the same study, with the `mesh` settings of `INPUT` (triangles, smoother, upsampling, sharp edges…), without optimizing again: tune the mesh in seconds. The densities must come from the same domain and grid settings (the grid is checked). `--trim PERCENT` trims thin dead weight costing at most that much stiffness ([§4.10](#410-trimming-thin-dead-weight)), overriding `mesh.trim_percent`. A JSON report goes to stdout:

```json
{"triangles": 30000, "vertices": 14998, "volume_mm3": 80043.9, "genus": 31,
 "pinned_vertices": 12879, "edge_vertices": 663, "corner_vertices": 6,
 "sharp_edges": true, "smooth_iterations": 20,
 "trim": {"candidates": 21, "removed_pieces": 6, "removed_elements": 26, "removed_fraction": 0.0029,
          "compliance_before": 105.03, "compliance_after": 105.52},
 "warnings": []}
```

`trim` is `null` without trimming.

Exit codes as for `solve` (1 when the densities do not fit the input or the mesh cannot be made).

### 3.4. `example`

Prints a complete spec with every field at its default. It is the best template: delete any field to get its default.

### 3.5. `watch`

```
topopt-solver watch [FOLDER] [--once] [--threads N]
```

Watches the folder the TopOpt plugin uses as `results` and solves its jobs as they come ([§9](#9-watch-mode-the-plugin-protocol)). `--once` solves the jobs already there and exits (for scripts and tests).

## 4. The spec format

A spec describes a study without a CAD application. The design domain is a box, an STL file or a triangle mesh inside the JSON. Supports and loads act on **regions of space** (boxes, spheres, planes) instead of face numbers, so specs are easy to generate.

Units: **mm, N, MPa, g/cm³**. No axis is special: give forces in the direction you need.

### 4.1. Structure

```jsonc
{
  "name": "bracket",                                   // optional: for reports and file names
  "domain": { "stl": { "path": "bracket.stl" } },     // required: the design domain
  "keep_out": [ { "box": { "min": [0,0,0], "max": [10,10,40] } } ],     // solids that stay empty
  "keep_in":  [ { "sphere": { "center": [60,20,20], "radius": 5 } } ],  // regions that stay solid
  "material": { "preset": "aluminium-6061" },
  "settings": { "volfrac": 0.3, "max_elements": 60000 },
  "mesh":     { "triangles": 30000 },
  "supports": [ { "region": { "plane": { "point": [0,0,0], "normal": [1,0,0] } } } ],
  "loads":    [ { "region": { "sphere": { "center": [100,20,0], "radius": 6 } }, "force": [0,0,-1000] } ]
}
```

Unknown fields at the top level and in `domain`, `keep_out`, `keep_in`, regions, `material`, `supports` and `loads` are **errors**, so a typo never passes silently. `settings` and `mesh` ignore unknown fields and fill in missing ones with defaults.

### 4.2. Solids (`domain`, `keep_out`)

| Shape | JSON | Notes |
|---|---|---|
| Box | `{"box": {"min": [x,y,z], "max": [x,y,z]}}` | Axis-aligned |
| STL | `{"stl": {"path": "part.stl", "scale": 1}}` | Binary or ASCII. `path` is relative to the spec. `scale` multiplies the file's units into mm: `1000` for a file in metres, `25.4` for inches. Equal corners are welded |
| Mesh | `{"mesh": {"positions": [[x,y,z], …], "triangles": [[a,b,c], …]}}` | Zero-based indices |

A solid must be **closed**. Triangle winding does not matter: inside and outside come from the parity of ray crossings. An open mesh is reported ("the domain mesh is not closed…", "rays crossed the domain an odd number of times") and may voxelize wrongly near its holes.

`keep_out` is subtracted from the domain: material never goes there (room for a bolt, a cable, a neighbouring part).

### 4.3. Regions (`supports`, `loads`, `keep_in`)

| Region | JSON | Selects |
|---|---|---|
| Box | `{"box": {"min": [..], "max": [..], "tolerance": t}}` | Points inside the box grown by `t`. It may be flat: `min` and `max` equal along one axis (a face) or two (an edge) |
| Sphere | `{"sphere": {"center": [..], "radius": r, "tolerance": t}}` | Points within `r + t` of the centre |
| Plane | `{"plane": {"point": [..], "normal": [..], "tolerance": t, "within": [[..],[..]]}}` | The slab of half-thickness `t` around the plane; `within` (optional) limits it to a box `[min, max]` |

`tolerance` is optional and defaults to **0.6 element edges**, so a region drawn exactly on a face of the part catches the grid nodes at that face wherever the grid falls.

How a region becomes a boundary condition:
1. Take the grid nodes that touch material and lie in the region.
2. If some of them are **on the surface** of the material, keep only those; otherwise keep all (a region buried inside the part).
3. No node: the study is refused, e.g. `load 2 : its region holds no node of the part (grid step 3.34 mm)`. Grow the region or its `tolerance`, or use smaller elements.
4. Cells within `settings.pad` of the chosen nodes (default 1.5 elements) become **solid pads**: the optimizer never removes the material right under a support or load.

`keep_in` makes solid the elements whose centres lie in the region (no tolerance), e.g. a sleeve around a hole.

### 4.4. Supports

```json
{ "region": { … }, "dirs": [true, true, true] }
```

`dirs` are the held directions X, Y, Z (default: all). `[false, false, true]` holds only along Z and lets the part slide in the plane. The supports together must stop rigid-body motion, or the solver cannot converge.

### 4.5. Loads

```json
{ "region": { … }, "force": [0, 0, -1000], "case": 1 }
```

- `force` is the **total** force, N, spread evenly over the region's nodes.
- `case` is the load case (default 1). Loads with the same case act together. Different cases are separate situations: the optimizer minimizes the sum of their compliances, so the part must withstand each.

### 4.6. Material

```json
{ "preset": "steel" }                         // a preset
{ "preset": "pla", "e": 3000 }                // a preset with some values replaced
{ "e": 70000, "nu": 0.33, "density": 2.7 }    // your own
```

| Preset | E, MPa | ν | ρ, g/cm³ |
|---|---|---|---|
| `aluminium-6061` | 68 900 | 0.33 | 2.70 |
| `aluminium-7075` | 71 700 | 0.33 | 2.81 |
| `steel` | 200 000 | 0.29 | 7.85 |
| `stainless-316` | 193 000 | 0.30 | 8.00 |
| `titanium` | 113 800 | 0.342 | 4.43 |
| `pla` | 3 500 | 0.36 | 1.24 |
| `petg` | 2 100 | 0.38 | 1.27 |
| `abs` | 2 200 | 0.35 | 1.05 |
| `nylon` | 1 700 | 0.39 | 1.01 |

Without `material`, `aluminium-6061` is used. Requirements: E > 0, 0 ≤ ν < 0.5, ρ > 0. The optimal shape does not depend on E and depends little on ν; the material mostly changes the compliance and mass in the summary.

### 4.7. `settings`

The same fields and defaults as the plugin:

| Field | Default | Meaning |
|---|---|---|
| `volfrac` | 0.3 | Material to keep, as a fraction of the domain's volume (0.05…0.95) |
| `element_size` | 0 | Element edge, mm; 0 = from `max_elements` |
| `max_elements` | 40 000 | Element budget over the domain's bounding box: the main knob for detail and time |
| `min_member` | 0 | Thinnest member, mm; 0 = about three elements. Filter radius = max(1.5 h, min_member / 2) |
| `pad` | 0 | Thickness of the solid pads under supports and loads, mm; 0 = 1.5 elements |
| `max_iterations` | 150 | Iteration limit |
| `time_limit` | 480 | Optimization time limit, s |
| `penal` | 3 | SIMP penalty: higher gives fewer grey densities |
| `beta_max` | 16 | Final sharpness of the Heaviside projection |
| `move_limit` | 0.2 | Largest density change of an element per iteration |
| `use_solver` | true | Not used by specs (a plugin flag) |

Elements are cubes. Grid dimensions are rounded up to multiples of 4 (of 2 when small) so multigrid can coarsen, which can leave the element count well below the budget; `check` shows the real count.

### 4.8. `mesh`

| Field | Default | Meaning |
|---|---|---|
| `triangles` | 30 000 | Triangle budget after decimation |
| `upsample` | 2 | Sample the density field this many times finer (1…3, cubic interpolation) before making the surface: 2 gives round members instead of faceted ones, for about 8× the meshing work |
| `smoothing` | 20 | Smoothing passes |
| `keep_faces` | true | Vertices near the domain's surface are pulled onto it, so mounting faces stay exact |
| `remove_islands` | true | Remove pieces not connected to a support or load |
| `fill_cavities` | true | Fill closed internal voids (they cannot be manufactured, and CAD dislikes them) |
| `iso` | 0.5 | Density threshold of the surface |
| `smoother` | `"bilateral"` | `"bilateral"` (bilateral normal filtering: the most accurate surface, keeps creases), `"taubin"` (fastest) or `"hc"` (HC Laplacian, little shrinkage) |
| `preserve_volume` | true | Offset the free surface along its normals so smoothing does not change the volume |
| `sharp_edges` | true | With `keep_faces`: keep the domain's sharp edges and corners exactly (crease vertices inserted on them, corners fixed) |
| `fillets` | true | Round the joints where the optimized material meets the pads under supports and loads and `keep_in` regions with fillets ([topopt-mesh.md](topopt-mesh.md#rounded-joints-blend)) |
| `fillet_radius` | 0 | Fillet radius, mm; 0 = 1.5 elements. Larger radii round more but can close small openings near the joints |
| `trim_percent` | 0 | Trim thin dead weight as long as the stiffness drops by at most this many percent; 0 = off ([§4.10](#410-trimming-thin-dead-weight)) |
| `convert_solid` | true | Not used by specs (a plugin flag) |

The surface is made from the **filtered** densities (before the Heaviside projection), sampled `upsample` times finer: the projected densities are nearly 0 or 1 per element, so a surface through them follows the voxels, while the filtered field varies smoothly and puts the same surface (the projection's threshold is the iso value) between the element centres. On the bridge of `examples/`, cross-sections of the struts became rounder (fit to an ellipse 0.956 → 0.984) and the surface smoother (mean angle between neighbouring triangles 15.7° → 10.4° at the same triangle count); meshing takes ~20 s instead of ~1 s natively.

### 4.9. Full example: a bracket from STL with a bolt clearance

```json
{
  "name": "bracket",
  "domain": { "stl": { "path": "parts/bracket_envelope.stl", "scale": 1 } },
  "keep_out": [ { "stl": { "path": "parts/bolt_clearance.stl" } } ],
  "material": { "preset": "aluminium-7075" },
  "settings": { "volfrac": 0.25, "max_elements": 80000, "min_member": 3, "time_limit": 300 },
  "mesh": { "triangles": 40000 },
  "supports": [
    { "region": { "plane": { "point": [0, 0, 0], "normal": [1, 0, 0], "within": [[0, 0, 0], [0, 60, 15]] } } }
  ],
  "loads": [
    { "region": { "sphere": { "center": [140, 30, 40], "radius": 8 } }, "force": [0, 0, -2500], "case": 1 },
    { "region": { "sphere": { "center": [140, 30, 40], "radius": 8 } }, "force": [0, 1200, 0],  "case": 2 }
  ]
}
```

Two load cases: the part must carry the vertical and the side force, but not at the same time.

### 4.10. Trimming thin dead weight

Results often keep thin pieces that look useless: spikes, hooks and strands hanging from a member ("stalactites"). With `trim_percent` above 0 (or `remesh --trim`), after optimizing:
1. the part is analysed as it will be made: the design thresholded at `iso`, solid or empty;
2. candidates are its thin parts — what a morphological opening with a 3×3×3 cube removes (pieces up to about two elements thick) — grouped into connected pieces, never the pads under supports and loads;
3. each piece, smallest first, is removed on trial and the compliance recomputed; it goes if it costs less than half the budget, and as long as all removals together stay within `trim_percent`. A piece whose removal leaves part of the structure floating is kept.

The summary's `warnings` say what was trimmed (`trimmed 5 thin pieces (0.5 % of the material) for 0.78 % of the stiffness`) and `trimmed` gives the share of material.

Expect little: such pieces mostly **do** carry load. On the cantilever of `examples/`, 14 thin pieces (0.9 % of the material) cost 1.4 % of the stiffness — more than their share of the mass. On the bridge, the "stalactites" under the deck are short ribs stiffening the thin loaded deck (a pad the optimizer cannot thicken): with a 6 mm pad instead they vanish, but the deck then takes most of the material and the bridge is 6.6× less stiff. Trimming is a cosmetic choice with a known price, so it is off by default.

## 5. Outputs

### 5.1. Mesh (`.ply` / `.stl`)

Closed, 2-manifold, outward-oriented, without self-intersections (checked by `topopt-mesh` on every run, not assumed). Millimetres, in the spec's coordinates. PLY is binary little-endian: `vertex x y z` (float) and `face vertex_indices` (uchar + int). STL is binary, unindexed.

SolidCraft (`closed: true`, `mesh.convert` to a solid), trimesh (`is_watertight`), Blender and slicers accept these meshes.

### 5.2. Densities (`density.json`)

```jsonc
{
  "grid":   { "nx": 36, "ny": 12, "nz": 12, "h": 3.34, "origin": [x0, y0, z0] },
  "values": [ … ],   // element density 0…1, nx·ny·nz numbers
  "pads":   [ … ],   // true: the element was a solid pad
  "inside": [ … ],   // true: inside the domain (false: outside it or in a keep_out)
  "filtered": [ … ]  // the filtered densities before projection (pads 1, outside 0): the mesh is made from these
}
```

Element order: index `i + nx·(j + ny·k)`; the element's centre is `origin + h·(i+0.5, j+0.5, k+0.5)`. `values` are the physical densities after filtering and projection, i.e. what the analysis used; `filtered` (empty in files from older versions, which then mesh from `values`) is the smooth field the surface is made from. Trimmed pieces are 0 in both. Re-mesh with other settings from this file (as the plugin's Re-mesh does), plot sections of it or hand the field to another tool.

### 5.3. Summary (`result.json`)

```json
{
  "outcome": "iteration limit",
  "iterations": 150,
  "compliance": 109.79,
  "volume_fraction": 0.2395,
  "mass_g": 124.1,
  "mesh_volume_mm3": 45977.5,
  "triangles": 18000,
  "genus": 18,
  "elements": 5184,
  "grid": [36, 12, 12],
  "h": 3.344,
  "seconds": 10.0,
  "stl": "out/cantilever.ply",
  "inserted": false,
  "warnings": [],
  "trimmed": 0.0,
  "solver": "CPU, 6 threads"
}
```

| Field | Meaning |
|---|---|
| `outcome` | `converged`; or `iteration limit` / `time limit`: stopped by a limit, the mesh made from the last design |
| `compliance` | Work of the external forces, N·mm, summed over load cases. Lower is stiffer; compare it only between variants of one study |
| `volume_fraction` | Mesh volume / domain volume |
| `mass_g` | From the mesh volume and the material's density |
| `genus` | Number of through-holes in the part: a rough measure of shape complexity |
| `elements`, `grid`, `h` | The analysis grid |
| `seconds` | Total time: voxelization, optimization, mesh |
| `stl` | Where the mesh was written (a historical name: it may be PLY). Empty with `--no-files` |
| `inserted` | Used by the plugin; always `false` here |
| `trimmed` | Share of the material trimmed as dead weight (`mesh.trim_percent`), 0 when off |
| `solver` | What solved the linear systems: `CPU, 6 threads`, or a GPU's name with `topopt-cuda-solver` |
| `warnings` | Open domain, too few elements across, time limit reached, sharp edges that had to be rounded… **Show these to the pipeline's user** |

## 6. Progress, exit codes, cancelling

### 6.1. Streams

- **stdout**: only the final summary (JSON) on success, or the `check` report, or the `example` spec. Pipe it straight into `jq` or `json.loads`.
- **stderr**: progress and error messages.

### 6.2. `--progress json`

One event per line, always in this order: one `started`, any number of `progress`, then `finished` **or** `failed`.

```json
{"event":"started","name":"beam","threads":6,"version":"0.1.0"}
{"event":"progress","fraction":0.0,"message":"voxelizing (36×12×12, 3.34 mm)","elapsed":0.0}
{"event":"progress","fraction":0.0557,"message":"iteration 1 — compliance …, β 1","elapsed":0.07}
{"event":"progress","fraction":0.9,"message":"making the mesh","elapsed":9.6}
{"event":"progress","fraction":1.0,"message":"done","elapsed":10.0}
{"event":"finished","result":{…summary…},"files":["out/beam.ply","out/beam.density.json","out/beam.result.json"]}
```

`fraction` rises monotonically from 0 to 1: 0–0.05 voxelization, 0.05–0.9 iterations (relative to `max_iterations`, so it jumps when the study converges early), 0.9–1 the mesh. `elapsed` is seconds since the start. `message` is for people and may change: do not parse it.

On failure: `{"event":"failed","message":"load 1 : its region holds no node of the part (grid step 3.34 mm)"}`.

### 6.3. Exit codes

| Code | When |
|---|---|
| 0 | Done: summary on stdout, files written |
| 1 | The run failed: no convergence, empty mesh, a file could not be written |
| 2 | Bad arguments or input: missing file, not JSON, unknown field, a region without nodes (in `check`) |
| 3 | Cancelled through `--cancel-file` |

`solve` reports study errors found while it runs (a region without nodes…) with code 1; use `check` to validate input.

### 6.4. Cancelling

- `--cancel-file PATH`: the solver looks for the file every 0.2 s. When it appears, the current iteration finishes, no files are written, and the exit code is 3. Use it from orchestrators that cannot send signals.
- Ctrl+C or killing the process is safe too: files are written to a temporary name and renamed, so no half-written mesh is left.
- `--time-limit` is the soft stop: you get a result, just a less converged one.

## 7. Pipeline examples

### 7.1. PowerShell: a batch of parts

```powershell
# One study per STL in parts\, from a template
$template = Get-Content template.json -Raw | ConvertFrom-Json
foreach ($stl in Get-ChildItem parts\*.stl) {
    $spec = $template.PSObject.Copy()
    $spec.name = $stl.BaseName
    $spec.domain = @{ stl = @{ path = $stl.FullName } }
    $spec | ConvertTo-Json -Depth 20 | Set-Content "jobs\$($stl.BaseName).json"
    topopt-solver solve "jobs\$($stl.BaseName).json" --out results --progress none
    if ($LASTEXITCODE -ne 0) { Write-Warning "$($stl.Name): exit code $LASTEXITCODE" }
}
```

### 7.2. Python: run with progress and read the result

```python
import json, subprocess, sys

def optimize(spec: dict, out_dir: str, on_progress=None) -> dict:
    """Solves a spec and returns the summary. Raises RuntimeError with the solver's message."""
    proc = subprocess.Popen(
        ["topopt-solver", "solve", "-", "--out", out_dir, "--progress", "json"],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8")
    proc.stdin.write(json.dumps(spec))
    proc.stdin.close()
    error = None
    for line in proc.stderr:                      # events as the run goes
        event = json.loads(line)
        if event["event"] == "progress" and on_progress:
            on_progress(event["fraction"], event["message"])
        elif event["event"] == "failed":
            error = event["message"]
    result = proc.stdout.read()
    if proc.wait() != 0:
        raise RuntimeError(error or f"topopt-solver exit code {proc.returncode}")
    return json.loads(result)

spec = json.loads(subprocess.check_output(["topopt-solver", "example"]))
spec["settings"]["volfrac"] = 0.2
summary = optimize(spec, "out", lambda f, m: print(f"{f:5.0%} {m}", file=sys.stderr))
print(summary["mass_g"], "g;", summary["warnings"])
```

A spec read from standard input resolves relative paths from the process's current folder: pass absolute paths or set `cwd=`.

A parameter sweep (say `volfrac` from 0.15 to 0.4) is a loop around `optimize()`; compare `compliance` and `mass_g` in the summaries.

### 7.3. Python: check a result

```python
import trimesh
m = trimesh.load("out/cantilever.ply")
assert m.is_watertight and m.is_winding_consistent
print(m.volume, m.bounds)
```

### 7.4. CI

```yaml
- run: topopt-solver check specs/bracket.json          # quick validation, under a second
- run: topopt-solver solve specs/bracket.json --out artifacts --progress json --time-limit 600
```

### 7.5. A job prepared by the plugin

A job written by the plugin (`results/topopt/solver/jobs/<id>.part0…`, concatenated into one JSON file) can be solved on its own:

```sh
topopt-solver solve job.json --out out
```

This moves a heavy run to another machine: SolidCraft has already voxelized the part, so that machine needs only the executable.

## 8. Rust API

```toml
[dependencies]
topopt-solver = { path = "…/tools/topopt-solver" }   # library and binary in one package
```

### 8.1. Essentials

```rust
use std::sync::{Arc, atomic::AtomicBool};
use topopt_solver::{Event, Input, RunOptions, run, check, spec::Spec};

let spec: Spec = serde_json::from_str(&std::fs::read_to_string("bracket.json")?)?;
let input = Input::Spec { spec, base: "specs".into() };   // base: where STL paths resolve from

// Only check the study:
let report = check(&input)?;                 // CheckReport: grid, nodes, time estimate
println!("{} elements, ~{:.0} s", report.elements, report.estimate_seconds);

// Solve in memory: no files, with progress and a cancel flag
let cancel = Arc::new(AtomicBool::new(false));
let opts = RunOptions { cancel: Some(cancel.clone()), time_limit: Some(120.0), ..Default::default() };
let out = run(&input, &opts, &mut |e: &Event| {
    if let Event::Progress { fraction, message, .. } = e {
        eprintln!("{:3.0}% {message}", fraction * 100.0);
    }
})?;
let mesh = out.mesh;           // topopt_study::Mesh: positions (mm), triangles
let density = out.density;     // the density field
let summary = out.result;      // RunResult, as in result.json
```

### 8.2. Types

| Item | Purpose |
|---|---|
| `Input::Spec { spec, base }` | A study from a spec; `base` is the folder relative paths resolve from |
| `Input::Job(Box<Job>)` | A job prepared by the plugin |
| `Input::from_file(path)`, `Input::from_json(bytes, base)` | Read either kind, detected from the content |
| `RunOptions` | `mesh_out`, `density_out`, `result_out` (`Option<PathBuf>`); `time_limit`, `max_iterations`; `cancel: Option<Arc<AtomicBool>>`; `solver: Option<SolverFactory>` (`topopt_cuda::factory(…)` solves on a GPU) |
| `serve_with`, `remesh_with` | `serve` and `remesh` with a `SolverFactory` |
| `cli::main_with(App, args)` | The whole command line as a function: `topopt-cuda-solver` is this with a GPU factory |
| `run(&Input, &RunOptions, &mut dyn FnMut(&Event)) -> Result<RunOutput, String>` | A whole run. The error `"cancelled"` means it was cancelled |
| `RunOutput` | `result: RunResult`, `mesh: Mesh`, `density: Density`, `files: Vec<PathBuf>` |
| `Event` | `Started`, `Progress`, `Finished`, `Failed`; its `Serialize` is the `--progress json` format |
| `prepare(&Input, on_event) -> Result<Job, String>` | Only build the job (voxelize a spec) |
| `remesh(&Input, &Density) -> Result<(Mesh, RemeshReport), String>` | Mesh saved densities again with the input's mesh settings |
| `check(&Input) -> Result<CheckReport, String>` | Build and report |
| `mesh_bytes(&Mesh, path, comment)` | PLY or STL by extension |
| `write_atomic(path, bytes)` | Write through a temporary file and a rename |
| `threads()` | Worker threads |
| `Folder`, `Heartbeat`, `run_job`, `serve` | Watch mode ([§9](#9-watch-mode-the-plugin-protocol)) |
| `spec` | Re-export of `topopt_study::spec`: `Spec`, `Shape`, `Region`, `SupportSpec`, `LoadSpec`, `MaterialSpec`, `box_mesh`, `MeshCaster` |

The progress closure runs on the thread that called `run`. The computation inside is parallel (rayon), but events arrive one at a time. Set the thread count with `RAYON_NUM_THREADS` **before the first** run (the pool is created once per process).

### 8.3. A spec from code

```rust
use topopt_solver::spec::{Region, Shape, Spec, SupportSpec, LoadSpec};

let mut spec = Spec::example();
spec.domain = Shape::Mesh { positions, triangles };        // your own mesh, no files
spec.supports = vec![SupportSpec {
    region: Region::Box { min: [0.0, 0.0, 0.0], max: [0.0, 40.0, 40.0], tolerance: None },
    dirs: [true; 3],
}];
spec.loads = vec![LoadSpec {
    region: Region::Sphere { center: [120.0, 20.0, 0.0], radius: 5.0, tolerance: None },
    force: [0.0, 0.0, -500.0],
    case: 1,
}];
spec.settings.max_elements = 60_000;
```

One level down are `topopt_study::prepare_spec` and `topopt_study::solve` ([topopt-study](topopt-study.md)), or the core `topopt_core::Optimizer` itself ([topopt-core](topopt-core.md)) if you build the element grid and boundary conditions yourself.

## 9. Watch mode: the plugin protocol

The SolidCraft sandbox gives a plugin one thread, so the plugin hands heavy runs to the solver. All they share is the plugin's `results` folder (granted at install, Preferences › Plugins › TopOpt).

| File (inside `results`) | Written by | Content |
|---|---|---|
| `topopt/solver/alive.json` | solver, every second | Heartbeat: `protocol`, `version`, `time` (unix, s), `threads`, `job`, `progress`, `message` |
| `topopt/solver/jobs/<id>.partN` | plugin | The job (`Job`) as JSON, in parts of at most 1 000 000 bytes |
| `topopt/solver/jobs/<id>.json` | plugin | Header `{protocol, id, parts, bytes}` |
| `topopt/solver/jobs/<id>.cancel` | plugin | Asks to cancel |
| `topopt/<body>.ply` | solver | The mesh (path from `job.mesh_path`) |
| `topopt/solver/results/<id>.density.partN` | solver | Densities in parts (for Re-mesh) |
| `topopt/solver/results/<id>.json` | solver, **last** | `{protocol, id, threads, result \| error, density_parts, density_bytes}` |

The parts exist because the host lets a plugin read or write at most 1 MiB per call.

The cycle:
1. The solver writes `alive.json`; the plugin takes the solver as running while the heartbeat is at most 10 s old.
2. The plugin writes the job's parts and header. A job is complete when its parts add up to the header's `bytes`.
3. The solver takes complete jobs oldest first, reports `job` and `progress` in the heartbeat and solves.
4. It writes the mesh, the densities and, last, the result file; then it deletes the job's files.
5. The plugin reads the result on Refresh or Create body.

Every file the solver writes is written atomically (`.tmp`, then rename).

Another client can speak the same protocol (write a job and its header, wait for `results/<id>.json`; constants and types are in `topopt_study::solver`), but for external pipelines `solve` is simpler.

For the plugin: double-click `build/solver/start-topopt-solver.cmd`. It watches `build/results` unless given another folder; that folder must be the plugin's `results` folder. Minimize the window; closing it stops the solver (a running job stops, and the plugin then solves by itself).

## 10. Performance and grid size

Time grows roughly linearly with the number of elements and of iterations. Measured on the development machine (Intel i5-9400F, 6 threads):

| Grid | Per iteration | 150 iterations |
|---|---|---|
| 5 184 elements (36×12×12) | 0.05–0.07 s | ≈ 10 s |
| 16 000 elements (40×20×20) | ≈ 0.17 s | ≈ 25 s |
| 65 536 elements (64×32×32) | ≈ 0.5 s | ≈ 75 s |

Tips:
- Start at 20–40 thousand elements. Doubling the budget makes the element 1.26× smaller and the run twice as long.
- Give the thinnest part of the domain **at least 8 elements**, or multigrid converges slowly (you get a warning).
- Tune parameters on a coarse grid with `--max-iterations 40`; run the final study fine.
- `--threads` below the core count keeps the machine responsive; independent studies can run side by side with the cores split between them.

How the solver works inside (multigrid, element colouring, direct-solve thresholds): [topopt-core.md](topopt-core.md#performance).

With an NVIDIA GPU, `topopt-cuda-solver` takes the same arguments and runs the linear solves 3–4× faster (whole runs about 2× on the examples, since filtering, the update and meshing stay on the CPU): [topopt-cuda.md](topopt-cuda.md).

## 11. Troubleshooting

| Message | Cause and fix |
|---|---|
| `not a TopOpt spec: unknown field ...` | A typo in a field name; compare with `topopt-solver example` |
| `the spec has no supports: the part would float` | Add at least one support |
| `load N : its region holds no node of the part` | The region misses the material: check the coordinates (mm!), grow `radius` or `tolerance` |
| `the domain produced no elements: is it a closed solid, in millimetres?` | The STL is in metres (use `"scale": 1000`) or open |
| `the domain mesh is not closed` (warning) | The STL has holes; inside/outside may be wrong near them. Repair the mesh |
| `supports, loads and keep_in regions cover the whole part` | The pads took the whole part: reduce `pad` or the regions |
| `the body is only N elements thick somewhere` | Thin part or large element: lower `element_size` |
| `stopped at the N s time limit before converging` | Out of time: the result is there but could be better. Raise `time_limit` or lower the budget |
| Result volume well below `volfrac` | On a coarse grid, keeping the faces exact trims half-cells; finer grids reduce it |
| The plugin does not see the solver | Different folders: the folder the solver watches and the plugin's `results` folder must be the same |
