# TopOpt documentation

TopOpt finds the stiffest part that fits a design domain, its supports and loads, and a share of material, and delivers it as a **closed, manifold mesh** ready to become a solid. It runs as a SolidCraft plugin, as a standalone program, and as a Rust library.

## Where to start

| You want to… | Read |
|---|---|
| Optimize a part in SolidCraft | [plugin.md](plugin.md): the plugin's user guide |
| Pick the settings for a study | [settings-guide.md](settings-guide.md): what each study and mesh setting changes, measured on one bracket |
| Put optimization into a pipeline (scripts, CI, batches, your own software) | [topopt-solver.md](topopt-solver.md): CLI, spec format, outputs, progress, exit codes, Rust API, examples |
| Solve on an NVIDIA GPU | [topopt-cuda.md](topopt-cuda.md): TopOpt CUDA Solver and the topopt-cuda crate |
| Understand the algorithm or use the core directly | [topopt-core.md](topopt-core.md) |
| Mesh your own density fields or SDFs | [topopt-mesh.md](topopt-mesh.md) |
| See how geometry becomes a problem; plug in your own ray caster | [topopt-study.md](topopt-study.md) |
| Build, test, release | [development.md](development.md) |
| Know why it is built this way, and what is next | [design-notes.md](design-notes.md) |

## Components

| Component | Where | What it does |
|---|---|---|
| **topopt-core** | `crates/topopt-core` | SIMP on a Hex8 grid: matrix-free FE, multigrid-CG, density filter, Heaviside projection, OC; single- or multi-threaded (rayon) |
| **topopt-mesh** | `crates/topopt-mesh` | Density field → marching tetrahedra → smoothing pinned to the original surface → QEM → validation; STL and PLY in and out |
| **topopt-study** | `crates/topopt-study` | Studies (body, faces, material, settings), ray voxelization, problem assembly, the spec format, the plugin ↔ solver protocol, panel texts |
| **TopOpt plugin** | `plugins/topopt-src` → `plugins/topopt` | WebAssembly component for SolidCraft 0.49: dialogs, panel, ribbon, result insertion |
| **TopOpt Solver** | `tools/topopt-solver` | Native program and library: `solve` and `check` for pipelines, `watch` for the plugin |
| **topopt-cuda** | `crates/topopt-cuda` | The core's multigrid-CG on an NVIDIA GPU (mixed precision, CUDA driver API, kernels checked in as PTX) |
| **TopOpt CUDA Solver** | `tools/topopt-cuda-solver` | TopOpt Solver with the linear solves on the GPU; the same commands and protocol |
| References | `tools/reference` | Python: an independent optimizer (scipy) and a mesh checker (trimesh) |

## How they fit

```
                         ┌───────────────────── SolidCraft ─────────────────────┐
  user ── GUI ─────────► │ TopOpt plugin (wasm, one thread)                     │
                         │   prepare: view.pick → voxelization → Job            │
                         │   solve (when no solver runs) ─┐                     │
                         │   Create body ◄── PLY ─────────┤                     │
                         └──────────────┬─────────────────┼─────────────────────┘
                                        │ Job in parts    │ PLY, densities, summary
                                        ▼                 │
                              results folder (file protocol)
                                        │                 ▲
                                        ▼                 │
                         ┌──────────── TopOpt Solver (native, all cores) ───────┐
  pipeline ─ spec.json ─►│ solve / check / watch                                │──► mesh.ply|stl
                         │   prepare_spec: STL → MeshCaster → voxelization      │    density.json
                         │   solve: topopt-core (rayon) → topopt-mesh           │    result.json
                         └──────────────────────────────────────────────────────┘

  shared library:  topopt-core ◄── topopt-mesh ◄── topopt-study ◄── {plugin, solver}
                   topopt-core ◄── topopt-cuda ◄── TopOpt CUDA Solver (= TopOpt Solver + GPU)
```

A run has two halves: **prepare** (geometry → a `Job`) and **solve** (job → densities → mesh). Prepare needs geometry: SolidCraft's rays in the plugin, an STL in a spec. Solve is pure computation and identical everywhere, so the solver and the plugin give the same result on one study, and a job prepared in one place can be solved in another.

## What a result guarantees

Every output mesh is checked (`Mesh::validate`):
- closed (no boundary edges) and 2-manifold (no edge with three faces, no pinch vertices);
- outward-oriented, positive volume;
- no self-intersections, no degenerate triangles;
- no two vertices that merge when rounded to f32 (STL).

Checked independently with trimesh (`is_watertight`, consistent winding) and SolidCraft (`closed: true`, `mesh.convert` to a solid).

Core accuracy: agrees with an independent Python implementation to 1.45·10⁻⁸; cantilever deflection 0.982 of Timoshenko's; multigrid-CG equals plain CG to 10⁻¹².
