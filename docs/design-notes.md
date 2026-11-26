# Design notes

Why the system is built the way it is: what the SolidCraft host allows, the choices made, the measurements behind them, and what is left to do. Measurements are from the development machine (Intel i5-9400F, 6 cores, 16 GB, Windows 11) and SolidCraft 0.49.1.

## 1. The SolidCraft host

Plugins are WebAssembly components (`wasm32-unknown-unknown`, no WASI) built against the WIT world `craft:plugin@0.2.0`. Everything below was verified by running it, headless (`solidcraft-cli run`) or in the app.

| Question | Answer |
|---|---|
| `view.pick {rays}` | **One hit per body**, the nearest; `all` lists bodies, not entry and exit points |
| A ray started just past an entry hit | Finds the **exit** (the inside of the surface). Iterative ray voxelization works: one call of up to 1024 rays per crossing |
| `inspect.classify_points` | Exact, usable from a read entry |
| `model.mesh` | Shared vertices; at most 8192 vertices / 16384 triangles |
| `model.faces` | Name, normal, point, area, centroid per face |
| Writing a file and `file.insert_mesh` **in one run** | Impossible: a run's writes are staged and published only when it ends, so the app does not see the file yet |
| The same in **two runs** (write, then insert) | Works; `grant:results/...` resolves to an absolute path |
| Mesh body from a closed STL | `closed: true` only if shared corners are **bit-identical** (e.g. `-0.0` vs `0.0` breaks it): the mesher must produce an indexed mesh |
| `mesh.convert` | 1000 triangles → solid in 156 ms; 30 000 in ≈ 5 s (one planar face per triangle) |
| `storage.put-blob` / `get-blob` | 8 MiB round trip in 168 ms |
| `clock.now()` | Seconds, not milliseconds |
| `chart` panel widget | A placeholder in this version |
| Panels | The host does not scroll plugin panels; only tables and lists scroll |
| Ribbon | Every entry appears in Model › Plugins. The label is the `ribbon` slot's name, else the `menu` slot's, else the entry name; the icon is `x-solidcraft-icon` from the input schema |

Limits that shaped the design:

| Limit | Normal entry | `"job": true` entry |
|---|---|---|
| Time | 5 s default, ≤ 30 s | **10 minutes** |
| Host calls | **256** | 1 048 576 |
| Fuel | 10 M default; app ceiling 100 G; up to **10 T** when approved at install | same |
| Memory | 64 MiB default; up to 4096 MiB when approved | same |
| File read / write | 1 MiB per call | same |
| Manifest | ≤ 16 KiB | |

Fuel, not time, is the real limit: 100 G of fuel is about 28 s of computation, and inside the app fuel runs out faster than headless. The sandbox has one thread, no GPU, no processes.

The matrix-free Hex8 product runs at about 2 GFLOP/s in the sandbox (Cranelift): ~150 ms per product on 64³ elements. The first run after an install compiles with the baseline compiler (Winch) and is about twice as slow; the compiled code is cached.

## 2. Choices

**Structured grid, Hex8, SIMP + OC, density filter, Heaviside projection with β-continuation.** The classic, well-understood formulation (top3d, Liu & Tovar 2014; top99neo, Ferrari & Sigmund 2020; robust projection, Wang, Lazarov & Sigmund 2011). A structured grid makes the stiffness operator matrix-free and multigrid geometric.

**Matrix-free geometric multigrid-preconditioned CG** (Amir, Aage & Lazarov 2014). No sparse assembly, small memory, and the iteration count barely grows with the grid. Element-wise Galerkin coarsening gives the exact coarse operator without a global matrix. λmax-based Jacobi damping was needed: a fixed ω diverged on high-contrast density fields.

**Marching tetrahedra instead of marching cubes.** On the Freudenthal split every tetrahedron yields a flat triangle or a convex quad, so the surface is watertight and manifold by construction, with no ambiguous cases. Component clean-up uses the same 14-neighbour connectivity, so after it the mesh is exactly one shell without voids.

**Our own decimation instead of the host's `mesh.reduce`**: it respects pinned vertices and checks manifoldness.

**Pinning to the original surface.** Mounting faces and holes must stay exact, so vertices near the original surface are projected onto it during smoothing. Pinning blindly folded triangles at edges (both sides of an edge projected onto different faces); now a vertex is pinned only when the projection direction agrees with the normal of a freely smoothed copy, and vertices of offending triangles are unpinned selectively.

**Blur off by default.** It moves the iso-surface inwards on thin convex members (a torus loses 7 % of its volume against 1.8 % without), i.e. thinner, weaker struts. Taubin smoothing removes the voxel steps instead.

**Nudge 0.02 off the iso value.** With 1e-3, two vertices near one sample could land closer than the f32 resolution of STL (≈ 8 µm at 100 mm) and SolidCraft welded them.

**PLY for the plugin's output.** The 1 MiB file-write limit: indexed binary PLY is about half the size of STL; meshes are capped at 50 000 triangles.

**Two runs for a result:** Optimize writes the mesh, Create body inserts it (the staged-write rule).

**A native solver beside the app.** The sandbox has one thread, about 17 % of a 6-core CPU. TopOpt Solver is the same code built natively with rayon; the plugin hands jobs over through its `results` folder ([topopt-solver.md §9](topopt-solver.md#9-watch-mode-the-plugin-protocol)). The plugin still voxelizes (only it can cast rays at the model) and solves by itself when the solver is not running.

**Direct-solve thresholds of 4000 dof everywhere.** Multigrid-CG beats faer's sparse Cholesky well below 10 000 dof natively too, also on thin domains (measurements in [topopt-core.md](topopt-core.md#performance)); in wasm the factorization is another 8× slower.

**Toolchain.** The development machine has no MSVC linker, so everything builds with the self-contained `stable-x86_64-pc-windows-gnu` toolchain; `wasm-tools` does not build under it (`windows-sys` needs `dlltool`), so `xtask` componentizes with `wit-component` directly and refuses imports other than `craft:plugin/*`.

## 3. Measurements

| What | Result |
|---|---|
| Core against an independent Python implementation (assembled K, scipy), 40×16×8 | compliance agrees to 1.45·10⁻⁸ every iteration |
| Cantilever against Timoshenko | 0.982 |
| MG-PCG against Jacobi-CG, 32×16×16 | equal to 10⁻¹²; 15 / 33 iterations against 245 / 659 |
| 65 536 elements, per iteration | sandbox ≈ 2.2 s; native 1 thread 1.4 s; native 6 threads ≈ 0.48 s |
| Bracket, 21 448 elements, 60 iterations, whole Optimize | in the plugin 66 s; with TopOpt Solver 1.6 s hand-over + 20 s; identical volume 18 898.9 mm³ |
| Cantilever 48×24×24 → mesh | 0.36 s from densities to a 30 000-triangle mesh; trimesh watertight; SolidCraft `closed: true`, converted to a solid in 5.1 s with the volume unchanged to 1e-10 |
| Example spec, 5 184 elements, 150 iterations, native | 56 s with the old 40 000-dof direct threshold, 10 s now |

## 4. GPU

Inside the plugin: impossible (no GPU access; plugin WGSL shaders cannot use storage buffers or compute). Outside, through the solver protocol, it is possible. On the development machine's Quadro P620 (FP64 at 1/32 rate) a CUDA FP64 port would match the 6-thread CPU solver, so the multi-threaded CPU solver came first. A GPU pays off only in FP32 (mixed precision: FP32 V-cycle and products, FP64 outer CG) and for grids beyond ~250k elements. Vendor-neutral `wgpu` or `cudarc` would be the routes; NVIDIA's `cuda-bindings` needs the CUDA Toolkit and MSVC at build time, and its Rust kernels (`cuda-oxide`) require Ampere or newer.

## 5. Licences of the sources consulted

- **ninja7v/TopoptComec** (MIT): reference for the algorithms (Hex8, OC, filter); not run, the independent Python reference checks the core instead.
- **dapsvi/phoenix** (GPLv3): ideas only (problem model, post-processing order); no code taken.
- **SolidCraft** `plugin.wit` (MIT OR Apache-2.0): copied unchanged into `wit/`.

## 6. Roadmap

Quality:
- [ ] Account for the fill fraction of boundary voxels (on coarse grids the result is 2–8 % below the target volume).
- [ ] Sharp edges of the original body when projecting vertices.
- [ ] `plan_grid`: rounding dimensions to multiples of 4 can waste half the budget on small grids (a 4000 budget gives 1536 elements on 120×40×40); choose the dimensions closest to the budget from below.
- [ ] Recalibrate `estimate_seconds` for the new multigrid thresholds (native estimates are high).
- [ ] Semi-coarsening for thin domains.

Plugin:
- [ ] Solver progress in the panel without Refresh, and a Cancel for solver jobs.
- [ ] Convergence chart once the host's `chart` widget works.

Extensions, by priority:
1. MMA (Svanberg) for several constraints.
2. Overhang filter for 3D printing (Langelaar).
3. Stress constraints (p-norm, MMA).
4. Optional GPU backend behind the solver protocol.
5. Eigenfrequencies and self-weight.
6. Symmetry planes.
