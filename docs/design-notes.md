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

**Our own decimation instead of the host's `mesh.reduce`**: it respects pinned vertices and checks manifoldness. It must be deterministic: it once iterated a hash set, so the same input sometimes left a self-intersection and fell back to a 1.5× larger budget; neighbourhoods are now sorted, and crossings lock their neighbourhood and decimate again before the budget grows.

**Sharp edges by inserting vertices, not by snapping them.** Marching tetrahedra chamfers or rounds the creases of the part (on two of a box's four parallel edges the split's diagonals cut across them, by up to 0.5 h). Snapping the nearest contour vertices onto the crease collapsed triangles (two chamfer vertices land on the same point) and was undone by the self-intersection fallback. Splitting every mesh edge that crosses the crease's bisector plane at the crease puts vertices exactly on it while the chamfer's vertices stay on their faces; the box comes out exact, and a real cantilever's clamp edges went from 0.66 mm to 0.00 mm off.

**Bilateral filtering on normals, not on vertices.** Fleishman et al.'s vertex filter, on marching-tetrahedra triangles, pulled the walls of thin members together and crossed the surface hundreds of times. Filtering triangle normals (Zheng et al.) and moving each vertex towards its own triangles' planes (Sun et al.) is stable and gave the most accurate surfaces in the comparisons ([topopt-mesh.md](topopt-mesh.md#measurements)).

**Round members: mesh the filtered field, finer.** The optimizer's projected densities are nearly 0 or 1 per element, so any surface through them is voxels smoothed afterwards; struts three or four elements thick come out beaded. The filtered field before the projection carries the same surface (the projection threshold is the iso value) but varies smoothly; sampled ×2 by cubic interpolation and smoothed bilaterally, strut cross-sections on a bridge fit an ellipse to 0.984 instead of 0.956. What did not work: smoothing a distance field (it closed 26 of 31 holes), mean-curvature flow (it crossed the surface on thin parts and fell back to no smoothing), stronger smoothing of the projected field alone (0.971).

**Rounded joints by a rounded union of distance fields.** As implicit CAD kernels do, the optimized material and the kept parts (pads, keep-in) get exact signed distance fields, and the field is raised to their rounded union only in the fillet zone. The circular round union is used, not the polynomial smooth minimum: the latter's fillet is about 0.6 × its parameter and it adds a quarter of it along the whole contact. Without smoothing the raised field, the half-cell steps of voxel distances made the bridge's smoothing back off entirely. The default radius (1.5 elements) changes little on our examples, because the filtered field already rounds joints by about that much.

**Folds are checked.** Upsampling plus crease insertion once produced triangles folded flat onto each other in a planar face: no crossing, but SolidCraft's `mesh.convert` refused the mesh ("a planar group's boundary does not enclose it"). Folds are now part of `Mesh::validate` and of the smoothing safety loop.

**Thin plates along a face keep their thickness.** A face pin used to accept a vertex whose normal was parallel to the face either way, so the inner side of a one-element plate was pinned onto the face too, flattening the plate (a bridge deck: 31 % of the part's volume). A vertex is now pinned only if its outward normal agrees with the face's.

**"Stalactites" are mostly load-carrying.** Thin hooks and strands hanging from members looked like dead ends. Selecting elements by low strain energy found almost none (the lowest 1 % still had 20 % of the mean energy); removing thin pieces one by one with re-analysis showed they cost more stiffness than their share of the mass (cantilever: 0.9 % of the material, 1.4 % of the stiffness). Under the bridge deck they are ribs stiffening the thin loaded deck: a thicker load pad removes them, but changes the whole design. Hence trimming exists with an explicit stiffness budget and is off by default.

**Volume preservation restores the contour's volume, with a bounded offset.** Smoothing shrinks; offsetting the free vertices along their normals (Newton on the exact volume derivative) undoes that for any smoother. When almost every vertex is pinned, the few free ones would have to carry the whole correction and grew 0.15 h bumps, so the offset is limited by the free share of the surface.

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
- [x] Sharp edges and corners of the original body (crease insertion, ranked pins in decimation).
- [x] Volume-preserving smoothing; HC and bilateral smoothers.
- [ ] Sharpen creases the optimizer makes inside the domain (feature detection on the result itself).
- [x] Round members: filtered field, ×2 cubic upsampling, bilateral smoothing.
- [x] Optional trimming of thin dead weight within a stiffness budget.
- [x] Rounded joints between the optimized material and the kept parts.
- [ ] Skeleton-based reconstruction: struts as exact cylinders (capsules) blended at the nodes.
- [ ] Conformal lattices and ribs on surfaces.
- [ ] Faster ×2 meshing (the smoothing and decimation of 4× more triangles dominate).
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
