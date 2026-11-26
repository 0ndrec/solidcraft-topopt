# topopt-core — the optimization core

`crates/topopt-core` does density-based topology optimization (SIMP) on a structured 3D grid of Hex8 elements. It is pure Rust without I/O; the same core runs natively (the solver, the tests) and inside the SolidCraft plugin sandbox (`wasm32-unknown-unknown`).

The problem: **minimum compliance under a volume constraint**. Find the material layout in the design domain that makes the part stiffest under the given loads and supports while using at most a given fraction of material.

## Modules

| Module | What it does |
|---|---|
| `grid` | Numbering of nodes, elements and degrees of freedom; `coarsen()` for multigrid |
| `ke` | The Hex8 stiffness matrix (24×24, 2×2×2 quadrature) of a cube of side `h` |
| `problem` | The problem: grid, material, element kinds, supports, load cases |
| `fem` | Matrix-free `K·u` and diagonal, element energies |
| `mg` | Geometric multigrid: grid hierarchy, Jacobi smoothing, direct solve of the coarsest level |
| `pcg` | Preconditioned conjugate gradients |
| `filter` | Density filter (a cone of radius `rmin`), Heaviside projection and its derivative |
| `opt` | The optimization loop: `Optimizer`, `Settings`, `Iteration`, `Outcome` |
| `par` | Parallel loops (rayon, feature `parallel`) or plain loops |

## Units and grid

Millimetres, newtons, megapascals (N/mm²): displacements come out in mm, compliance in N·mm.

The grid has `nx × ny × nz` cubic elements of side `h`, starting at `origin`:
- node `(i, j, k)`, `0 ≤ i ≤ nx`, has index `i + (nx+1)·(j + (ny+1)·k)` and degrees of freedom `3·node + d`, `d ∈ {x, y, z}`;
- element `(i, j, k)` has index `i + nx·(j + ny·k)`;
- `node_pos(n)` and `elem_center(e)` give world coordinates.

## The problem: `Problem`

```rust
use topopt_core::{ElementKind, Grid, Material, Problem};

let grid = Grid::new(60, 20, 20, 2.0);            // 120×40×40 mm, h = 2 mm, origin (0, 0, 0)
let mut p = Problem::new(grid);
p.material = Material { e: 68_900.0, nu: 0.33 };

// Helpers by box (world coordinates, bounds inclusive):
p.fix_box([0.0, 0.0, 0.0], [0.0, 40.0, 40.0], [true; 3]);                 // clamp the x = 0 face
p.load_box([120.0, 0.0, 0.0], [120.0, 40.0, 0.0], [0.0, 0.0, -1000.0]);   // 1000 N on an edge
p.set_kind_box([50.0, 10.0, 10.0], [70.0, 30.0, 30.0], ElementKind::Void); // a cut-out
```

| Field | Meaning |
|---|---|
| `grid` | The grid |
| `material` | `e` (MPa), `nu`; default aluminium-like 70 000 / 0.33 |
| `kinds` | One kind per element: `Design` (the optimizer chooses its density), `Solid` (always 1), `Void` (always 0). Empty means all `Design` |
| `supports` | `NodeFix { node, dirs: [bool; 3] }`: a node held along the chosen axes |
| `load_cases` | `LoadCase { weight, forces: Vec<NodeForce { node, f: [f64; 3] }> }`; the objective is the weighted sum of compliances |

`validate()` checks the length of `kinds` and the node indices, that there are supports and loads, non-negative weights, E > 0 and at least one `Design` element. `Optimizer::new` calls it.

`Problem` is serde-serializable; that is how the plugin hands a problem to the solver.

## Settings

| Field | Default | Meaning |
|---|---|---|
| `volfrac` | 0.3 | Material fraction of the `Design` elements |
| `penal` | 3.0 | SIMP penalty `p`: `E(x) = E_min + x^p·(E − E_min)` |
| `emin` | 1e-9 | Void stiffness relative to solid |
| `rmin` | 1.5 | Filter radius, **in elements** |
| `projection` | true | Heaviside projection with β-continuation |
| `eta` | 0.5 | Projection threshold |
| `beta_start`, `beta_max` | 1, 16 | Initial and final sharpness |
| `beta_every` | 40 | Double β at least this often (earlier when the design settles) |
| `move_limit` | 0.2 | Largest density change per iteration (OC) |
| `max_iter` | 300 | Iteration limit |
| `tol` | 0.01 | Converged when the largest density change at the final β is below this |
| `compliance_tol` | 1e-3 | …or when compliance changed by less than this (relative) in each of the last five iterations at the final β |
| `pcg_tol`, `pcg_max_iter` | 1e-6, 200 | CG tolerance and iteration limit |
| `mg` | below | Multigrid settings |

## The optimizer

```rust
use topopt_core::{Optimizer, Settings};

let mut opt = Optimizer::new(p, Settings { volfrac: 0.25, ..Default::default() })?;
let outcome = opt.run(|it| {
    println!("{} c={:.4} v={:.3} β={} cg={}", it.iter, it.compliance, it.volume, it.beta, it.cg_iterations);
    true                      // false stops the run (Outcome::Cancelled)
})?;
let densities = &opt.x_phys; // physical densities 0…1 per element
```

| Method or field | Purpose |
|---|---|
| `new(problem, settings)` | Checks, element matrix, filter, multigrid hierarchy |
| `run(callback) -> Outcome` | Iterates until converged (`Converged`), the limit (`MaxIterations`) or `false` from the callback (`Cancelled`) |
| `step() -> Iteration` | One iteration: FE solve, sensitivities, filter, OC |
| `evaluate() -> Evaluation` | Compliance and gradients `dc`, `dv` at the current design, without a step |
| `set_design(&x)` | Set the densities of the `Design` elements (warm start, gradient checks) |
| `x`, `x_filtered`, `x_phys` | Variables, filtered, projected |
| `u` | Displacements per load case (CG warm starts between iterations) |
| `history: Vec<Iteration>` | `iter`, `compliance`, `volume`, `change`, `beta`, `grayness`, `cg_iterations`, `cg_worst_residual` |
| `mg_levels()` | The multigrid level grids (diagnostics) |

### One iteration

1. Variables `x` → density filter → `x̃` → Heaviside projection → `x̄`; `Solid` and `Void` elements are forced to 1 and 0.
2. Element moduli `E_e = E_min + x̄_e^p·(1 − E_min)` (a factor on the element matrix).
3. For each load case, solve `K(E)·u = f` with MG-PCG, warm-started from the previous `u`.
4. Compliance `c = Σ w·fᵀu`; sensitivity `∂c/∂x̄_e = −p·x̄^(p−1)·u_eᵀ k u_e`.
5. Back through the projection and filter (chain rule) to `∂c/∂x` and `∂V/∂x`.
6. Optimality-criteria update with bisection on the volume multiplier, within `move_limit`.
7. β-continuation: double β when the design settles or after `beta_every` iterations.

## Linear solver

**Matrix-free `K·u`.** No global matrix is assembled: for each element, gather its 24 displacements, multiply by `E_e·k`, scatter back. With the `parallel` feature the elements are coloured by the parity of `(i, j, k)` into 8 colours whose elements share no node, so threads add into `y` without locks.

**Geometric multigrid as the CG preconditioner** (`mg::Hierarchy`):
- levels halve the grid along each axis while a level has more than `max_coarse_dofs` degrees of freedom;
- coarse-level matrices come from **element-wise Galerkin coarsening** `Pᵀ K P` (a coarse element is eight fine ones), stored in f32: half the memory and memory traffic, the same CG iteration count;
- damped Jacobi smoothing with step `ω·2/λ_max(D⁻¹A)`, `λ_max` estimated per level by power iteration (warm-started between updates). Without it, multigrid diverged on high-contrast density fields;
- the coarsest level is solved **directly** (faer's sparse Cholesky, symbolic analysis done once).

| `MgSettings` | Default | Meaning |
|---|---|---|
| `max_coarse_dofs` | 4000 | Coarsen while a level has more degrees of freedom than this |
| `direct_dofs` | 4000 | Solve a grid up to this size directly, without multigrid |
| `smooth_steps` | 2 | Jacobi sweeps before and after the coarse correction |
| `omega` | 0.6 | Fraction of the Jacobi stability limit |
| `power_steps` | 3 | Power-iteration steps per update |

## Performance

Measured with `examples/profile`, native, 6 threads (Intel i5-9400F). Time per optimization iteration:

| Grid | Direct solve | Multigrid-CG |
|---|---|---|
| 36×12×12, 18 759 dof | 370 ms | **52 ms** |
| 20×10×10, 7 623 dof | 82 ms | **31 ms** (with a 10×5×5 coarse level) |
| 64×32×4 thin plate, 160 iterations | 84 s | **47 s**, same compliance |
| 24×8×8, 6 075 dof | 43 ms | 44 ms |

So the direct-solve thresholds are the same natively and in wasm (4000 dof). In wasm a factorization is another 8× slower: a 4000-element study took 65 s with direct solves and 7.7 s without.

Other measurements:
- 65 536 elements: ≈ 2.2 s per iteration in the plugin sandbox (one thread), ≈ 0.48 s natively on 6 threads;
- matrix-free product on 65k elements: 7.7 ms natively on one thread (after rewriting the kernel: masked copy of `x`, incremental indices, no bounds checks), 14.3 ms in the sandbox. `+simd128` gave no measurable gain: the kernel is bound by gather and scatter;
- on thin domains (4 elements across) CG grows to ~100 iterations at high β: the coarse grid loses members one or two elements thick. Semi-coarsening (not coarsening the thin axis) would fix it; until then, give the thinnest direction ≥ 8 elements.

## Verification

| Check | Result |
|---|---|
| Independent Python implementation (`tools/reference/top3d_ref.py`: assembled K, scipy `spsolve`, its own Hex8 and numbering), 40×16×8, 30 iterations without projection | compliance agrees to **1.45·10⁻⁸** at every iteration |
| Cantilever deflection against Timoshenko (`tests/solver.rs`) | ratio 0.982 (fully integrated Hex8 is slightly stiff) |
| MG-PCG against Jacobi-CG, 32×16×16, uniform / random contrast 1e-3…1 | solutions agree to 10⁻¹²; 15 / 33 iterations against 245 / 659 |
| Gradient against finite differences (`tests/optimize.rs`), with and without projection, with passive elements | relative error < 1e-4 |
| Small cantilever converges to a crisp (0/1) design | low `grayness` |
| Coloured parallel `K·u` against the sequential one; assembled matrix symmetric | equal |
| Hex8: symmetric, zero energy for rigid modes, exact uniaxial-stretch energy | holds |
| Element-wise Galerkin against assembled `Pᵀ K P` | exact |

```sh
cargo +stable-x86_64-pc-windows-gnu test --release -p topopt-core
```

## Examples

```sh
# The top3d benchmark: x = 0 clamped, line load along y on the free bottom edge; one JSON line per iteration
cargo +stable-x86_64-pc-windows-gnu run --release -p topopt-core --features parallel --example cantilever -- 60 20 4 100

# Where an iteration's time goes: NX NY NZ [SMOOTH_STEPS [DIRECT_DOFS [STEPS [MAX_COARSE_DOFS]]]]
cargo +stable-x86_64-pc-windows-gnu run --release -p topopt-core --features parallel --example profile -- 36 12 12 2 4000 20
```

## Features and dependencies

- `parallel`: rayon for the heavy loops (`K·u`, filter, energies, Galerkin, transfers, vector operations). Without it the code is single-threaded, as built for the plugin. `par::threads()` gives the thread count; `RAYON_NUM_THREADS` sets it.
- `faer` 0.24 (only `std`, `linalg`, `sparse`, `sparse-linalg`): sparse Cholesky.
- `serde`: serialization of `Problem`, `Settings`, `Iteration`.

## Out of scope

- Surfaces: that is [topopt-mesh](topopt-mesh.md).
- CAD, STL and regions: [topopt-study](topopt-study.md) turns geometry into a `Problem`.
- Only linear elasticity, isotropic material, compliance with a volume constraint. No stress, frequency or overhang constraints yet; they would go here, as a term in `evaluate()` with its sensitivity.
