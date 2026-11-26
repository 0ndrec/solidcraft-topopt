# topopt-study — studies: from geometry to a problem and back

`crates/topopt-study` connects part geometry to the core ([topopt-core](topopt-core.md)) and the mesher ([topopt-mesh](topopt-mesh.md)). It holds everything that does not depend on SolidCraft but knows about a study:

- the **study** (`Study`): design body, keep-outs, material, settings, supports and loads on faces, the last result. The plugin stores it in the design;
- **voxelization** of a body by ray casting through an abstract `Caster` (`view.pick` in the plugin, parry3d in tests and specs);
- **problem assembly** (`Problem`): support and load nodes, solid pads, voids;
- **solving**: optimization and mesh, with progress, cancel and clock through `Hooks`;
- the **spec**: a study without CAD (regions instead of faces), the input of `topopt-solver solve`;
- the plugin ↔ TopOpt Solver **file protocol**;
- **panel texts**: the "next step" hint, status rows, time estimate.

The same code is compiled to wasm in the plugin and natively in the solver, so both give the same result on one study (a test bracket: 18 898.9 mm³ in both).

## Two halves of a run

```
             ┌──────────── prepare (needs CAD: rays) ─────────────┐   ┌────────── solve (pure computation) ──────────┐
Study + body │ plan_grid → voxelize → build → finish_job → Job    │ → │ Optimizer → Density → make_mesh → Mesh       │
             └────────────────────────────────────────────────────┘   │ RunResult (summary)                          │
Spec         │ shape_mesh → voxelize(MeshCaster) → assemble → Job │   └──────────────────────────────────────────────┘
             └────────────────── prepare_spec ────────────────────┘
```

`Job` is the boundary. It serializes to JSON and holds everything `solve` needs: the study, the body's box and volume, the `Problem` with element kinds, the core settings, warnings, the original surface for vertex pinning and the mesh path. The plugin runs `prepare` itself, as only it can cast rays at the model; `solve` runs either in the plugin or in TopOpt Solver.

## API

### Entry points

| Function | Purpose |
|---|---|
| `prepare(study, body, design_caster, keep_out_casters, surface, hooks) -> Job` | Voxelize the body and keep-outs, assemble the problem from the study's faces |
| `prepare_spec(spec, read_file, hooks) -> Job` | The same from a spec; `read_file(path)` returns a file's bytes (STL), so the caller decides where files come from |
| `solve(job, hooks) -> Outcome { mesh, density, result }` | Optimize and mesh. Multigrid settings are the defaults **of the platform it runs on**, not taken from the job |
| `run(...)` | `prepare` + `solve` in one call |
| `make_mesh(study, density, surface) -> (Mesh, MeshReport)` | Only the mesh, from saved densities (Re-mesh) |

### `Hooks`

```rust
pub trait Hooks {
    fn progress(&mut self, fraction: f32, message: &str);  // 0…1 and a text for people
    fn cancelled(&mut self) -> bool;                        // polled every iteration
    fn elapsed(&mut self) -> f64;                           // seconds since the start (0 without a clock)
}
```

Progress: 0–0.04 voxelization, 0.05–0.9 iterations, 0.9 mesh, 1.0 done. Cancelling returns `Err("cancelled")`. `settings.time_limit` is compared with `elapsed()`: optimization stops, and the mesh is still made.

### `Caster` and voxelization

```rust
pub trait Caster {
    /// For each ray (origin, unit direction): the nearest hit on the body ahead of the origin.
    fn cast(&mut self, rays: &[([f64; 3], [f64; 3])]) -> Result<Vec<Option<Hit>>, String>;
}
pub struct Hit { pub point: [f64; 3], pub normal: [f64; 3], pub face: u32 }
```

- `plan_grid(lo, hi, element_size, max_elements)`: the grid over the body's box. The element edge is given, or derived from the budget; dimensions are multiples of 4 (of 2 when small) so multigrid can coarsen, and the grid is centred on the box. Rounding can leave the element count well below the budget.
- `voxelize(caster, grid, progress) -> Voxels { inside, hits, unbalanced_rays }`. Along each axis, one ray runs through every row of cell centres, hit by hit (the host reports only the first hit; the next ray starts just past it). Entering hits open an inside interval, exiting hits close it. A cell is inside when at least two of the three axes say so, which forgives a ray grazing an edge. Every hit is kept with its face: supports, loads and pads are placed from them.
- `MeshCaster::new(&Mesh)`: a `Caster` over a triangle mesh (parry3d). Inside and outside come from the parity of crossings along the whole line, so triangle winding does not matter.

### The study (`study`)

| Type | Purpose |
|---|---|
| `Study` | `body`, `keep_out`, `material`, `settings`, `mesh`, `conditions`, `result`, `pending` (a job at the solver), `last_error` |
| `Condition::Support { id, faces, dirs }` | Faces held along the chosen axes |
| `Condition::Load { id, faces, mode, force, magnitude, case }` | `mode = "vector"`: total force `force`; `"normal"`: `magnitude` N into the faces (along their mean normal) |
| `FaceRef { body, index, point, label }` | A face: body, face number (`model.faces`), the picked point, a label |
| `Settings`, `MeshSettings` | Fields and defaults in [topopt-solver.md §4.7–4.8](topopt-solver.md#47-settings) |
| `Material`, `MATERIALS` | 9 presets (`aluminium-6061` … `nylon`): E, ν, density |
| `RunResult` | The run summary ([topopt-solver.md §5.3](topopt-solver.md#53-summary-resultjson)) |
| `Study::add / remove / set_option / problems` | Edits from dialogs and panel events; what is missing before a run |

`set_option(key, value)` takes the panel's keys (`settings/volume_percent`, `settings/use_solver`, `material/label`, `mesh/triangles`…). Values may be bare JSON scalars or objects like `{"number": 0.3}`, as the host sends them, and are range-checked.

### How faces become boundary conditions (`build`)

1. Ray hits are grouped by face; each condition takes the points of its faces (plus the picked point).
2. Nodes within one element of those points that touch material become support or load nodes; a force is spread evenly over them.
3. Cells within `pad` (default 1.5 elements) become **solid pads** (`ElementKind::Solid`): real mounting surfaces, and the anchors of mesh clean-up.
4. Cells outside the body or inside a keep-out are `Void`.
5. `finish_job`: the volume fraction is converted from a fraction of the body to a fraction of the design elements (accounting for the pads) and clamped to 0.05…0.95; the filter radius is `rmin = max(1.5, min_member / 2h)`.

Specs (`spec::assemble`) work the same way, except that nodes are chosen by **regions**, preferably on the surface of the material ([topopt-solver.md §4.3](topopt-solver.md#43-regions-supports-loads-keep_in)).

### Specs (`spec`)

`Spec`, `Shape` (`Box`, `Stl`, `Mesh`), `Region` (`Box`, `Sphere`, `Plane`), `SupportSpec`, `LoadSpec`, `MaterialSpec`, `Spec::example()`, `box_mesh(lo, hi)`. The format: [topopt-solver.md §4](topopt-solver.md#4-the-spec-format).

### Solver protocol (`solver`)

Path constants (`ALIVE`, `JOBS`, `RESULTS`, `job_header`, `job_part`, `job_cancel`, `result_file`, `density_part`), `PART_BYTES = 1 000 000`, `ALIVE_SECONDS = 10`, `PROTOCOL = 1`. Types `Alive` (with `online(now)`), `JobHeader`, `ResultFile`, `Pending`, `SolverView`; `split(bytes)` cuts data into parts. The protocol: [topopt-solver.md §9](topopt-solver.md#9-watch-mode-the-plugin-protocol).

### Panel state (`status`)

- `next_step_with(study, solver)`: the one-line "what next" above the panel's tabs.
- `rows_with(study, body, solver)`: the Study table on the Run tab.
- `estimate_seconds(study, elements, threads)`: a rough time estimate; `threads = 0` means the plugin sandbox.
- `fmt_duration`, `fmt_n` (numbers with thin-space grouping).

## Example: your own application with your own ray caster

```rust
use topopt_study::{BodyInfo, Caster, Condition, FaceRef, Hit, Hooks, Study, prepare, solve};

struct MyCaster { /* your BVH */ }
impl Caster for MyCaster {
    fn cast(&mut self, rays: &[([f64; 3], [f64; 3])]) -> Result<Vec<Option<Hit>>, String> {
        Ok(rays.iter().map(|(o, d)| self.first_hit(*o, *d)).collect())   // point, outward normal, face id
    }
}

let mut study = Study { body: Some("Bracket".into()), ..Default::default() };
let face = |index| FaceRef { body: "Bracket".into(), index, point: None, label: None };
study.add(Condition::Support { id: 0, faces: vec![face(5)], dirs: [true; 3] }, None);
study.add(Condition::Load { id: 0, faces: vec![face(3)], mode: "vector".into(),
                            force: [0.0, 0.0, -1000.0], magnitude: 0.0, case: 1 }, None);
let body = BodyInfo { lo: [0.0; 3], hi: [80.0, 40.0, 20.0], volume: 64_000.0 };
let job = prepare(&study, body, &mut MyCaster { .. }, &mut [], None, &mut hooks)?;
let out = solve(&job, &mut hooks)?;
```

Without face numbers (just a mesh or an STL), a spec is simpler: `prepare_spec` builds its own `MeshCaster`.

## Tests

| File | Checks |
|---|---|
| `tests/run.rs` | Box voxelization (volume within 4 %, all 6 faces sampled, no unbalanced rays); an end-to-end "support −X, load +X" study; a "normal" load equals the vector |
| `tests/spec.rs` | The example spec solves; box, inline mesh and STL in metres (`scale` 1000) voxelize alike; keep-outs; clear region errors; JSON round trip |
| unit | `set_option` and study JSON, "next step" hints through the steps, protocol parts and heartbeat |

```sh
cargo +stable-x86_64-pc-windows-gnu test --release -p topopt-study
```
