# Development: build, test, release

## Tools

| Tool | Why |
|---|---|
| Rust stable, toolchain **`stable-x86_64-pc-windows-gnu`** | The development machine has no MSVC `link.exe`; the GNU toolchain is self-contained. Every cargo command uses `+stable-x86_64-pc-windows-gnu` |
| Target `wasm32-unknown-unknown` | The plugin. `.cargo/config.toml` enables `+simd128` (Wasmtime with Cranelift supports it) |
| `solidcraft-cli` (SolidCraft 0.49.1) on `PATH` | `plugin check`, headless runs, installs |
| Python 3 | Manifest and icon generators, test scripts; run with `-I` |
| `tools/reference/.venv` (numpy, scipy, trimesh) | Independent checks of the core and of meshes |
| Git Bash | The `scripts/*.sh` |

```sh
rustup toolchain install stable-x86_64-pc-windows-gnu --profile minimal --target wasm32-unknown-unknown
python -m venv tools/reference/.venv && tools/reference/.venv/Scripts/pip install numpy scipy trimesh
```

## Layout

```
Cargo.toml                 workspace: shared serde_json and wit-bindgen; release profile opt 3, LTO, 1 codegen unit
wit/plugin.wit             SolidCraft's WIT (craft:plugin@0.2.0)
xtask/                     componentize (instead of `wasm-tools component new`), list imports
crates/topopt-core/        the core: FE, multigrid-CG, SIMP                 → docs/topopt-core.md
crates/topopt-mesh/        field → closed mesh; STL/PLY                       → docs/topopt-mesh.md
crates/topopt-study/       studies, voxelization, specs, protocol             → docs/topopt-study.md
crates/topopt-cuda/        the multigrid-CG on an NVIDIA GPU                  → docs/topopt-cuda.md
  kernels/                 its GPU kernels: Rust for nvptx64 (nightly), outside the workspace
  src/kernels.ptx          the built kernels, checked in (scripts/build-ptx.sh)
plugins/topopt-src/        the plugin: entries, host glue, panel
plugins/topopt/            the plugin package: plugin.json (generated), icons/, topopt.wasm (built)
tools/topopt-solver/       TopOpt Solver: library + CLI + launcher            → docs/topopt-solver.md
tools/topopt-cuda-solver/  TopOpt CUDA Solver: the CLI with a GPU factory + launcher
tools/reference/           Python references: top3d_ref.py, check_stl.py
examples/                  example specs for topopt-solver
scripts/                   builds, generators, headless tests, render_mesh.py (PNG renders of results)
build/                     build output (ignored): topopt.zip, solver/, results/
```

Crate dependencies:

```
topopt-core ← topopt-mesh ← topopt-study ← topopt-plugin (wasm)
     ↑                            ↑
     ├──── (feature parallel) ── topopt-solver (native) ◄── topopt-cuda-solver
     └──── topopt-cuda ─────────────────────────────────────────┘
```

## Build

### The plugin

```sh
python -I scripts/gen_topopt_icons.py
python -I scripts/gen_topopt_manifest.py
scripts/build-plugin.sh topopt
```

`build-plugin.sh`:
1. `cargo build --release --target wasm32-unknown-unknown -p topopt-plugin`;
2. `xtask component` turns the core module into a component and refuses any import other than `craft:plugin/*`;
3. `solidcraft-cli plugin check plugins/topopt` checks the manifest, icons and imports;
4. packs `plugins/topopt` into `build/topopt.zip`.

Do **not** edit `plugin.json` by hand: `gen_topopt_manifest.py` writes it, and the plugin's `ui-provider.describe()` reads the same file, so the two cannot drift. The generator writes compact JSON and asserts the host's 16 KiB limit. The plugin's version is set there too.

### TopOpt Solver

```sh
scripts/build-solver.sh
```

Builds both release binaries, `topopt-solver` and `topopt-cuda-solver`, and copies them with their launchers into `build/solver/`. A running solver keeps its exe locked on Windows, so the script renames it out of the way first.

### The GPU kernels

```sh
rustup toolchain install nightly --profile minimal --component llvm-tools,llvm-bitcode-linker,rust-src --target nvptx64-nvidia-cuda
scripts/build-ptx.sh
```

Only after changing `crates/topopt-cuda/kernels/lib.rs`: the PTX is checked in, so normal builds need no nightly and no CUDA toolkit. The test `ptx_matches_the_kernel_source` catches a stale PTX. Details: [topopt-cuda.md](topopt-cuda.md#6-rebuilding-the-kernels).

## Tests

```sh
# Everything native: core, mesher, studies, solver protocol and CLI
cargo +stable-x86_64-pc-windows-gnu test --workspace --release --exclude topopt-plugin
```

The plugin builds only for wasm and is tested by headless runs (below). Tests run an order of magnitude faster in release, hence `--release`.

| Suite | Checks |
|---|---|
| `topopt-core` (unit, `tests/solver.rs`, `tests/optimize.rs`) | Hex8, indexing, filter, Galerkin, MG-CG against CG, Timoshenko, gradient against finite differences, convergence |
| `topopt-mesh` (`tests/pipeline.rs`, unit) | Sphere, torus, random fields, pinning, clean-up, STL/PLY |
| `topopt-study` (`tests/run.rs`, `tests/spec.rs`, unit) | Voxelization, an end-to-end study, specs, panel options |
| `topopt-solver` (`tests/protocol.rs`) | A job in parts → `serve(once)` → summary, PLY and densities as the plugin reads them |
| `topopt-solver` (`tests/cli.rs`) | The binary: `check`, `example`, `solve` with an STL spec, JSON progress, files, stdin, exit codes 2 and 3; the library `run` in memory |

### Independent checks

```sh
# The optimizer in Python (scipy, assembled matrix); compare with the core's cantilever example
tools/reference/.venv/Scripts/python tools/reference/top3d_ref.py 12 4 4 20

# An STL or PLY through trimesh: watertight, consistent winding, volume
tools/reference/.venv/Scripts/python tools/reference/check_stl.py out/cantilever.ply
```

### Headless runs in SolidCraft

```sh
# The plugin's whole cycle on a test bracket (80×40×20 with a Ø10 hole), no GUI
python -I scripts/topopt_headless.py OUT.json build/topopt.zip APPROVAL RESULTS_DIR 8000 40

# End to end with TopOpt Solver: WORKDIR ZIP [max_elements] [iterations] [use_solver 0|1]
python -I scripts/topopt_solver_e2e.py WORKDIR build/topopt.zip 30000 60 1
```

`topopt_solver_e2e.py` expects a solver already watching `WORKDIR/results` (`topopt-solver watch WORKDIR/results`).

## Installing into the app

See [plugin.md](plugin.md#install): `solidcraft-cli --keep-plugin-approvals exec plugin.install` with `dry_run`, then with the returned `approval`, `"replace": true` and `approve`. **Restart SolidCraft** afterwards: it caches the plugin's buttons.

Driving the app from an agent over MCP:
- The MCP connection may touch only its own folders, so `plugin.install` through MCP is refused; install with `solidcraft-cli --keep-plugin-approvals`.
- Call `plugin.tools` before `execute topopt.*`.
- Entries run by the agent cannot write to a `results` folder outside the agent's file scope; from the GUI everything works.

## Host limits that are easy to trip on

| Limit | What we do |
|---|---|
| One thread in the sandbox | Heavy runs go to TopOpt Solver |
| Normal entry: 30 s and 256 host calls; `"job": true`: 10 min and 1M calls | Optimize and Re-mesh are job entries |
| Fuel is approved at install; we ask for 10 T | Fuel burns faster in the app than headless |
| 1 MiB per plugin file read or write | PLY output, up to 50k triangles; the solver protocol in parts of 1 000 000 bytes |
| A run's file writes are visible only after it ends | Optimize writes the file, Create body inserts it in another run |
| Manifest ≤ 16 KiB | Compact JSON, checked by the generator |
| The host does not scroll panels (only tables do) | Tabs; each page fits |
| A ribbon button's label is the `ribbon` slot's name, else the `menu` slot's, else the entry's | Every entry has a readable `menu` slot and an `x-solidcraft-icon` |
| `forms.open` works only in interactive runs | Panel buttons run the dialog entries, which open their form |
| No `user-storage@0.1` in 0.49.1 | Delete by number, or through a row's double-click |
| No underscores in entry names | `design-body`, `add-support`… |
| `model.mesh` ≤ 8192 vertices / 16384 triangles | The pinning surface is requested at 0.05 mm tolerance, then coarser (0.2, 0.5, 2 mm) until it fits |
| `model.mesh` has no face numbers | With kept faces, one `view.pick` ray per triangle (from just outside its centroid, against its normal) names its face: at most 16 calls of 1024 rays |

## How to add…

- **A plugin setting:** a field in `Settings` or `MeshSettings` (`crates/topopt-study/src/study.rs`), its key in `Study::set_option`, a form field in `gen_topopt_manifest.py`, its value in `panel.rs`. Specs get it for free: they use the same structs with `#[serde(default)]`.
- **A spec region:** a `Region` variant and its branch in `Region::contains` (`spec.rs`), a test in `tests/spec.rs`, a row in [topopt-solver.md §4.3](topopt-solver.md#43-regions-supports-loads-keep_in).
- **A mesh format:** a writer in `topopt-mesh/src/mesh.rs`, its extension in `topopt_solver::mesh_bytes`.
- **A constraint or objective:** a term in `Optimizer::evaluate` with its sensitivity (`topopt-core/src/opt.rs`) and a finite-difference test.
- **A solver-protocol field:** if an older plugin would misread it, bump `PROTOCOL` (or `JOB_VERSION` for the job format).

## Versions

- Plugin: the manifest version (0.12.1).
- Crates and solver: the workspace `version` (0.1.0); the solver reports it in its heartbeat, `version` and the `started` event.
- Plugin ↔ solver protocol: `PROTOCOL = 1`; job format: `JOB_VERSION = 1`.
