# TopOpt CUDA Solver and topopt-cuda

**TopOpt CUDA Solver** (`tools/topopt-cuda-solver`) is TopOpt Solver with the finite-element solves on an NVIDIA GPU. Everything else is the same: commands, spec format, outputs, progress lines, exit codes, and the folder protocol of the SolidCraft plugin. Results agree with the CPU solver to the solver tolerance.

The GPU part is the crate **topopt-cuda** (`crates/topopt-cuda`). It runs topopt-core's multigrid-preconditioned conjugate gradients with CUDA kernels and plugs into any program built on the core through `Optimizer::with_solver` or `RunOptions::solver`.

| | |
|---|---|
| Needs | An NVIDIA GPU of compute capability 5.2 or newer (second-generation Maxwell, Pascal and later) and its driver (CUDA 11.0 / R450 or later). No CUDA toolkit, either to build or to run |
| Without a usable GPU | It solves on the CPU and says why in the run's warnings; `--require-gpu` makes that an error |
| Speed (Quadro P620 against 6 CPU threads) | Linear solves 3–4× faster; whole runs of the examples 1.8–2.1× faster |
| Accuracy | The same iterations, compliance (to 6 digits), mass and topology as the CPU solver on the examples |

## 1. Use

```sh
scripts/build-solver.sh                       # build/solver/topopt-solver.exe and topopt-cuda-solver.exe

build/solver/topopt-cuda-solver devices       # the GPUs the driver sees
build/solver/topopt-cuda-solver solve examples/cantilever.json --out out
build/solver/topopt-cuda-solver watch build/results     # for the SolidCraft plugin
```

With the plugin, start `build/solver/start-topopt-cuda-solver.cmd` instead of `start-topopt-solver.cmd`. The plugin does not care which solver answers; start only one of them on a folder.

All commands and options of [TopOpt Solver](topopt-solver.md) work. The GPU options may appear anywhere on the line:

| Option | Meaning |
|---|---|
| `devices` | List the GPUs: index, name, memory, compute capability |
| `--device N` | The GPU to use (default 0) |
| `--cpu` | Solve on the CPU only, exactly as `topopt-solver` |
| `--require-gpu` | Fail (exit code 1) instead of falling back to the CPU when the GPU cannot be used |
| `--gpu-min-elements N` | Problems with fewer elements stay on the CPU (default 5000; below that, kernel launch latency eats the gain) |

`result.json` names what solved the systems in its `solver` field, for example `"Quadro P620 (CUDA, sm_61)"` or `"CPU, 6 threads"`. A fallback adds a warning such as `solved on the CPU: the GPU could not be used (no NVIDIA driver (nvcuda.dll not found))`. The progress stream says `optimizing on …` when the optimization starts.

## 2. How it works

The optimizer, filter, sensitivities, trimming and meshing stay on the CPU; the GPU does what dominates an iteration on the CPU, the linear solves `K(E)·u = f`.

**Mixed precision with iterative refinement.** Consumer and workstation GPUs run double precision 32–64 times slower than single (the P620: 1/32). So the expensive parts run in single precision: the fine operator, the whole V-cycle and the CG vectors. Dot products still accumulate in double. An outer loop in double precision keeps the answer as accurate as the CPU's:

1. `r = f − K·u` in double (one double-precision product);
2. stop when `‖r‖ ≤ tol·‖f‖`;
3. solve `K·δ = r/‖r‖` with single-precision multigrid-CG to a relative `max(0.5·tol·‖f‖/‖r‖, 10⁻⁴)`;
4. `u += ‖r‖·δ`; repeat.

Two or three refinements per solve are typical. CG needs about 15 % more iterations than the CPU's double-precision CG, because restarting loses the Krylov space.

**The same multigrid as the CPU's** (`mg::Hierarchy`): the same levels, damped Jacobi smoothing with the step from a power iteration, element-wise Galerkin coarse operators, and the coarsest level solved by the CPU's sparse Cholesky. Per optimization iteration:

| Where | What |
|---|---|
| GPU | Fine operator `Σ E_e·k` (matrix-free), its inverse diagonal; the first coarse level's element matrices (`Σ_c E_child·P_cᵀkP_c`, the largest level, never leaves the GPU); smoothers' diagonals and damping (power iteration); every V-cycle except the coarsest solve; CG; refinement |
| CPU | The deeper levels' Galerkin matrices (1/64 of the elements and fewer), uploaded; the coarsest level's factorization, and its solve in each V-cycle (a few thousand unknowns, ≈ 0.1–0.3 ms with the transfers) |

**Kernels** (`crates/topopt-cuda/kernels/lib.rs`, 19 of them):
- **One thread per node.** Each thread gathers its three rows from the up to eight elements around it. Nothing is written twice and no atomics are needed, so every result is deterministic.
- **The element matrix as a kernel parameter.** `ke` goes in by value (2304 bytes in single precision), so it lives in the constant bank. The loops over corners are unrolled at compile time, so the multiply-adds take its entries as operands with no loads. In double precision it is 4.5 KB, over the 4 KB parameter limit, and is read through the read-only cache instead.
- **The 27 neighbours loaded once.** The single-precision fine operator gathers the neighbouring nodes up front, though up to eight elements share each one. In double precision and on coarse levels this cost more registers than it saved, so those kernels load where needed.
- **Read-only cache.** All operator inputs go through `ld.global.nc` (CUDA's `__ldg`). On Pascal, plain global loads are cached in L2 only.
- **Entry-major coarse matrices** (`mats[q·nel + e]`): neighbouring threads read neighbouring addresses. With element-major storage, a coarse product cost 9 ms instead of 1.1 ms.
- **Reductions** go through a per-block scratch area and a fixed number of blocks, so they are deterministic.

**The driver is loaded at run time** (`nvcuda.dll` / `libcuda.so.1`) through a small safe layer (`driver.rs`: context, module, buffers, launches). The kernels are checked in as PTX (`src/kernels.ptx`), and the driver compiles them for the GPU at hand when the solver starts (and caches the result).

## 3. Accuracy

`cargo test -p topopt-cuda -p topopt-cuda-solver --release` (the GPU tests pass with a "skipped" note on machines without a GPU):

| Check | Result |
|---|---|
| Fine operator against `fem::apply` (random moduli 10⁻⁹…1, supports) | double 1.4·10⁻¹⁶, single 5.8·10⁻⁸ (relative) |
| One V-cycle against `Hierarchy::vcycle` | 1.0·10⁻⁷ |
| Solve at `tol` 10⁻⁸, 32×16×16 | residual 4·10⁻⁹ (checked on the CPU in double); 1.1·10⁻¹⁰ from the CPU's answer; 32 iterations against 27 |
| 25 optimization iterations, 24×12×12 | compliance within 10⁻⁴ every iteration; densities within 5.7·10⁻⁷ |
| The CLI on `examples/cantilever.json`, 6 iterations | the CPU's compliance to 10⁻⁵, the same iterations and triangles |
| The PTX against its source | the hash of `kernels/lib.rs` recorded in `kernels.ptx` |

Whole runs of the examples (80 000 elements, defaults):

| Example | CPU, 6 threads | Quadro P620 | Result |
|---|---|---|---|
| cantilever | 144 s, 169 iterations, C = 105.989 | **79 s**, 169 iterations, C = 105.989 | 227.3 g, genus 31 on both |
| bridge | 231 s, 200 iterations, C = 136.410 | **108 s**, 200 iterations, C = 136.409 | 506.7 g, genus 21 on both |

## 4. Performance

`cargo run --release -p topopt-cuda --example bench` takes one optimization iteration (update and a warm-started solve to 10⁻⁶) on a cantilever with a converged-looking design (solid, void at 10⁻⁹, some grey). Intel i5-9400F (6 threads) against a Quadro P620 (512 cores, 96 GB/s, 2 GB):

| Grid | Elements | CPU update + solve | GPU update + solve | Speed-up |
|---|---|---|---|---|
| 32×16×16 | 8 192 | 5 + 68 ms | 4 + 24 ms | 2.6× |
| 40×20×20 | 16 000 | 9 + 125 ms | 7 + 43 ms | 2.7× |
| 64×32×32 | 65 536 | 25 + 507 ms | 17 + 137 ms | 3.5× |
| 96×48×48 | 221 184 | 85 + 2243 ms | 53 + 546 ms | 3.9× |
| 128×64×64 | 524 288 | 193 + 6482 ms | 120 + 1544 ms | 4.0× |

`--example gpu_profile NX NY NZ` breaks the GPU work down. At 96×48×48:

| Operation | ms |
|---|---|
| fine product, single | 1.5 |
| fine product, double | 6.5 |
| first coarse level product (27 648 elements) | 1.1 |
| V-cycle | 12.7 |

The P620 is an entry-level card. A GPU with more memory bandwidth gains in proportion, since the products are bound by memory and loads, not arithmetic.

Whole runs gain less than the solves: of the 0.33 s an iteration takes at 80 000 elements, the GPU solve is about 0.13 s. The rest is CPU work shared with the CPU solver: the optimality-criteria bisection, which filters the design up to 200 times, plus the sensitivities. Meshing adds about 24 s per run. Those are the next things to speed up.

**Memory.** About 550 bytes per element: the fine level's vectors, single and double, plus the first coarse level's matrices (2304 bytes per coarse element). 400 000 elements, the study limit, need about 220 MB. The solver checks free memory first; when it is short, it falls back to the CPU with a warning saying how much was needed.

## 5. Rust API

```toml
[dependencies]
topopt-cuda = { path = "…/crates/topopt-cuda" }
```

```rust
use topopt_cuda::{CudaOptions, factory};

// In a pipeline: TopOpt Solver's run() with the GPU.
let opts = topopt_solver::RunOptions { solver: Some(factory(CudaOptions::default())), ..Default::default() };
let out = topopt_solver::run(&input, &opts, &mut |_| {})?;
println!("solved on {}", out.result.solver);

// With the core alone.
let make = factory(CudaOptions { device: 0, min_elements: 0, fallback: false });
let mut opt = topopt_core::Optimizer::with_solver(problem, settings, Some(&make))?;
```

| Item | Purpose |
|---|---|
| `factory(CudaOptions) -> SolverFactory` | Makes a `CudaSolver` per problem, or a CPU solver (small problems, or fallback). Keeps the CUDA context and compiled kernels between problems |
| `CudaOptions { device, min_elements, fallback }` | GPU index; problems smaller than `min_elements` stay on the CPU (default 5000); with `fallback`, use the CPU when the GPU cannot be used (default `true`) |
| `CudaSolver::new(gpu, grid, fixed, mg)` | The solver itself (`topopt_core::LinearSolver`) |
| `open(device) -> Arc<Gpu>`, `devices()` | Open a GPU with the kernels loaded; list the GPUs and the driver version |
| `memory_needed(grid, coarse_grids)` | Bytes of GPU memory a problem needs |

Other back-ends plug in the same way: implement `topopt_core::LinearSolver` (`update(ke, emod)`, `solve(…)`, `levels()`, `name()`, optionally `note()`) and pass a `SolverFactory`. `topopt-study` takes it through `Hooks::solver`, and TopOpt Solver through `RunOptions::solver`, `serve_with` and `remesh_with`. Its command line is the library function `topopt_solver::cli::main_with(App, args)`: `topopt-cuda-solver` is that function with a GPU factory and a few more options.

## 6. Rebuilding the kernels

The kernels are Rust (`#![no_std]`, `nvptx64-nvidia-cuda`) in `crates/topopt-cuda/kernels`, a separate crate outside the workspace. Building them needs a nightly toolchain; building and running the solver does not.

```sh
rustup toolchain install nightly --profile minimal --component llvm-tools,llvm-bitcode-linker,rust-src --target nvptx64-nvidia-cuda
scripts/build-ptx.sh            # PTX_TOOLCHAIN=nightly-x86_64-pc-windows-gnu to choose one
cargo test -p topopt-cuda --release
```

- **Target retargeting.** Rust's nvptx target emits PTX for sm_70 (Volta) and later. The kernels use no instruction newer than sm_32, so `build-ptx.sh` retargets the PTX to sm_52. That lets Maxwell and Pascal GPUs load it, and the driver still checks every instruction when it compiles the PTX.
- **Source hash.** The PTX's second line records an FNV-1a hash of `kernels/lib.rs`. The test `ptx_matches_the_kernel_source` fails when the source changed without a rebuild.
- **Nightly features in use:**
  - `abi_ptx`: kernels;
  - `stdarch_nvptx`: thread indices and barriers;
  - `core_float_math`: fused multiply-add;
  - `asm_experimental_arch`: `ld.global.nc`.

## 7. Limits

- NVIDIA only: CUDA through the driver API. One GPU per run.
- Tested on a Pascal GPU (Quadro P620, sm_61, driver CUDA 13.0). Newer GPUs compile the same PTX.
- The plugin itself cannot use a GPU (the SolidCraft sandbox has no GPU compute). It gets the GPU through the solver protocol, like the CPU cores.
- On Windows (WDDM) every kernel launch costs microseconds. Small problems therefore stay on the CPU (`--gpu-min-elements`), and the coarsest levels are better on the CPU anyway.
