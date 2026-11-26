//! The linear solver of `topopt-core` on an NVIDIA GPU: the same multigrid-preconditioned
//! conjugate gradients as [`topopt_core::CpuSolver`], run by CUDA kernels.
//!
//! - **Mixed precision.** The fine operator, the V-cycle and the CG vectors are single precision
//!   (consumer GPUs run double precision 32–64 times slower); dot products accumulate in double.
//!   An outer iterative refinement in double precision keeps the result as accurate as the CPU's:
//!   `r = f − K·u` in double, solve `K·δ = r` in single, `u += δ`, until `‖r‖ ≤ tol·‖f‖`.
//! - **What stays on the CPU.** The coarse levels' Galerkin matrices (cheap, once per
//!   optimization iteration) and the factorized coarsest level (a few thousand dofs) come from
//!   [`topopt_core::mg::Hierarchy`]; the smoothers' diagonals and damping (power iteration) are
//!   computed on the GPU.
//! - **No toolkit.** The kernels (`kernels/`, Rust for `nvptx64-nvidia-cuda`) are checked in as PTX
//!   and loaded through the driver, which is all a machine needs.
//!
//! [`factory`] plugs it into [`topopt_core::Optimizer::with_solver`]; it falls back to the CPU
//! when there is no usable GPU, or for problems too small to gain from one.

pub mod driver;

use driver::{Arg, Buffer, DeviceInfo, Driver, Gpu};
use std::sync::Arc;
use topopt_core::backend::{CpuSolver, LinearSolver, SolverFactory};
use topopt_core::ke::Ke;
use topopt_core::mg::{Hierarchy, MgSettings};
use topopt_core::pcg::PcgStats;
use topopt_core::Grid;

/// The kernels, built by `scripts/build-ptx.sh`.
pub const PTX: &str = include_str!("kernels.ptx");

const NAMES: [&str; 19] = [
    "apply_fine_f32",
    "apply_fine_f64",
    "apply_own_f32",
    "inv_diag_fine",
    "inv_diag_own",
    "jacobi",
    "restrict_residual",
    "prolong_add",
    "dot_f32",
    "dot_f64",
    "sum_partials",
    "cg_xr",
    "cg_p",
    "residual_f64",
    "to_f32",
    "add_f32_to_f64",
    "scale_f32",
    "scale_inv_diag",
    "galerkin_first",
];
const APPLY_FINE_F32: usize = 0;
const APPLY_FINE_F64: usize = 1;
const APPLY_OWN_F32: usize = 2;
const INV_DIAG_FINE: usize = 3;
const INV_DIAG_OWN: usize = 4;
const JACOBI: usize = 5;
const RESTRICT: usize = 6;
const PROLONG: usize = 7;
const DOT_F32: usize = 8;
const DOT_F64: usize = 9;
const SUM_PARTIALS: usize = 10;
const CG_XR: usize = 11;
const CG_P: usize = 12;
const RESIDUAL_F64: usize = 13;
const TO_F32: usize = 14;
const ADD_F32_TO_F64: usize = 15;
const SCALE_F32: usize = 16;
const SCALE_INV_DIAG: usize = 17;
const GALERKIN_FIRST: usize = 18;

const THREADS: u32 = 256;
const REDUCE_BLOCKS: u32 = 128;

/// Device scalars: `r·z` of two consecutive CG iterations, `p·q`, `r·r`, and a slot for host dots.
const RZ: [u32; 2] = [0, 1];
const PQ: u32 = 2;
const RR: u32 = 3;
const HOST: u32 = 4;

/// The inner (single-precision) CG reduces its residual at most this far per refinement: below,
/// rounding makes further iterations useless.
const INNER_FLOOR: f64 = 1e-4;
/// Refinements per solve before giving up.
const MAX_REFINEMENTS: usize = 30;

fn blocks(n: usize) -> u32 {
    n.div_ceil(THREADS as usize) as u32
}

/// The driver's version and its GPUs.
pub fn devices() -> Result<(i32, Vec<DeviceInfo>), String> {
    let d = Driver::load()?;
    Ok((d.version, d.devices()?))
}

/// "13.0" for 13000.
pub fn version_string(v: i32) -> String {
    format!("{}.{}", v / 1000, (v % 1000) / 10)
}

/// How [`factory`] chooses.
#[derive(Clone, Debug, PartialEq)]
pub struct CudaOptions {
    /// The GPU, as numbered by the driver.
    pub device: usize,
    /// Smaller problems are solved on the CPU, where the GPU's per-kernel latency would dominate.
    pub min_elements: usize,
    /// Use the CPU (with a note) when the GPU cannot be used; otherwise that is an error.
    pub fallback: bool,
}

impl Default for CudaOptions {
    fn default() -> Self {
        CudaOptions { device: 0, min_elements: 5_000, fallback: true }
    }
}

/// A [`SolverFactory`] that makes [`CudaSolver`]s.
pub fn factory(opts: CudaOptions) -> SolverFactory {
    let gpu: std::sync::Mutex<Option<Arc<Gpu>>> = std::sync::Mutex::new(None);
    Arc::new(move |grid: &Grid, fixed: &[bool], mg: &MgSettings| -> Result<Box<dyn LinearSolver>, String> {
        if grid.nel() < opts.min_elements {
            return Ok(Box::new(CpuSolver::new(*grid, fixed.to_vec(), *mg)));
        }
        // The context and the JIT-compiled kernels are kept for the next solver.
        let made = (|| {
            let mut cached = gpu.lock().map_err(|_| "GPU lock poisoned".to_string())?;
            if cached.is_none() {
                *cached = Some(open(opts.device)?);
            }
            CudaSolver::new(cached.clone().unwrap(), *grid, fixed.to_vec(), *mg)
        })();
        match made {
            Ok(s) => Ok(Box::new(s)),
            Err(e) if opts.fallback => Ok(Box::new(Fallback { cpu: CpuSolver::new(*grid, fixed.to_vec(), *mg), reason: e })),
            Err(e) => Err(e),
        }
    })
}

/// Opens GPU `device` with the kernels loaded.
pub fn open(device: usize) -> Result<Arc<Gpu>, String> {
    let driver = Driver::load()?;
    let count = driver.devices()?.len();
    if device >= count {
        return Err(if count == 0 { "no CUDA GPU found".into() } else { format!("no CUDA GPU {device} (there are {count})") });
    }
    Gpu::open(driver, device, PTX, &NAMES)
}

/// The CPU solver standing in for a GPU that could not be used.
struct Fallback {
    cpu: CpuSolver,
    reason: String,
}

impl LinearSolver for Fallback {
    fn update(&mut self, ke: &Ke, emod: &[f64]) -> Result<(), String> {
        self.cpu.update(ke, emod)
    }
    fn solve(&mut self, ke: &Ke, emod: &[f64], f: &[f64], u: &mut [f64], tol: f64, max_iter: usize) -> Result<PcgStats, String> {
        self.cpu.solve(ke, emod, f, u, tol, max_iter)
    }
    fn levels(&self) -> Vec<Grid> {
        self.cpu.levels()
    }
    fn name(&self) -> String {
        self.cpu.name()
    }
    fn note(&self) -> Option<String> {
        Some(format!("solved on the CPU: the GPU could not be used ({})", self.reason))
    }
}

struct Level {
    grid: Grid,
    nnodes: usize,
    ndof: usize,
    fixed: Vec<bool>,
    fix: Buffer<u8>,
    invd: Buffer<f32>,
    /// Right-hand side and solution of this level in a V-cycle (unused on the fine level, which
    /// works on CG's vectors).
    b: Buffer<f32>,
    x: Buffer<f32>,
    ax: Buffer<f32>,
    /// Power-iteration vector, kept as the warm start.
    eig: Buffer<f32>,
    eig_ready: bool,
    omega: f32,
    /// Element matrices (coarse levels).
    mats: Option<Buffer<f32>>,
}

/// Multigrid-preconditioned CG on a GPU, mixed precision (see the crate documentation).
pub struct CudaSolver {
    gpu: Arc<Gpu>,
    mg: MgSettings,
    hierarchy: Hierarchy,
    /// The fine level, then every coarse level; the last one is solved by the CPU's factorization.
    levels: Vec<Level>,
    ke32: Buffer<f32>,
    /// The same, passed by value to the fine operator (it lands in the constant bank).
    ke32_param: Box<[f32; 576]>,
    ke64: Buffer<f64>,
    /// The first coarse level's eight child Galerkin matrices.
    child_galerkin: Buffer<f32>,
    emod32: Buffer<f32>,
    emod64: Buffer<f64>,
    // Inner CG (single precision): solution, residual, preconditioned residual, direction, K·p.
    cx: Buffer<f32>,
    r: Buffer<f32>,
    z: Buffer<f32>,
    p: Buffer<f32>,
    q: Buffer<f32>,
    // Outer refinement (double precision).
    f: Buffer<f64>,
    u: Buffer<f64>,
    ax64: Buffer<f64>,
    r64: Buffer<f64>,
    scratch: Buffer<f64>,
    partial: Buffer<f64>,
    scal: Buffer<f64>,
    /// Milliseconds of the last update's phases: coarse operators on the CPU, their upload, the
    /// smoothers on the GPU.
    #[doc(hidden)]
    pub update_ms: [f64; 3],
}

/// Bytes of GPU memory a problem needs.
pub fn memory_needed(grid: &Grid, coarse: &[Grid]) -> usize {
    let fine = grid.ndof() * (8 * 4 + 4 * 8) + grid.nnodes() + grid.nel() * 12;
    let levels: usize = coarse.iter().map(|g| g.ndof() * 5 * 4 + g.nnodes() + g.nel() * 576 * 4).sum();
    fine + levels + (REDUCE_BLOCKS * THREADS) as usize * 8 + (1 << 20)
}

impl CudaSolver {
    pub fn new(gpu: Arc<Gpu>, grid: Grid, fixed: Vec<bool>, mg: MgSettings) -> Result<Self, String> {
        gpu.bind()?;
        let hierarchy = Hierarchy::new(grid, fixed.clone(), mg);
        let coarse: Vec<Grid> = hierarchy.level_grids()[1..].to_vec();
        let need = memory_needed(&grid, &coarse);
        let (free, _) = gpu.memory()?;
        if need > free {
            return Err(format!(
                "{} elements need about {} MB of GPU memory and {} MB are free",
                grid.nel(),
                need.div_ceil(1 << 20),
                free >> 20
            ));
        }
        let mut levels = Vec::new();
        for l in 0..=coarse.len() {
            let (g, fx) = if l == 0 { (grid, fixed.clone()) } else { let (g, f, _) = hierarchy.coarse_level(l - 1); (g, f.to_vec()) };
            let masks: Vec<u8> = fx.chunks(3).map(|c| c[0] as u8 | (c[1] as u8) << 1 | (c[2] as u8) << 2).collect();
            let n = g.ndof();
            levels.push(Level {
                grid: g,
                nnodes: g.nnodes(),
                ndof: n,
                fix: Buffer::from_slice(&gpu, &masks)?,
                invd: Buffer::new(&gpu, n)?,
                b: Buffer::new(&gpu, if l == 0 { 0 } else { n })?,
                x: Buffer::new(&gpu, if l == 0 { 0 } else { n })?,
                ax: Buffer::new(&gpu, n)?,
                eig: Buffer::new(&gpu, n)?,
                eig_ready: false,
                omega: 0.0,
                mats: if l == 0 { None } else { Some(Buffer::new(&gpu, 576 * g.nel())?) },
                fixed: fx,
            });
        }
        let n = grid.ndof();
        let nel = grid.nel();
        Ok(CudaSolver {
            mg,
            hierarchy,
            levels,
            ke32: Buffer::new(&gpu, 576)?,
            ke32_param: Box::new([0.0; 576]),
            ke64: Buffer::new(&gpu, 576)?,
            child_galerkin: Buffer::new(&gpu, 8 * 576)?,
            emod32: Buffer::new(&gpu, nel)?,
            emod64: Buffer::new(&gpu, nel)?,
            cx: Buffer::new(&gpu, n)?,
            r: Buffer::new(&gpu, n)?,
            z: Buffer::new(&gpu, n)?,
            p: Buffer::new(&gpu, n)?,
            q: Buffer::new(&gpu, n)?,
            f: Buffer::new(&gpu, n)?,
            u: Buffer::new(&gpu, n)?,
            ax64: Buffer::new(&gpu, n)?,
            r64: Buffer::new(&gpu, n)?,
            scratch: Buffer::new(&gpu, (REDUCE_BLOCKS * THREADS) as usize)?,
            partial: Buffer::new(&gpu, REDUCE_BLOCKS as usize)?,
            scal: Buffer::new(&gpu, 8)?,
            update_ms: [0.0; 3],
            gpu,
        })
    }

    pub fn device(&self) -> &DeviceInfo {
        &self.gpu.info
    }

    fn dims(g: &Grid) -> [Arg; 3] {
        [Arg::U32(g.nx as u32), Arg::U32(g.ny as u32), Arg::U32(g.nz as u32)]
    }

    /// `y = A_l·x` on level `l` (single precision).
    fn apply(&self, l: usize, x: &Buffer<f32>, y: &Buffer<f32>) -> Result<(), String> {
        let lv = &self.levels[l];
        let [nx, ny, nz] = Self::dims(&lv.grid);
        match &lv.mats {
            None => {
                let ke = Arg::Bytes(self.ke32_param.as_ptr() as *const std::ffi::c_void);
                self.gpu.launch(APPLY_FINE_F32, blocks(lv.nnodes), THREADS, &[nx, ny, nz, ke, self.emod32.arg(), lv.fix.arg(), x.arg(), y.arg()])
            }
            Some(m) => self.gpu.launch(APPLY_OWN_F32, blocks(lv.nnodes), THREADS, &[nx, ny, nz, m.arg(), lv.fix.arg(), x.arg(), y.arg()]),
        }
    }

    /// `y = K·x` on the fine level in double precision.
    fn apply64(&self, x: &Buffer<f64>, y: &Buffer<f64>) -> Result<(), String> {
        let lv = &self.levels[0];
        let [nx, ny, nz] = Self::dims(&lv.grid);
        self.gpu.launch(APPLY_FINE_F64, blocks(lv.nnodes), THREADS, &[nx, ny, nz, self.ke64.arg(), self.emod64.arg(), lv.fix.arg(), x.arg(), y.arg()])
    }

    /// `scal[slot] = a·b`.
    fn dot<T: Copy>(&self, kernel: usize, n: usize, a: &Buffer<T>, b: &Buffer<T>, slot: u32) -> Result<(), String> {
        self.gpu.launch(kernel, REDUCE_BLOCKS, THREADS, &[Arg::U32(n as u32), a.arg(), b.arg(), self.scratch.arg(), self.partial.arg()])?;
        self.gpu.launch(SUM_PARTIALS, 1, THREADS, &[Arg::U32(REDUCE_BLOCKS), self.partial.arg(), self.scratch.arg(), self.scal.arg(), Arg::U32(slot)])
    }

    fn host_dot32(&self, n: usize, a: &Buffer<f32>, b: &Buffer<f32>) -> Result<f64, String> {
        self.dot(DOT_F32, n, a, b, HOST)?;
        self.scal.get(HOST as usize)
    }

    fn host_dot64(&self, a: &Buffer<f64>, b: &Buffer<f64>) -> Result<f64, String> {
        self.dot(DOT_F64, a.len(), a, b, HOST)?;
        self.scal.get(HOST as usize)
    }

    fn jacobi(&self, l: usize, b: &Buffer<f32>, x: &Buffer<f32>, first: bool) -> Result<(), String> {
        let lv = &self.levels[l];
        self.gpu.launch(
            JACOBI,
            blocks(lv.ndof),
            THREADS,
            &[Arg::U32(lv.ndof as u32), Arg::F32(lv.omega), lv.invd.arg(), b.arg(), lv.ax.arg(), x.arg(), Arg::U32(first as u32)],
        )
    }

    /// `x ≈ A_l⁻¹·b`: one V-cycle from zero (as [`Hierarchy::vcycle`]).
    fn cycle(&self, l: usize, b: &Buffer<f32>, x: &Buffer<f32>) -> Result<(), String> {
        if l + 1 == self.levels.len() {
            return self.coarse_solve(l, b, x);
        }
        let lv = &self.levels[l];
        let steps = self.mg.smooth_steps;
        let smooth = |first_zero: bool| -> Result<(), String> {
            for s in 0..steps {
                if first_zero && s == 0 {
                    self.jacobi(l, b, x, true)?;
                } else {
                    self.apply(l, x, &lv.ax)?;
                    self.jacobi(l, b, x, false)?;
                }
            }
            Ok(())
        };
        if steps == 0 {
            x.zero()?;
        }
        smooth(true)?;
        self.apply(l, x, &lv.ax)?;
        let next = &self.levels[l + 1];
        let [cnx, cny, cnz] = Self::dims(&next.grid);
        let [fnx, fny, fnz] = Self::dims(&lv.grid);
        self.gpu.launch(RESTRICT, blocks(next.nnodes), THREADS, &[cnx, cny, cnz, next.fix.arg(), fnx, fny, fnz, lv.fix.arg(), b.arg(), lv.ax.arg(), next.b.arg()])?;
        self.cycle(l + 1, &next.b, &next.x)?;
        self.gpu.launch(PROLONG, blocks(lv.nnodes), THREADS, &[cnx, cny, cnz, next.fix.arg(), fnx, fny, fnz, lv.fix.arg(), next.x.arg(), x.arg()])?;
        smooth(false)
    }

    /// The coarsest level on the CPU: download, factorized solve, upload.
    fn coarse_solve(&self, l: usize, b: &Buffer<f32>, x: &Buffer<f32>) -> Result<(), String> {
        let n = self.levels[l].ndof;
        let mut b32 = vec![0.0f32; n];
        b.download(&mut b32)?;
        let b64: Vec<f64> = b32.iter().map(|&v| v as f64).collect();
        let mut x64 = vec![0.0; n];
        self.hierarchy.coarse_solve(&b64, &mut x64);
        let x32: Vec<f32> = x64.iter().map(|&v| v as f32).collect();
        x.upload(&x32)
    }

    /// Largest eigenvalue of `D⁻¹·A` on level `l` by power iteration, as the CPU's `lambda_max`.
    fn lambda_max(&mut self, l: usize) -> Result<f64, String> {
        let steps = if self.levels[l].eig_ready { self.mg.power_steps.max(1) } else { 12 };
        if !self.levels[l].eig_ready {
            let lv = &self.levels[l];
            let v: Vec<f32> = (0..lv.ndof).map(|i| if lv.fixed[i] { 0.0 } else { 1.0 + ((i * 7919) % 13) as f32 / 13.0 }).collect();
            lv.eig.upload(&v)?;
            self.levels[l].eig_ready = true;
        }
        let mut lambda = 1.0;
        for _ in 0..steps {
            let lv = &self.levels[l];
            let n = lv.ndof;
            let norm = self.host_dot32(n, &lv.eig, &lv.eig)?.sqrt();
            if norm == 0.0 || !norm.is_finite() {
                return Ok(1.0);
            }
            self.gpu.launch(SCALE_F32, blocks(n), THREADS, &[Arg::U32(n as u32), Arg::F32((1.0 / norm) as f32), lv.eig.arg()])?;
            self.apply(l, &lv.eig, &lv.ax)?;
            self.gpu.launch(SCALE_INV_DIAG, blocks(lv.nnodes), THREADS, &[Arg::U32(lv.nnodes as u32), lv.fix.arg(), lv.invd.arg(), lv.ax.arg()])?;
            lambda = self.host_dot32(n, &lv.eig, &lv.ax)?;
            let lv = &mut self.levels[l];
            std::mem::swap(&mut lv.eig, &mut lv.ax);
        }
        Ok((lambda * 1.1).max(1.0))
    }

    /// Inner CG in single precision on `K·cx = r` from `cx = 0`, until `‖r‖` falls by `tol` or
    /// `max_iter` iterations. Returns the iterations.
    fn pcg32(&self, tol: f64, max_iter: usize) -> Result<usize, String> {
        let n = self.levels[0].ndof;
        let nu = Arg::U32(n as u32);
        self.cx.zero()?;
        let r0 = self.host_dot32(n, &self.r, &self.r)?.sqrt();
        if r0 == 0.0 {
            return Ok(0);
        }
        self.cycle(0, &self.r, &self.z)?;
        self.dot(DOT_F32, n, &self.r, &self.z, RZ[0])?;
        self.p.copy_from(&self.z)?;
        let mut cur = 0;
        let mut s = [0.0f64; 8];
        for it in 1..=max_iter {
            self.apply(0, &self.p, &self.q)?;
            self.dot(DOT_F32, n, &self.p, &self.q, PQ)?;
            self.gpu.launch(CG_XR, blocks(n), THREADS, &[nu, self.scal.arg(), Arg::U32(RZ[cur]), Arg::U32(PQ), self.cx.arg(), self.r.arg(), self.p.arg(), self.q.arg()])?;
            self.dot(DOT_F32, n, &self.r, &self.r, RR)?;
            self.scal.download_prefix(&mut s[..4])?;
            let (pq, rr) = (s[PQ as usize], s[RR as usize]);
            if !(pq > 0.0) || !rr.is_finite() || rr.sqrt() <= tol * r0 {
                return Ok(it);
            }
            self.cycle(0, &self.r, &self.z)?;
            let next = 1 - cur;
            self.dot(DOT_F32, n, &self.r, &self.z, RZ[next])?;
            self.gpu.launch(CG_P, blocks(n), THREADS, &[nu, self.scal.arg(), Arg::U32(RZ[next]), Arg::U32(RZ[cur]), self.p.arg(), self.z.arg()])?;
            cur = next;
        }
        Ok(max_iter)
    }

    /// Diagnostics: one V-cycle of `r` on the GPU, as [`Hierarchy::vcycle`] (after an update).
    #[doc(hidden)]
    pub fn vcycle(&self, r: &[f64]) -> Result<Vec<f64>, String> {
        self.gpu.bind()?;
        self.r.upload(&r.iter().map(|&v| v as f32).collect::<Vec<_>>())?;
        self.cycle(0, &self.r, &self.z)?;
        let mut z = vec![0.0f32; r.len()];
        self.z.download(&mut z)?;
        Ok(z.into_iter().map(|v| v as f64).collect())
    }

    /// Diagnostics: `K·x` on the GPU in single (`double = false`) or double precision.
    #[doc(hidden)]
    pub fn apply_fine(&self, x: &[f64], double: bool) -> Result<Vec<f64>, String> {
        self.gpu.bind()?;
        if double {
            self.u.upload(x)?;
            self.apply64(&self.u, &self.ax64)?;
            let mut y = vec![0.0; x.len()];
            self.ax64.download(&mut y)?;
            Ok(y)
        } else {
            self.p.upload(&x.iter().map(|&v| v as f32).collect::<Vec<_>>())?;
            self.apply(0, &self.p, &self.q)?;
            let mut y = vec![0.0f32; x.len()];
            self.q.download(&mut y)?;
            Ok(y.into_iter().map(|v| v as f64).collect())
        }
    }

    /// Diagnostics: milliseconds per operation (after an update), each timed over `reps` runs.
    #[doc(hidden)]
    pub fn profile(&self, reps: usize) -> Result<Vec<(String, f64)>, String> {
        self.gpu.bind()?;
        let time = |f: &dyn Fn() -> Result<(), String>| -> Result<f64, String> {
            f()?;
            self.gpu.synchronize()?;
            let t = std::time::Instant::now();
            for _ in 0..reps {
                f()?;
            }
            self.gpu.synchronize()?;
            Ok(t.elapsed().as_secs_f64() * 1e3 / reps as f64)
        };
        let mut out = Vec::new();
        for (k, name) in NAMES.iter().enumerate().take(3) {
            let (regs, local) = self.gpu.kernel_resources(k)?;
            out.push((format!("{name}: {regs} registers, {local} B local"), 0.0));
        }
        let n = self.levels[0].ndof;
        out.push(("apply fine f32".into(), time(&|| self.apply(0, &self.p, &self.q))?));
        out.push(("apply fine f64".into(), time(&|| self.apply64(&self.u, &self.ax64))?));
        for l in 1..self.levels.len() {
            let lv = &self.levels[l];
            out.push((format!("apply level {l} ({} el)", lv.grid.nel()), time(&|| self.apply(l, &lv.x, &lv.ax))?));
        }
        out.push(("dot f32".into(), time(&|| self.dot(DOT_F32, n, &self.r, &self.z, HOST))?));
        out.push(("dot f32 + read".into(), time(&|| self.host_dot32(n, &self.r, &self.z).map(|_| ()))?));
        out.push(("jacobi".into(), time(&|| self.jacobi(0, &self.r, &self.z, false))?));
        let last = self.levels.len() - 1;
        out.push(("coarse solve".into(), time(&|| self.coarse_solve(last, &self.levels[last].b, &self.levels[last].x))?));
        out.push(("V-cycle".into(), time(&|| self.cycle(0, &self.r, &self.z))?));
        Ok(out)
    }

    /// Diagnostics: the smoothers' damping per level (after an update).
    #[doc(hidden)]
    pub fn omegas(&self) -> Vec<f32> {
        self.levels.iter().map(|l| l.omega).collect()
    }
}

impl LinearSolver for CudaSolver {
    fn update(&mut self, ke: &Ke, emod: &[f64]) -> Result<(), String> {
        self.gpu.bind()?;
        *self.ke32_param = ke.map(|v| v as f32);
        self.ke32.upload(&self.ke32_param[..])?;
        self.ke64.upload(ke)?;
        self.emod32.upload(&emod.iter().map(|&v| v as f32).collect::<Vec<_>>())?;
        self.emod64.upload(emod)?;
        let t = std::time::Instant::now();
        self.hierarchy.update_coarse(ke, emod)?;
        self.update_ms[0] = t.elapsed().as_secs_f64() * 1e3;
        let t = std::time::Instant::now();
        for l in 1..self.levels.len() {
            let lv = &self.levels[l];
            let dst = lv.mats.as_ref().unwrap();
            if l == 1 {
                // The largest coarse level is formed on the GPU from the fine moduli.
                let g = self.hierarchy.child_galerkin().ok_or("no child Galerkin matrices")?;
                let g32: Vec<f32> = g.as_flattened().iter().map(|&v| v as f32).collect();
                self.child_galerkin.upload(&g32)?;
                let [cnx, cny, cnz] = Self::dims(&lv.grid);
                self.gpu.launch(GALERKIN_FIRST, blocks(lv.grid.nel()), THREADS, &[cnx, cny, cnz, self.child_galerkin.arg(), self.emod32.arg(), dst.arg()])?;
                continue;
            }
            // Entry-major on the GPU (`[q·nel + e]`), for coalesced reads.
            let (_, _, mats) = self.hierarchy.coarse_level(l - 1);
            let nel = mats.len();
            let mut soa = vec![0.0f32; 576 * nel];
            for (e, m) in mats.iter().enumerate() {
                for (q, &v) in m.iter().enumerate() {
                    soa[q * nel + e] = v;
                }
            }
            dst.upload(&soa)?;
        }
        self.update_ms[1] = t.elapsed().as_secs_f64() * 1e3;
        let t = std::time::Instant::now();
        // Smoothers on every level but the factorized one.
        for l in 0..self.levels.len() - 1 {
            let lv = &self.levels[l];
            let [nx, ny, nz] = Self::dims(&lv.grid);
            match &lv.mats {
                None => self.gpu.launch(INV_DIAG_FINE, blocks(lv.nnodes), THREADS, &[nx, ny, nz, self.ke32.arg(), self.emod32.arg(), lv.fix.arg(), lv.invd.arg()])?,
                Some(m) => self.gpu.launch(INV_DIAG_OWN, blocks(lv.nnodes), THREADS, &[nx, ny, nz, m.arg(), lv.fix.arg(), lv.invd.arg()])?,
            }
            let lambda = self.lambda_max(l)?;
            self.levels[l].omega = (self.mg.omega * 2.0 / lambda) as f32;
        }
        self.update_ms[2] = t.elapsed().as_secs_f64() * 1e3;
        Ok(())
    }

    fn solve(&mut self, _ke: &Ke, _emod: &[f64], f: &[f64], u: &mut [f64], tol: f64, max_iter: usize) -> Result<PcgStats, String> {
        self.gpu.bind()?;
        let n = f.len();
        self.f.upload(f)?;
        let fnorm = self.host_dot64(&self.f, &self.f)?.sqrt();
        if fnorm == 0.0 {
            u.iter_mut().for_each(|v| *v = 0.0);
            return Ok(PcgStats { iterations: 0, rel_residual: 0.0, converged: true });
        }
        self.u.upload(u)?;
        let nu = Arg::U32(n as u32);
        let mut iterations = 0;
        let mut rel = f64::INFINITY;
        let mut converged = false;
        for refinement in 0..=MAX_REFINEMENTS {
            self.apply64(&self.u, &self.ax64)?;
            self.gpu.launch(RESIDUAL_F64, blocks(n), THREADS, &[nu, self.f.arg(), self.ax64.arg(), self.r64.arg()])?;
            let s = self.host_dot64(&self.r64, &self.r64)?.sqrt();
            let prev = rel;
            rel = s / fnorm;
            if rel <= tol {
                converged = true;
                break;
            }
            // Out of iterations, or single precision stopped making progress.
            if iterations >= max_iter || refinement == MAX_REFINEMENTS || !rel.is_finite() || rel > 0.9 * prev {
                break;
            }
            self.gpu.launch(TO_F32, blocks(n), THREADS, &[nu, self.r64.arg(), Arg::F64(1.0 / s), self.r.arg()])?;
            let inner_tol = (0.5 * tol / rel).clamp(INNER_FLOOR, 0.5);
            iterations += self.pcg32(inner_tol, max_iter - iterations)?;
            self.gpu.launch(ADD_F32_TO_F64, blocks(n), THREADS, &[nu, Arg::F64(s), self.cx.arg(), self.u.arg()])?;
        }
        self.u.download(u)?;
        Ok(PcgStats { iterations, rel_residual: rel, converged })
    }

    fn levels(&self) -> Vec<Grid> {
        self.levels.iter().map(|l| l.grid).collect()
    }

    fn name(&self) -> String {
        let i = &self.gpu.info;
        format!("{} (CUDA, sm_{}{})", i.name, i.compute_capability.0, i.compute_capability.1)
    }
}
