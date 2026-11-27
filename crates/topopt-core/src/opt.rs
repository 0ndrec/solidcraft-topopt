//! The optimization loop: minimum compliance under a volume constraint, SIMP interpolation,
//! density filter, Heaviside projection with β-continuation, optimality-criteria update.
//!
//! Chain per iteration: design `x` → filtered `x̃` → projected `x̄` (passive elements forced to
//! 0 or 1) → modulus `E = E_min + x̄ᵖ·(1 − E_min)` → FE solves → compliance and sensitivities →
//! back through projection and filter → OC update with a bisection on the volume multiplier.

use crate::backend::{CpuSolver, LinearSolver, SolverFactory};
use crate::fem;
use crate::filter::{DensityFilter, project, project_derivative};
use crate::ke::{Ke, hex8_ke};
use crate::mg::MgSettings;
use crate::pcg::{PcgStats, dot};
use crate::problem::{ElementKind, Problem};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Settings {
    /// Target material fraction of the design elements.
    pub volfrac: f64,
    /// SIMP penalty.
    pub penal: f64,
    /// Void stiffness relative to solid.
    pub emin: f64,
    /// Filter radius in element lengths.
    pub rmin: f64,
    /// Heaviside projection on (off: physical density = filtered density, no continuation).
    pub projection: bool,
    /// Projection threshold.
    pub eta: f64,
    pub beta_start: f64,
    pub beta_max: f64,
    /// Double β after this many iterations at one β, or earlier when the design settles.
    pub beta_every: usize,
    /// OC move limit.
    pub move_limit: f64,
    pub max_iter: usize,
    /// Stop when the largest design change is below this at the final β.
    pub tol: f64,
    /// Or when, at the final β, compliance changed by less than this (relative) in each of the
    /// last five iterations. OC with a sharp projection often oscillates in a few elements, so
    /// the design change alone may never fall below `tol`.
    pub compliance_tol: f64,
    pub pcg_tol: f64,
    pub pcg_max_iter: usize,
    pub mg: MgSettings,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            volfrac: 0.3,
            penal: 3.0,
            emin: 1e-9,
            rmin: 1.5,
            projection: true,
            eta: 0.5,
            beta_start: 1.0,
            beta_max: 16.0,
            beta_every: 40,
            move_limit: 0.2,
            max_iter: 300,
            tol: 0.01,
            compliance_tol: 1e-3,
            pcg_tol: 1e-6,
            pcg_max_iter: 200,
            mg: MgSettings::default(),
        }
    }
}

/// What one iteration reports.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Iteration {
    pub iter: usize,
    /// Weighted compliance, N·mm.
    pub compliance: f64,
    /// Material fraction of the design elements (projected densities).
    pub volume: f64,
    /// Largest change of a design variable.
    pub change: f64,
    pub beta: f64,
    /// Measure of non-discreteness: mean of 4·x̄·(1−x̄) over design elements (0 = crisp).
    pub grayness: f64,
    /// CG iterations of this iteration's solves, summed over load cases.
    pub cg_iterations: usize,
    pub cg_worst_residual: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Outcome {
    Converged,
    MaxIterations,
    Cancelled,
}

/// Objective and gradients at the current design.
pub struct Evaluation {
    pub compliance: f64,
    /// ∂c/∂x for every design variable (0 for passive elements).
    pub dc: Vec<f64>,
    /// ∂V/∂x of the summed projected density of design elements.
    pub dv: Vec<f64>,
    pub stats: Vec<PcgStats>,
}

/// One load case of [`Optimizer::check`].
#[derive(Clone, Debug)]
pub struct CaseCheck {
    /// F·u, N·mm.
    pub compliance: f64,
    /// The largest nodal displacement, mm.
    pub max_displacement: f64,
    /// Von Mises stress per element, MPa.
    pub von_mises: Vec<f64>,
}

pub struct Optimizer {
    pub problem: Problem,
    pub settings: Settings,
    ke: Ke,
    forces: Vec<Vec<f64>>,
    weights: Vec<f64>,
    filter: DensityFilter,
    solver: Box<dyn LinearSolver>,
    design: Vec<bool>,
    n_design: usize,
    /// Design variables.
    pub x: Vec<f64>,
    /// Filtered densities.
    pub x_filtered: Vec<f64>,
    /// Physical (projected, passive-forced) densities.
    pub x_phys: Vec<f64>,
    /// Displacements per load case (warm starts).
    pub u: Vec<Vec<f64>>,
    pub beta: f64,
    pub history: Vec<Iteration>,
    iters_at_beta: usize,
}

impl Optimizer {
    pub fn new(problem: Problem, settings: Settings) -> Result<Self, String> {
        Self::with_solver(problem, settings, None)
    }

    /// Like [`Optimizer::new`], with the linear systems solved by what `factory` makes (the
    /// CPU multigrid-CG when `None`).
    pub fn with_solver(problem: Problem, settings: Settings, factory: Option<&SolverFactory>) -> Result<Self, String> {
        problem.validate()?;
        if !(settings.volfrac > 0.0 && settings.volfrac < 1.0) {
            return Err("volume fraction must be between 0 and 1".into());
        }
        let grid = problem.grid;
        let mut ke = hex8_ke(problem.material.nu, grid.h);
        ke.iter_mut().for_each(|v| *v *= problem.material.e);
        let fixed = problem.fixed_dofs();
        let forces = problem.force_vectors(&fixed);
        let weights = problem.load_cases.iter().map(|c| c.weight).collect();
        let filter = DensityFilter::new(grid, settings.rmin);
        let solver: Box<dyn LinearSolver> = match factory {
            Some(make) => make(&grid, &fixed, &settings.mg)?,
            None => Box::new(CpuSolver::new(grid, fixed, settings.mg)),
        };
        let design: Vec<bool> = (0..grid.nel()).map(|e| problem.kind(e) == ElementKind::Design).collect();
        let n_design = design.iter().filter(|&&d| d).count();
        let x: Vec<f64> = (0..grid.nel())
            .map(|e| match problem.kind(e) {
                ElementKind::Design => settings.volfrac,
                ElementKind::Solid => 1.0,
                ElementKind::Void => 0.0,
            })
            .collect();
        let nel = grid.nel();
        let ncases = forces.len();
        let mut opt = Optimizer {
            ke,
            forces,
            weights,
            filter,
            solver,
            design,
            n_design,
            x,
            x_filtered: vec![0.0; nel],
            x_phys: vec![0.0; nel],
            u: vec![vec![0.0; grid.ndof()]; ncases],
            beta: settings.beta_start,
            history: Vec::new(),
            iters_at_beta: 0,
            problem,
            settings,
        };
        opt.update_physical();
        Ok(opt)
    }

    /// Levels of the multigrid hierarchy (finest first).
    pub fn mg_levels(&self) -> Vec<crate::Grid> {
        self.solver.levels()
    }

    /// What solves the linear systems ([`LinearSolver::name`]).
    pub fn solver_name(&self) -> String {
        self.solver.name()
    }

    /// What the user should know about the linear solver ([`LinearSolver::note`]).
    pub fn solver_note(&self) -> Option<String> {
        self.solver.note()
    }

    fn physical_of(&self, x: &[f64], filtered: &mut [f64], phys: &mut [f64]) {
        self.filter.apply(x, filtered);
        for e in 0..x.len() {
            phys[e] = match self.problem.kind(e) {
                ElementKind::Design => self.proj(filtered[e]),
                ElementKind::Solid => 1.0,
                ElementKind::Void => 0.0,
            };
        }
    }

    #[inline]
    fn proj(&self, x: f64) -> f64 {
        if self.settings.projection { project(x, self.beta, self.settings.eta) } else { x }
    }

    #[inline]
    fn proj_derivative(&self, x: f64) -> f64 {
        if self.settings.projection { project_derivative(x, self.beta, self.settings.eta) } else { 1.0 }
    }

    /// Replaces the design variables (passive elements keep their forced values).
    pub fn set_design(&mut self, x: &[f64]) {
        for e in 0..self.x.len() {
            if self.design[e] {
                self.x[e] = x[e].clamp(0.0, 1.0);
            }
        }
        self.update_physical();
    }

    fn update_physical(&mut self) {
        let (mut f, mut p) = (vec![0.0; self.x.len()], vec![0.0; self.x.len()]);
        self.physical_of(&self.x, &mut f, &mut p);
        self.x_filtered = f;
        self.x_phys = p;
    }

    fn design_volume(&self, phys: &[f64]) -> f64 {
        phys.iter().zip(&self.design).filter(|(_, d)| **d).map(|(v, _)| v).sum::<f64>() / self.n_design as f64
    }

    /// Element moduli relative to the solid's.
    pub fn moduli(&self, phys: &[f64]) -> Vec<f64> {
        let (p, emin) = (self.settings.penal, self.settings.emin);
        phys.iter().map(|&x| emin + x.powf(p) * (1.0 - emin)).collect()
    }

    /// Solves every load case for the current physical densities; returns per-case stats.
    pub fn solve(&mut self) -> Result<Vec<PcgStats>, String> {
        let emod = self.moduli(&self.x_phys);
        self.solver.update(&self.ke, &emod)?;
        let mut stats = Vec::new();
        for c in 0..self.forces.len() {
            let (tol, max_iter) = (self.settings.pcg_tol, self.settings.pcg_max_iter);
            stats.push(self.solver.solve(&self.ke, &emod, &self.forces[c], &mut self.u[c], tol, max_iter)?);
        }
        Ok(stats)
    }

    /// Analysis of the current physical densities `x_phys` (set them directly to analyse any
    /// design): the weighted compliance and each element's strain energy, `E_e · u_eᵀ k u_e`
    /// summed over load cases (relative units; compare elements with each other).
    pub fn evaluate_energies(&mut self) -> Result<(f64, Vec<f64>), String> {
        let nel = self.problem.grid.nel();
        self.solve()?;
        let emod = self.moduli(&self.x_phys);
        let grid = self.problem.grid;
        let mut compliance = 0.0;
        let mut energy = vec![0.0; nel];
        let mut ce = vec![0.0; nel];
        for c in 0..self.forces.len() {
            compliance += self.weights[c] * dot(&self.forces[c], &self.u[c]);
            fem::element_energies(&grid, &self.ke, &self.u[c], &mut ce);
            for e in 0..nel {
                energy[e] += self.weights[c] * emod[e] * ce[e];
            }
        }
        Ok((compliance, energy))
    }

    /// Analyses a given design (physical densities, one per element; usually the result made
    /// 0 or 1) and reports, per load case, its compliance, largest displacement and von Mises
    /// stresses (see [`crate::stress`]). The design stays in `x_phys`.
    pub fn check(&mut self, phys: &[f64]) -> Result<Vec<CaseCheck>, String> {
        if phys.len() != self.problem.grid.nel() {
            return Err("check: one density per element".into());
        }
        self.x_phys = phys.to_vec();
        self.solve()?;
        let emod = self.moduli(&self.x_phys);
        let m = self.problem.material;
        Ok((0..self.forces.len())
            .map(|c| CaseCheck {
                compliance: dot(&self.forces[c], &self.u[c]),
                max_displacement: crate::stress::max_displacement(&self.u[c]),
                von_mises: crate::stress::von_mises(&self.problem.grid, m.e, m.nu, &emod, &self.u[c]),
            })
            .collect())
    }

    /// Compliance and its gradients with respect to the design variables at the current design.
    pub fn evaluate(&mut self) -> Result<Evaluation, String> {
        let nel = self.problem.grid.nel();
        let stats = self.solve()?;
        let grid = self.problem.grid;
        let p = self.settings.penal;
        let emin = self.settings.emin;

        // Compliance and its sensitivity to the physical densities.
        let mut compliance = 0.0;
        let mut energy = vec![0.0; nel];
        let mut ce = vec![0.0; nel];
        for c in 0..self.forces.len() {
            compliance += self.weights[c] * dot(&self.forces[c], &self.u[c]);
            fem::element_energies(&grid, &self.ke, &self.u[c], &mut ce);
            for e in 0..nel {
                energy[e] += self.weights[c] * ce[e];
            }
        }
        let mut dc = vec![0.0; nel];
        let mut dv = vec![0.0; nel];
        for e in 0..nel {
            if self.design[e] {
                let xp = self.x_phys[e];
                let dproj = self.proj_derivative(self.x_filtered[e]);
                dc[e] = -p * xp.powf(p - 1.0) * (1.0 - emin) * energy[e] * dproj;
                dv[e] = dproj;
            }
        }
        let (mut dc_x, mut dv_x) = (vec![0.0; nel], vec![0.0; nel]);
        self.filter.apply_transpose(&dc, &mut dc_x);
        self.filter.apply_transpose(&dv, &mut dv_x);
        for e in 0..nel {
            if !self.design[e] {
                dc_x[e] = 0.0;
                dv_x[e] = 0.0;
            }
        }
        Ok(Evaluation { compliance, dc: dc_x, dv: dv_x, stats })
    }

    /// One optimization iteration. Returns its report.
    pub fn step(&mut self) -> Result<Iteration, String> {
        let nel = self.problem.grid.nel();
        let Evaluation { compliance, dc: dc_x, dv: dv_x, stats } = self.evaluate()?;

        // OC update with bisection on the multiplier so the projected volume meets the target.
        let target = self.settings.volfrac;
        let mv = self.settings.move_limit;
        let candidate = |lambda: f64, out: &mut [f64]| {
            for e in 0..nel {
                if !self.design[e] {
                    out[e] = self.x[e];
                    continue;
                }
                let b = (-dc_x[e]).max(0.0) / (lambda * dv_x[e].max(1e-12));
                let xe = self.x[e];
                out[e] = (xe * b.sqrt()).clamp((xe - mv).max(0.0), (xe + mv).min(1.0));
            }
        };
        let mut xnew = vec![0.0; nel];
        let (mut f, mut ph) = (vec![0.0; nel], vec![0.0; nel]);
        let (mut lo, mut hi) = (1e-30f64, 1e30f64);
        for _ in 0..200 {
            let mid = (lo * hi).sqrt();
            candidate(mid, &mut xnew);
            self.physical_of(&xnew, &mut f, &mut ph);
            if self.design_volume(&ph) > target {
                lo = mid;
            } else {
                hi = mid;
            }
            if hi / lo < 1.0 + 1e-6 {
                break;
            }
        }
        candidate(hi, &mut xnew);
        let change = (0..nel).filter(|&e| self.design[e]).map(|e| (xnew[e] - self.x[e]).abs()).fold(0.0, f64::max);
        self.x = xnew;
        self.update_physical();

        let volume = self.design_volume(&self.x_phys);
        let grayness = (0..nel).filter(|&e| self.design[e]).map(|e| 4.0 * self.x_phys[e] * (1.0 - self.x_phys[e])).sum::<f64>()
            / self.n_design as f64;
        let it = Iteration {
            iter: self.history.len() + 1,
            compliance,
            volume,
            change,
            beta: self.beta,
            grayness,
            cg_iterations: stats.iter().map(|s| s.iterations).sum(),
            cg_worst_residual: stats.iter().map(|s| s.rel_residual).fold(0.0, f64::max),
        };
        self.history.push(it);

        // β-continuation.
        self.iters_at_beta += 1;
        if self.settings.projection && self.beta < self.settings.beta_max && (self.iters_at_beta >= self.settings.beta_every || change < self.settings.tol) {
            self.beta = (self.beta * 2.0).min(self.settings.beta_max);
            self.iters_at_beta = 0;
            self.update_physical();
        }
        Ok(it)
    }

    /// Compliance moved less than `compliance_tol` in each of the last five iterations, all at the
    /// final β.
    fn compliance_settled(&self) -> bool {
        let h = &self.history;
        if h.len() < 6 {
            return false;
        }
        let tail = &h[h.len() - 6..];
        tail.iter().all(|i| i.beta == tail[5].beta)
            && tail.windows(2).all(|w| ((w[1].compliance - w[0].compliance) / w[1].compliance).abs() < self.settings.compliance_tol)
    }

    /// Runs until converged, out of iterations, or `progress` returns false.
    pub fn run(&mut self, mut progress: impl FnMut(&Iteration) -> bool) -> Result<Outcome, String> {
        while self.history.len() < self.settings.max_iter {
            let it = self.step()?;
            if !progress(&it) {
                return Ok(Outcome::Cancelled);
            }
            let final_beta = !self.settings.projection || it.beta >= self.settings.beta_max;
            if final_beta && (it.change < self.settings.tol || self.compliance_settled()) {
                return Ok(Outcome::Converged);
            }
        }
        Ok(Outcome::MaxIterations)
    }
}
