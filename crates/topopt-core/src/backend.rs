//! Pluggable linear solvers: what solves `K(E)·u = f` for the optimizer.
//!
//! The optimizer only needs, once per iteration, new element moduli ([`LinearSolver::update`])
//! and then one solve per load case ([`LinearSolver::solve`]). [`CpuSolver`] is the
//! multigrid-preconditioned CG of [`crate::mg`] and [`crate::pcg`]; other crates add other
//! backends (the `topopt-cuda` crate runs the same method on an NVIDIA GPU) through a
//! [`SolverFactory`] given to [`crate::Optimizer::with_solver`].

use crate::fem::{self, Elements};
use crate::grid::Grid;
use crate::ke::Ke;
use crate::mg::{Hierarchy, MgSettings};
use crate::pcg::{PcgStats, pcg};
use std::sync::Arc;

/// Solves the stiffness system of one grid with fixed supports.
pub trait LinearSolver: Send {
    /// Prepares for new element moduli: `emod[e]` scales `ke` (the element matrix of the solid,
    /// with the material's modulus) for element `e`.
    fn update(&mut self, ke: &Ke, emod: &[f64]) -> Result<(), String>;

    /// Solves `K·u = f` starting from `u` (a warm start) until `‖f − K·u‖ ≤ tol·‖f‖` or
    /// `max_iter` iterations. `ke` and `emod` are those of the last [`LinearSolver::update`].
    fn solve(&mut self, ke: &Ke, emod: &[f64], f: &[f64], u: &mut [f64], tol: f64, max_iter: usize) -> Result<PcgStats, String>;

    /// Grids of the multigrid levels, finest first.
    fn levels(&self) -> Vec<Grid>;

    /// What runs the solves, for reports ("CPU, 6 threads", a GPU's name).
    fn name(&self) -> String;

    /// Something the user should know about this solver (for example why a GPU was not used).
    fn note(&self) -> Option<String> {
        None
    }
}

/// Makes a solver for a grid, its fixed dofs and the multigrid settings.
pub type SolverFactory = Arc<dyn Fn(&Grid, &[bool], &MgSettings) -> Result<Box<dyn LinearSolver>, String> + Send + Sync>;

/// The factory of [`CpuSolver`].
pub fn cpu_factory() -> SolverFactory {
    Arc::new(|grid: &Grid, fixed: &[bool], mg: &MgSettings| Ok(Box::new(CpuSolver::new(*grid, fixed.to_vec(), *mg)) as Box<dyn LinearSolver>))
}

/// Multigrid-preconditioned CG on the CPU (all cores with the `parallel` feature).
pub struct CpuSolver {
    grid: Grid,
    fixed: Vec<bool>,
    hierarchy: Hierarchy,
}

impl CpuSolver {
    pub fn new(grid: Grid, fixed: Vec<bool>, mg: MgSettings) -> Self {
        let hierarchy = Hierarchy::new(grid, fixed.clone(), mg);
        CpuSolver { grid, fixed, hierarchy }
    }
}

impl LinearSolver for CpuSolver {
    fn update(&mut self, ke: &Ke, emod: &[f64]) -> Result<(), String> {
        self.hierarchy.update(ke, emod)
    }

    fn solve(&mut self, ke: &Ke, emod: &[f64], f: &[f64], u: &mut [f64], tol: f64, max_iter: usize) -> Result<PcgStats, String> {
        let (grid, fixed, h) = (&self.grid, &self.fixed, &self.hierarchy);
        let elems = Elements::Scaled { ke, emod };
        Ok(pcg(|x, y| fem::apply(grid, &elems, fixed, x, y), |r, z| h.vcycle(ke, emod, r, z), f, u, tol, max_iter))
    }

    fn levels(&self) -> Vec<Grid> {
        self.hierarchy.level_grids()
    }

    fn name(&self) -> String {
        let t = crate::par::threads();
        format!("CPU, {t} thread{}", if t == 1 { "" } else { "s" })
    }
}
