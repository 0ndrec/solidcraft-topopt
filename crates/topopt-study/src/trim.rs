//! Trimming dead weight: the "stalactites and stalagmites" an optimization result often keeps —
//! thin pieces sticking out of the structure that add mass but almost no stiffness.
//!
//! The criterion is the one that matters for the part: **removing the piece does not make the
//! part noticeably less stiff**. It is tested directly, by analysis, on the part as it will be
//! made (the design thresholded at the mesh's iso value):
//! 1. candidates come from geometry: the thin parts of the solid, those a morphological opening
//!    with a 3×3×3 cube removes (struts and spikes up to about two elements thick), grouped into
//!    face-connected pieces of at least `min_cells` cells. Pads under supports and loads and
//!    `keep_in` regions are never candidates;
//! 2. each piece, smallest first, is removed on trial and the compliance recomputed (one FE
//!    solve, warm-started). The piece goes if it raises the compliance by less than `piece_loss`,
//!    and as long as all removals together stay below `max_loss`.
//!
//! A first version selected elements by low strain energy. On real results it found almost
//! nothing: the dangling-looking strands do reach other members through thin necks and carry some
//! load. What matters is how much, which the trial removal measures.

use serde::{Deserialize, Serialize};
use topopt_core::{ElementKind, Grid, Optimizer, Problem};

/// Settings of [`trim_dead_ends`].
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct TrimSettings {
    /// A piece may raise the compliance by at most this (relative).
    pub piece_loss: f64,
    /// All removals together may raise it by at most this (relative).
    pub max_loss: f64,
    /// Smallest piece considered, in cells.
    pub min_cells: usize,
}

impl Default for TrimSettings {
    fn default() -> Self {
        TrimSettings { piece_loss: 0.001, max_loss: 0.005, min_cells: 3 }
    }
}

/// What [`trim_dead_ends`] did.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct TrimReport {
    /// Thin pieces found, and those removed.
    pub candidates: usize,
    pub removed_pieces: usize,
    pub removed_elements: usize,
    /// Removed share of the solid elements.
    pub removed_fraction: f64,
    pub compliance_before: f64,
    pub compliance_after: f64,
}

fn neighbours6(grid: &Grid, e: usize) -> impl Iterator<Item = usize> + '_ {
    let [i, j, k] = grid.elem_ijk(e);
    let (i, j, k) = (i as i64, j as i64, k as i64);
    [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)].into_iter().filter_map(move |(di, dj, dk)| {
        let (a, b, c) = (i + di, j + dj, k + dk);
        (a >= 0 && b >= 0 && c >= 0 && a < grid.nx as i64 && b < grid.ny as i64 && c < grid.nz as i64).then(|| grid.elem(a as usize, b as usize, c as usize))
    })
}

/// Erosion (`grow` false) or dilation (`grow` true) with a 3×3×3 cube.
fn morph(grid: &Grid, cells: &[bool], grow: bool) -> Vec<bool> {
    (0..grid.nel())
        .map(|e| {
            let [i, j, k] = grid.elem_ijk(e);
            let mut any = false;
            let mut all = true;
            for dk in -1i64..=1 {
                for dj in -1i64..=1 {
                    for di in -1i64..=1 {
                        let (a, b, c) = (i as i64 + di, j as i64 + dj, k as i64 + dk);
                        let v = a >= 0
                            && b >= 0
                            && c >= 0
                            && a < grid.nx as i64
                            && b < grid.ny as i64
                            && c < grid.nz as i64
                            && cells[grid.elem(a as usize, b as usize, c as usize)];
                        any |= v;
                        all &= v;
                    }
                }
            }
            if grow { any } else { all }
        })
        .collect()
}

/// Empties in `density` (values and filtered field) the cells that [`trim_dead_ends`] emptied in
/// `trimmed` (its output for `density.values`).
pub fn apply(density: &mut crate::Density, trimmed: &[f64]) {
    for (e, &t) in trimmed.iter().enumerate() {
        if t == 0.0 && density.values[e] != 0.0 {
            density.values[e] = 0.0;
            if let Some(f) = density.filtered.get_mut(e) {
                *f = 0.0;
            }
        }
    }
}

/// Removes thin pieces that add (almost) no stiffness from physical densities `x` (one per
/// element, 0…1) of `problem`, as meshed at `iso`. Returns the densities with those pieces emptied
/// (everything else unchanged, grey values included) and a report. `settings` are the optimizer's
/// (material interpolation, solver).
pub fn trim_dead_ends(problem: &Problem, settings: &topopt_core::Settings, x: &[f64], iso: f64, trim: &TrimSettings) -> Result<(Vec<f64>, TrimReport), String> {
    trim_dead_ends_with(problem, settings, x, iso, trim, None)
}

/// [`trim_dead_ends`] with the analyses solved by what `solver` makes (the CPU when `None`).
pub fn trim_dead_ends_with(
    problem: &Problem,
    settings: &topopt_core::Settings,
    x: &[f64],
    iso: f64,
    trim: &TrimSettings,
    solver: Option<&topopt_core::SolverFactory>,
) -> Result<(Vec<f64>, TrimReport), String> {
    let grid = problem.grid;
    let nel = grid.nel();
    let mut opt = Optimizer::with_solver(problem.clone(), *settings, solver)?;
    let solid: Vec<bool> = x.iter().map(|&v| v >= iso).collect();
    let mut cur: Vec<f64> = solid.iter().map(|&s| if s { 1.0 } else { 0.0 }).collect();
    let analyse = |opt: &mut Optimizer, d: &[f64]| -> Result<f64, String> {
        opt.x_phys = d.to_vec();
        Ok(opt.evaluate_energies()?.0)
    };
    let c0 = analyse(&mut opt, &cur)?;
    let mut report = TrimReport { compliance_before: c0, compliance_after: c0, ..Default::default() };

    // Thin parts: solid minus its opening, outside pads and keep-in regions.
    let opened = morph(&grid, &morph(&grid, &solid, false), true);
    let thin: Vec<bool> = (0..nel).map(|e| solid[e] && !opened[e] && problem.kind(e) == ElementKind::Design).collect();
    let mut seen = vec![false; nel];
    let mut pieces: Vec<Vec<usize>> = Vec::new();
    for start in 0..nel {
        if !thin[start] || seen[start] {
            continue;
        }
        seen[start] = true;
        let (mut stack, mut members) = (vec![start], Vec::new());
        while let Some(e) = stack.pop() {
            members.push(e);
            for f in neighbours6(&grid, e) {
                if thin[f] && !seen[f] {
                    seen[f] = true;
                    stack.push(f);
                }
            }
        }
        if members.len() >= trim.min_cells {
            pieces.push(members);
        }
    }
    pieces.sort_by_key(|p| p.len());
    report.candidates = pieces.len();

    let mut c_cur = c0;
    for p in &pieces {
        let mut trial = cur.clone();
        for &e in p {
            trial[e] = 0.0;
        }
        // A removal that leaves part of the structure floating makes the system singular: that
        // piece holds something, so it stays.
        let Ok(c1) = analyse(&mut opt, &trial) else { continue };
        if c1 - c_cur <= trim.piece_loss * c0 && c1 <= c0 * (1.0 + trim.max_loss) {
            cur = trial;
            c_cur = c1;
            report.removed_pieces += 1;
            report.removed_elements += p.len();
        }
    }
    report.compliance_after = c_cur;
    report.removed_fraction = report.removed_elements as f64 / solid.iter().filter(|&&s| s).count().max(1) as f64;
    let out = x.iter().zip(&cur).map(|(&v, &c)| if v >= iso && c == 0.0 { 0.0 } else { v }).collect();
    Ok((out, report))
}
