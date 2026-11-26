//! Geometric multigrid V-cycle used as the CG preconditioner.
//!
//! Level 0 is the fine grid, applied matrix-free (`emod[e] · ke`). Each coarser level halves the
//! grid and stores one 24×24 matrix per element: the exact Galerkin product `Pᵀ·K·P`, formed
//! element by element. Trilinear prolongation from a coarse element only reaches the 27 fine nodes
//! inside it, and agrees on shared faces, so summing `P_cᵀ·K_child·P_c` over the eight children of
//! each coarse element gives the assembled `Pᵀ·K·P` exactly. The coarsest level is factorized
//! (sparse Cholesky); the pattern is fixed, so its symbolic analysis is done once.
//!
//! Supports: a coarse dof is fixed when any fine node near the coarse node is fixed in that
//! direction. Restriction and prolongation mask fixed dofs on both sides (`R = Pᵀ` of the masked
//! `P`), so the V-cycle stays symmetric positive definite.

use crate::fem::{self, Elements, Ke32};
use crate::grid::{CORNERS, Grid};
use crate::ke::Ke;
use faer::Side;
use faer::linalg::solvers::SolveCore;
use faer::sparse::linalg::solvers::{Llt, SymbolicLlt};
use faer::sparse::{Argsort, Pair, SparseColMat, SymbolicSparseColMat};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct MgSettings {
    /// Stop coarsening once a level has at most this many dofs (it is then factorized).
    pub max_coarse_dofs: usize,
    /// Damped-Jacobi sweeps before and after the coarse correction.
    pub smooth_steps: usize,
    /// Jacobi damping as a fraction of the stable limit: each level uses
    /// `omega · 2 / λ_max(D⁻¹A)`, with `λ_max` estimated by power iteration.
    pub omega: f64,
    /// Power-iteration steps for `λ_max` after the first update (which uses 12), warm-started.
    pub power_steps: usize,
    /// A fine grid with at most this many dofs is factorized directly, without coarsening. Thin
    /// domains coarsen badly (members one element thick vanish on a grid twice as coarse), and a
    /// direct solve of a small problem is cheap anyway.
    pub direct_dofs: usize,
}

impl Default for MgSettings {
    fn default() -> Self {
        // Factorize only small levels, natively as in the sandbox. Multigrid-CG beats faer's
        // sparse Cholesky well before 10 000 dofs, also on thin domains where its CG count grows
        // (native, 6 threads, `examples/profile`): 36×12×12 (18 759 dofs) 52 ms against 370 ms
        // per iteration; a 64×32×4 plate over 160 iterations 47 s against 84 s, the same
        // compliance; 20×10×10 with a 10×5×5 coarse level 31 ms against 82 ms. In wasm the
        // factorization is 8× slower still (a 4000-element study: 65 s direct, 7.7 s without).
        MgSettings { max_coarse_dofs: 4_000, smooth_steps: 2, omega: 0.6, power_steps: 3, direct_dofs: 4_000 }
    }
}

struct Level {
    grid: Grid,
    fixed: Vec<bool>,
    /// Element matrices (levels ≥ 1).
    mats: Vec<Ke32>,
    inv_diag: Vec<f64>,
    /// Jacobi step for this level.
    omega: f64,
    /// Power-iteration vector, kept as the warm start for the next update.
    eig: Vec<f64>,
}

struct Coarse {
    /// Lower-triangle (row ≥ col) entries in pattern order: (row, col).
    pattern: Vec<(usize, usize)>,
    symbolic: SymbolicSparseColMat<usize>,
    argsort: Argsort<usize>,
    llt_symbolic: SymbolicLlt<usize>,
    llt: Option<Llt<usize, f64>>,
    /// For each element and each lower-triangle (r, c) of its 24×24 matrix with both dofs free,
    /// the slot in `pattern` it adds to.
    scatter: Vec<u32>,
}

pub struct Hierarchy {
    pub settings: MgSettings,
    fine_grid: Grid,
    fine_fixed: Vec<bool>,
    fine_inv_diag: Vec<f64>,
    fine_omega: f64,
    fine_eig: Vec<f64>,
    levels: Vec<Level>,
    /// `P_cᵀ · ke · P_c` for the eight child positions of a level-1 element (unscaled).
    child_galerkin: Option<[Ke; 8]>,
    coarse: Option<Coarse>,
}

/// Scalar trilinear weights from a coarse element's 8 corners to child `c`'s 8 corners:
/// `p[a][b]` is the weight of coarse corner `b` at fine corner `a` of the child.
fn child_prolongation(c: usize) -> [[f64; 8]; 8] {
    let off = CORNERS[c];
    let mut p = [[0.0; 8]; 8];
    for a in 0..8 {
        // Position of the fine corner in the coarse element's 3×3×3 node lattice, as 0, ½ or 1.
        let t: Vec<f64> = (0..3).map(|d| (CORNERS[a][d] + off[d]) as f64 / 2.0).collect();
        for b in 0..8 {
            p[a][b] = (0..3).map(|d| if CORNERS[b][d] == 1 { t[d] } else { 1.0 - t[d] }).product();
        }
    }
    p
}

/// `out += Pᵀ · m · P` for the scalar 8×8 `p`, expanded per direction (identity 3×3 blocks).
fn galerkin_add(m: &Ke, p: &[[f64; 8]; 8], out: &mut Ke) {
    // t[(a,i),(B,j)] = Σ_b m[(a,i),(b,j)] · p[b][B]
    let mut t = [0.0f64; 576];
    for r in 0..24 {
        for bc in 0..8 {
            for j in 0..3 {
                let mut acc = 0.0;
                for b in 0..8 {
                    let w = p[b][bc];
                    if w != 0.0 {
                        acc += m[r * 24 + 3 * b + j] * w;
                    }
                }
                t[r * 24 + 3 * bc + j] = acc;
            }
        }
    }
    // out[(A,i),(B,j)] += Σ_a p[a][A] · t[(a,i),(B,j)]
    for ac in 0..8 {
        for a in 0..8 {
            let w = p[a][ac];
            if w == 0.0 {
                continue;
            }
            for i in 0..3 {
                let src = &t[(3 * a + i) * 24..(3 * a + i) * 24 + 24];
                let dst = &mut out[(3 * ac + i) * 24..(3 * ac + i) * 24 + 24];
                for c in 0..24 {
                    dst[c] += w * src[c];
                }
            }
        }
    }
}

/// For a fine node coordinate, the coarse coordinates it interpolates from, with weights.
#[inline]
fn axis_weights(i: usize) -> ([usize; 2], [f64; 2], usize) {
    if i % 2 == 0 { ([i / 2, 0], [1.0, 0.0], 1) } else { ([(i - 1) / 2, (i + 1) / 2], [0.5, 0.5], 2) }
}

/// `fine += P · coarse` (P masked by both fixed sets).
fn prolong_add(cg: &Grid, cfix: &[bool], fg: &Grid, ffix: &[bool], coarse: &[f64], fine: &mut [f64]) {
    let plane = 3 * (fg.nx + 1) * (fg.ny + 1);
    crate::par::chunks(fine, plane, |k, slab| {
        let base = k * plane;
        let (kz, wz, nz) = axis_weights(k);
        for j in 0..=fg.ny {
            let (jy, wy, ny) = axis_weights(j);
            for i in 0..=fg.nx {
                let (ix, wx, nx) = axis_weights(i);
                let fnode = fg.node(i, j, k);
                let mut v = [0.0f64; 3];
                for c in 0..nz {
                    for b in 0..ny {
                        for a in 0..nx {
                            let w = wx[a] * wy[b] * wz[c];
                            let cn = cg.node(ix[a], jy[b], kz[c]);
                            for d in 0..3 {
                                if !cfix[3 * cn + d] {
                                    v[d] += w * coarse[3 * cn + d];
                                }
                            }
                        }
                    }
                }
                for d in 0..3 {
                    if !ffix[3 * fnode + d] {
                        slab[3 * fnode + d - base] += v[d];
                    }
                }
            }
        }
    });
}

/// `coarse = Pᵀ · fine` (masked), the transpose of [`prolong_add`], written as a gather so each
/// coarse node is computed by one thread: coarse node `I` collects fine node `2I` with weight 1
/// and `2I ± 1` with weight ½, per axis.
fn restrict(cg: &Grid, cfix: &[bool], fg: &Grid, ffix: &[bool], fine: &[f64], coarse: &mut [f64]) {
    let plane = 3 * (cg.nx + 1) * (cg.ny + 1);
    crate::par::chunks(coarse, plane, |kc, slab| {
        let base = kc * plane;
        let span = |c: usize, n: usize| -> ([usize; 3], [f64; 3], usize) {
            let mut idx = [0usize; 3];
            let mut w = [0.0f64; 3];
            let mut m = 0;
            for (f, wt) in [(2 * c as i64 - 1, 0.5), (2 * c as i64, 1.0), (2 * c as i64 + 1, 0.5)] {
                if f >= 0 && f <= n as i64 {
                    idx[m] = f as usize;
                    w[m] = wt;
                    m += 1;
                }
            }
            (idx, w, m)
        };
        let (kf, wz, nz) = span(kc, fg.nz);
        for jc in 0..=cg.ny {
            let (jf, wy, ny) = span(jc, fg.ny);
            for ic in 0..=cg.nx {
                let (if_, wx, nx) = span(ic, fg.nx);
                let cn = cg.node(ic, jc, kc);
                let mut v = [0.0f64; 3];
                for c in 0..nz {
                    for b in 0..ny {
                        for a in 0..nx {
                            let w = wx[a] * wy[b] * wz[c];
                            let fnode = fg.node(if_[a], jf[b], kf[c]);
                            for d in 0..3 {
                                if !ffix[3 * fnode + d] {
                                    v[d] += w * fine[3 * fnode + d];
                                }
                            }
                        }
                    }
                }
                for d in 0..3 {
                    slab[3 * cn + d - base] = if cfix[3 * cn + d] { 0.0 } else { v[d] };
                }
            }
        }
    });
}

/// A coarse dof is fixed when any fine node within one fine element of the coarse node is.
fn coarse_fixed(cg: &Grid, fg: &Grid, ffix: &[bool]) -> Vec<bool> {
    let mut fixed = vec![false; cg.ndof()];
    for k in 0..=cg.nz {
        for j in 0..=cg.ny {
            for i in 0..=cg.nx {
                let cn = cg.node(i, j, k);
                let range = |c: usize, n: usize| (2 * c).saturating_sub(1)..=(2 * c + 1).min(n);
                for fk in range(k, fg.nz) {
                    for fj in range(j, fg.ny) {
                        for fi in range(i, fg.nx) {
                            let fnode = fg.node(fi, fj, fk);
                            for d in 0..3 {
                                if ffix[3 * fnode + d] {
                                    fixed[3 * cn + d] = true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    fixed
}

/// Largest eigenvalue of `D⁻¹·A` by power iteration (slightly overestimated, as it should be for
/// a damping bound). `v` is the previous update's vector: densities change little between
/// optimization iterations, so a few steps from it suffice.
fn lambda_max(grid: &Grid, elems: &Elements, fixed: &[bool], inv_diag: &[f64], v: &mut Vec<f64>, warm_steps: usize) -> f64 {
    let n = grid.ndof();
    let steps = if v.len() == n { warm_steps.max(1) } else { 12 };
    if v.len() != n {
        // A deterministic, non-smooth start vector.
        *v = (0..n).map(|i| if fixed[i] { 0.0 } else { 1.0 + ((i * 7919) % 13) as f64 / 13.0 }).collect();
    }
    let mut av = vec![0.0; n];
    let mut lambda = 1.0;
    for _ in 0..steps {
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm == 0.0 {
            return 1.0;
        }
        v.iter_mut().for_each(|x| *x /= norm);
        fem::apply(grid, elems, fixed, &v, &mut av);
        for i in 0..n {
            av[i] = if fixed[i] { 0.0 } else { av[i] * inv_diag[i] };
        }
        lambda = v.iter().zip(&av).map(|(a, b)| a * b).sum::<f64>();
        std::mem::swap(v, &mut av);
    }
    (lambda * 1.1).max(1.0)
}

fn invert(diag: Vec<f64>) -> Vec<f64> {
    diag.into_iter().map(|d| 1.0 / d).collect()
}

impl Hierarchy {
    /// Builds the grid levels and support masks; densities come with [`Hierarchy::update`].
    pub fn new(fine_grid: Grid, fine_fixed: Vec<bool>, settings: MgSettings) -> Self {
        let mut levels: Vec<Level> = Vec::new();
        let (mut g, mut fix) = (fine_grid, fine_fixed.clone());
        while g.ndof() > settings.max_coarse_dofs && fine_grid.ndof() > settings.direct_dofs {
            let Some(cg) = g.coarsen() else { break };
            let cfix = coarse_fixed(&cg, &g, &fix);
            levels.push(Level { grid: cg, fixed: cfix.clone(), mats: vec![[0.0f32; 576]; cg.nel()], inv_diag: Vec::new(), omega: 0.0, eig: Vec::new() });
            g = cg;
            fix = cfix;
        }
        Hierarchy {
            settings,
            fine_grid,
            fine_fixed,
            fine_inv_diag: Vec::new(),
            fine_omega: 0.0,
            fine_eig: Vec::new(),
            levels,
            child_galerkin: None,
            coarse: None,
        }
    }

    /// Number of levels including the fine one.
    pub fn depth(&self) -> usize {
        self.levels.len() + 1
    }

    /// Grid sizes of every level, finest first.
    pub fn level_grids(&self) -> Vec<Grid> {
        std::iter::once(self.fine_grid).chain(self.levels.iter().map(|l| l.grid)).collect()
    }

    /// Rebuilds coarse operators, smoother diagonals and the coarse factorization for new
    /// element moduli `emod` (scaling `ke`).
    pub fn update(&mut self, ke: &Ke, emod: &[f64]) -> Result<(), String> {
        let fine = Elements::Scaled { ke, emod };
        self.fine_inv_diag = invert(fem::diagonal(&self.fine_grid, &fine, &self.fine_fixed));
        self.fine_omega = self.settings.omega * 2.0 / lambda_max(&self.fine_grid, &fine, &self.fine_fixed, &self.fine_inv_diag, &mut self.fine_eig, self.settings.power_steps);
        self.galerkin(ke, emod);
        for level in &mut self.levels {
            level.inv_diag = invert(fem::diagonal(&level.grid, &Elements::Own(&level.mats), &level.fixed));
            level.omega = self.settings.omega * 2.0 / lambda_max(&level.grid, &Elements::Own(&level.mats), &level.fixed, &level.inv_diag, &mut level.eig, self.settings.power_steps);
        }
        self.factorize_coarsest(ke, emod)
    }

    /// Rebuilds only the coarse element matrices and the coarsest factorization, not the
    /// smoothers' diagonals and damping: for a backend that runs the V-cycle itself and computes
    /// those on its own (`topopt-cuda`).
    pub fn update_coarse(&mut self, ke: &Ke, emod: &[f64]) -> Result<(), String> {
        self.galerkin(ke, emod);
        self.factorize_coarsest(ke, emod)
    }

    /// `P_cᵀ·ke·P_c` for the eight child positions `c` of a first-level coarse element (with the
    /// unscaled `ke` of the last update): that element's matrix is `Σ_c emod[child c] · g[c]`.
    pub fn child_galerkin(&self) -> Option<&[Ke; 8]> {
        self.child_galerkin.as_ref()
    }

    /// Number of coarse levels (the levels below the fine grid).
    pub fn coarse_levels(&self) -> usize {
        self.levels.len()
    }

    /// Coarse level `l` (`0` is the level below the fine grid): its grid, fixed dofs and element
    /// matrices (after an update).
    pub fn coarse_level(&self, l: usize) -> (Grid, &[bool], &[Ke32]) {
        let lv = &self.levels[l];
        (lv.grid, &lv.fixed, &lv.mats)
    }

    /// The Galerkin element matrices of every coarse level.
    fn galerkin(&mut self, ke: &Ke, emod: &[f64]) {
        if !self.levels.is_empty() && self.child_galerkin.is_none() {
            let mut cg = [[0.0f64; 576]; 8];
            for (c, m) in cg.iter_mut().enumerate() {
                galerkin_add(ke, &child_prolongation(c), m);
            }
            self.child_galerkin = Some(cg);
        }
        let prolongs: Vec<[[f64; 8]; 8]> = (0..8).map(child_prolongation).collect();
        for l in 0..self.levels.len() {
            let (done, rest) = self.levels.split_at_mut(l);
            let level = &mut rest[0];
            let cgrid = level.grid;
            let fine_grid = self.fine_grid;
            let child_galerkin = &self.child_galerkin;
            let prolongs = &prolongs;
            let done: &[Level] = done;
            crate::par::chunks(&mut level.mats, 1, |e, slot| {
                // Summed in double precision, stored in single.
                let mut m = [0.0f64; 576];
                let [i, j, k] = cgrid.elem_ijk(e);
                for (c, off) in CORNERS.iter().enumerate() {
                    let (ci, cj, ck) = (2 * i + off[0], 2 * j + off[1], 2 * k + off[2]);
                    if l == 0 {
                        let s = emod[fine_grid.elem(ci, cj, ck)];
                        let g = &child_galerkin.as_ref().unwrap()[c];
                        for (v, gv) in m.iter_mut().zip(g.iter()) {
                            *v += s * gv;
                        }
                    } else {
                        let fine = &done[l - 1];
                        let child32 = &fine.mats[fine.grid.elem(ci, cj, ck)];
                        let child: Ke = std::array::from_fn(|q| child32[q] as f64);
                        galerkin_add(&child, &prolongs[c], &mut m);
                    }
                }
                for (dst, src) in slot[0].iter_mut().zip(m.iter()) {
                    *dst = *src as f32;
                }
            });
        }
    }

    fn factorize_coarsest(&mut self, ke: &Ke, emod: &[f64]) -> Result<(), String> {
        let (grid, fixed) = match self.levels.last() {
            Some(l) => (l.grid, l.fixed.clone()),
            None => (self.fine_grid, self.fine_fixed.clone()),
        };
        let n = grid.ndof();
        if self.coarse.is_none() {
            // Pattern: every lower-triangle (r ≥ c) pair of free dofs sharing an element, plus the
            // diagonal of every dof. Built once.
            let mut slot = std::collections::HashMap::<(usize, usize), u32>::new();
            let mut pattern = Vec::new();
            let mut add = |r: usize, c: usize, pattern: &mut Vec<(usize, usize)>| -> u32 {
                *slot.entry((r, c)).or_insert_with(|| {
                    pattern.push((r, c));
                    (pattern.len() - 1) as u32
                })
            };
            for d in 0..n {
                add(d, d, &mut pattern);
            }
            let mut scatter = Vec::new();
            fem::for_each_element(&grid, |_, dofs| {
                for r in 0..24 {
                    for c in 0..=r {
                        let (dr, dc) = (dofs[r], dofs[c]);
                        if fixed[dr] || fixed[dc] {
                            scatter.push(u32::MAX);
                        } else {
                            let (a, b) = if dr >= dc { (dr, dc) } else { (dc, dr) };
                            scatter.push(add(a, b, &mut pattern));
                        }
                    }
                }
            });
            let pairs: Vec<Pair<usize, usize>> = pattern.iter().map(|&(r, c)| Pair { row: r, col: c }).collect();
            let (symbolic, argsort) =
                SymbolicSparseColMat::try_new_from_indices(n, n, &pairs).map_err(|e| format!("coarse pattern: {e:?}"))?;
            let llt_symbolic =
                SymbolicLlt::try_new(symbolic.as_ref(), Side::Lower).map_err(|e| format!("coarse symbolic Cholesky: {e:?}"))?;
            self.coarse = Some(Coarse { pattern, symbolic, argsort, llt_symbolic, llt: None, scatter });
        }
        let coarse = self.coarse.as_mut().unwrap();
        let mut vals = vec![0.0f64; coarse.pattern.len()];
        let mut idx = 0;
        let mats = self.levels.last().map(|l| &l.mats);
        fem::for_each_element(&grid, |e, _| {
            let entry = |q: usize| -> f64 {
                match mats {
                    Some(m) => m[e][q] as f64,
                    None => emod[e] * ke[q],
                }
            };
            for r in 0..24 {
                for c in 0..=r {
                    let slot = coarse.scatter[idx];
                    idx += 1;
                    if slot != u32::MAX {
                        vals[slot as usize] += entry(r * 24 + c);
                    }
                }
            }
        });
        // Fixed dofs: identity. A tiny shift keeps the factorization defined if a region has no
        // support reaching it at this level.
        let max_diag = (0..n).filter(|&d| !fixed[d]).map(|d| vals[d]).fold(0.0f64, f64::max);
        for d in 0..n {
            if fixed[d] || vals[d] <= 0.0 {
                vals[d] = if fixed[d] { 1.0 } else { max_diag.max(1.0) * 1e-12 };
            } else {
                vals[d] += max_diag * 1e-12;
            }
        }
        let mat = SparseColMat::new_from_argsort(coarse.symbolic.clone(), &coarse.argsort, &vals)
            .map_err(|e| format!("coarse matrix: {e:?}"))?;
        let llt = Llt::try_new_with_symbolic(coarse.llt_symbolic.clone(), mat.as_ref(), Side::Lower)
            .map_err(|e| format!("coarse Cholesky: {e:?}"))?;
        coarse.llt = Some(llt);
        Ok(())
    }

    /// Solves the coarsest level's system (the fine one when there are no coarse levels) with
    /// its factorization: `x = A⁻¹·b`.
    pub fn coarse_solve(&self, b: &[f64], x: &mut [f64]) {
        x.copy_from_slice(b);
        let llt = self.coarse.as_ref().and_then(|c| c.llt.as_ref()).expect("update() before solving");
        let mat = faer::MatMut::from_column_major_slice_mut(x, b.len(), 1);
        llt.solve_in_place_with_conj(faer::Conj::No, mat);
    }

    /// Diagnostics: milliseconds of one operator application per level, one restriction and
    /// prolongation per level pair, and one coarse solve.
    #[doc(hidden)]
    pub fn timing_breakdown(&self, ke: &Ke, emod: &[f64]) -> Vec<(String, f64)> {
        let now = || std::time::Instant::now();
        let mut out = Vec::new();
        for l in 0..=self.levels.len() {
            let (grid, elems, fixed, _, _) = self.level_parts(l, ke, emod);
            let x = vec![1.0; grid.ndof()];
            let mut y = vec![0.0; grid.ndof()];
            let t = now();
            for _ in 0..5 {
                fem::apply(&grid, &elems, fixed, &x, &mut y);
            }
            out.push((format!("apply level {l} ({} el)", grid.nel()), t.elapsed().as_secs_f64() * 200.0));
            if l < self.levels.len() {
                let next = &self.levels[l];
                let mut c = vec![0.0; next.grid.ndof()];
                let t = now();
                restrict(&next.grid, &next.fixed, &grid, fixed, &x, &mut c);
                out.push((format!("restrict {l}->{}", l + 1), t.elapsed().as_secs_f64() * 1000.0));
                let t = now();
                prolong_add(&next.grid, &next.fixed, &grid, fixed, &c, &mut y);
                out.push((format!("prolong {}->{l}", l + 1), t.elapsed().as_secs_f64() * 1000.0));
            }
        }
        let n = self.levels.last().map(|l| l.grid.ndof()).unwrap_or(self.fine_grid.ndof());
        let b = vec![1.0; n];
        let mut x = vec![0.0; n];
        let t = now();
        self.coarse_solve(&b, &mut x);
        out.push((format!("coarse solve ({n} dofs)"), t.elapsed().as_secs_f64() * 1000.0));
        out
    }

    /// `z ≈ K⁻¹ · r`: one V-cycle from zero.
    pub fn vcycle(&self, ke: &Ke, emod: &[f64], r: &[f64], z: &mut [f64]) {
        self.cycle(0, ke, emod, r, z);
    }

    fn level_parts<'a>(&'a self, l: usize, ke: &'a Ke, emod: &'a [f64]) -> (Grid, Elements<'a>, &'a [bool], &'a [f64], f64) {
        if l == 0 {
            (self.fine_grid, Elements::Scaled { ke, emod }, &self.fine_fixed, &self.fine_inv_diag, self.fine_omega)
        } else {
            let lv = &self.levels[l - 1];
            (lv.grid, Elements::Own(&lv.mats), &lv.fixed, &lv.inv_diag, lv.omega)
        }
    }

    fn cycle(&self, l: usize, ke: &Ke, emod: &[f64], b: &[f64], x: &mut [f64]) {
        if l == self.levels.len() {
            self.coarse_solve(b, x);
            return;
        }
        let (grid, elems, fixed, inv_diag, omega) = self.level_parts(l, ke, emod);
        let n = b.len();
        let mut ax = vec![0.0; n];
        x.iter_mut().for_each(|v| *v = 0.0);
        let smooth = |x: &mut [f64], ax: &mut [f64]| {
            for _ in 0..self.settings.smooth_steps {
                fem::apply(&grid, &elems, fixed, x, ax);
                let ax: &[f64] = ax;
                crate::par::update(x, |i, v| v + omega * inv_diag[i] * (b[i] - ax[i]));
            }
        };
        smooth(x, &mut ax);
        // Residual, restricted.
        fem::apply(&grid, &elems, fixed, x, &mut ax);
        let mut res = vec![0.0; n];
        crate::par::fill(&mut res, |i| b[i] - ax[i]);
        let next = &self.levels[l];
        let mut rc = vec![0.0; next.grid.ndof()];
        restrict(&next.grid, &next.fixed, &grid, fixed, &res, &mut rc);
        let mut xc = vec![0.0; rc.len()];
        self.cycle(l + 1, ke, emod, &rc, &mut xc);
        prolong_add(&next.grid, &next.fixed, &grid, fixed, &xc, x);
        smooth(x, &mut ax);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn child_prolongation_partitions_unity() {
        for c in 0..8 {
            let p = child_prolongation(c);
            for a in 0..8 {
                let s: f64 = p[a].iter().sum();
                assert!((s - 1.0).abs() < 1e-15);
            }
        }
    }

    /// Element-wise Galerkin equals the assembled RAP of the fine operator (without supports).
    #[test]
    fn elementwise_galerkin_is_exact() {
        let fg = Grid::new(4, 2, 2, 1.0);
        let cg = fg.coarsen().unwrap();
        let ke = crate::ke::hex8_ke(0.3, 1.0);
        let emod: Vec<f64> = (0..fg.nel()).map(|e| 0.2 + 0.1 * e as f64).collect();
        // Supports only matter for the masks and the coarse factorization, not the element matrices.
        let mut fix = vec![false; fg.ndof()];
        for k in 0..=fg.nz {
            for j in 0..=fg.ny {
                for d in 0..3 {
                    fix[3 * fg.node(0, j, k) + d] = true;
                }
            }
        }
        let mut h = Hierarchy::new(fg, fix, MgSettings { max_coarse_dofs: 1, direct_dofs: 0, ..Default::default() });
        h.update(&ke, &emod).unwrap();
        let nofix_f = vec![false; fg.ndof()];
        let nofix_c = vec![false; cg.ndof()];
        let level = &h.levels[0];
        for col in 0..cg.ndof() {
            let mut ec = vec![0.0; cg.ndof()];
            ec[col] = 1.0;
            // RAP column
            let mut pf = vec![0.0; fg.ndof()];
            prolong_add(&cg, &nofix_c, &fg, &nofix_f, &ec, &mut pf);
            let mut kpf = vec![0.0; fg.ndof()];
            fem::apply(&fg, &Elements::Scaled { ke: &ke, emod: &emod }, &nofix_f, &pf, &mut kpf);
            let mut rap = vec![0.0; cg.ndof()];
            restrict(&cg, &nofix_c, &fg, &nofix_f, &kpf, &mut rap);
            // Element-wise column
            let mut ew = vec![0.0; cg.ndof()];
            fem::apply(&cg, &Elements::Own(&level.mats), &nofix_c, &ec, &mut ew);
            // Coarse matrices are stored in single precision.
            let scale = rap.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            for r in 0..cg.ndof() {
                assert!((rap[r] - ew[r]).abs() < 1e-6 * scale, "mismatch at ({r},{col}): {} vs {}", rap[r], ew[r]);
            }
        }
    }
}
