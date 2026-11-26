//! Matrix-free stiffness operators on a structured grid.
//!
//! The global matrix is never assembled. `K·u` loops over elements: gather the 24 element dofs,
//! multiply by the element matrix, scatter-add. Supports are applied symmetrically: fixed dofs are
//! taken as zero on input and the operator is the identity on them, so the system stays SPD.

use crate::grid::Grid;
use crate::ke::Ke;

/// Calls `f(e, dofs)` for every element in index order.
#[inline]
pub fn for_each_element(grid: &Grid, mut f: impl FnMut(usize, &[usize; 24])) {
    let mut e = 0;
    for k in 0..grid.nz {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                f(e, &grid.elem_dofs(i, j, k));
                e += 1;
            }
        }
    }
}

/// A coarse-level element matrix, stored in single precision: coarse levels only precondition,
/// and halving their size halves the memory each multigrid sweep streams.
pub type Ke32 = [f32; 576];

/// The element matrices of one grid level.
pub enum Elements<'a> {
    /// One shared matrix scaled per element (the finest level: `emod[e] · ke`).
    Scaled { ke: &'a Ke, emod: &'a [f64] },
    /// One matrix per element (coarse levels).
    Own(&'a [Ke32]),
}

impl Elements<'_> {
    /// Diagonal entry `a` of element `e`'s matrix.
    #[inline]
    fn diag(&self, e: usize, a: usize) -> f64 {
        match self {
            Elements::Scaled { ke, emod } => emod[e] * ke[a * 25],
            Elements::Own(m) => m[e][a * 25] as f64,
        }
    }
}

/// `y += s · M · x_e` for one element whose node 0 is `n0`: gather, 24×24 product, scatter.
///
/// # Safety
/// `n0 + off[a]` must be valid nodes of `xm` and `y` (3 dofs each), and no other thread may
/// write the same nodes concurrently.
#[inline(always)]
unsafe fn element<T: Copy + Into<f64>>(m: &[T; 576], s: f64, n0: usize, off: &[usize; 8], xm: &[f64], y: *mut f64) {
    unsafe {
        let mut xe = [0.0f64; 24];
        for a in 0..8 {
            let b = 3 * (n0 + off[a]);
            xe[3 * a] = *xm.get_unchecked(b);
            xe[3 * a + 1] = *xm.get_unchecked(b + 1);
            xe[3 * a + 2] = *xm.get_unchecked(b + 2);
        }
        let mut acc = [0.0f64; 24];
        for c in 0..24 {
            let xc = xe[c] * s;
            let col = m.get_unchecked(c * 24..c * 24 + 24);
            for r in 0..24 {
                acc[r] += (*col.get_unchecked(r)).into() * xc;
            }
        }
        for a in 0..8 {
            let b = 3 * (n0 + off[a]);
            *y.add(b) += acc[3 * a];
            *y.add(b + 1) += acc[3 * a + 1];
            *y.add(b + 2) += acc[3 * a + 2];
        }
    }
}

/// `y = K·x` with supports: `y[d] = x[d]` on fixed dofs.
///
/// The hot loop of the whole optimizer. Fixed dofs are zeroed once in a copy of `x`, so the
/// element loop has no branches; node indices advance incrementally along x; each node's three
/// dofs are contiguous.
///
/// With the `parallel` feature, elements are coloured by the parity of `(i, j, k)`: two elements
/// of one colour share no node, so the rows `(j, k)` of a colour are summed into `y` by different
/// threads without locks, one colour after another.
pub fn apply(grid: &Grid, elems: &Elements, fixed: &[bool], x: &[f64], y: &mut [f64]) {
    match elems {
        Elements::Scaled { ke, emod } => apply_with(grid, |e| (*ke, emod[e]), fixed, x, y),
        Elements::Own(m) => apply_with(grid, |e| (&m[e], 1.0), fixed, x, y),
    }
}

fn apply_with<'m, T: Copy + Into<f64> + Sync + 'm>(
    grid: &Grid,
    get: impl Fn(usize) -> (&'m [T; 576], f64) + Sync,
    fixed: &[bool],
    x: &[f64],
    y: &mut [f64],
) {
    let n = grid.ndof();
    assert!(x.len() == n && y.len() == n && fixed.len() == n);
    let mut xm = vec![0.0; n];
    crate::par::fill(&mut xm, |d| if fixed[d] { 0.0 } else { x[d] });
    crate::par::fill(y, |_| 0.0);
    let sx = grid.nx + 1;
    let sxy = sx * (grid.ny + 1);
    // Corner offsets in nodes, in CORNERS order.
    let off = [0, 1, 1 + sx, sx, sxy, 1 + sxy, 1 + sx + sxy, sx + sxy];
    let yp = crate::par::SharedMut(y.as_mut_ptr());
    if crate::par::threads() > 1 {
        for colour in 0..8usize {
            let (ci, cj, ck) = (colour & 1, (colour >> 1) & 1, (colour >> 2) & 1);
            let rows: Vec<(usize, usize)> =
                (ck..grid.nz).step_by(2).flat_map(|k| (cj..grid.ny).step_by(2).map(move |j| (j, k))).collect();
            let xm = &xm;
            crate::par::each(rows.len(), |r| {
                let (j, k) = rows[r];
                let row_base = grid.node(0, j, k);
                let e_base = grid.elem(0, j, k);
                let yp = yp;
                for i in (ci..grid.nx).step_by(2) {
                    let (m, s) = get(e_base + i);
                    // SAFETY: indices are in range (asserted sizes); elements of one colour share
                    // no node, and each row of this colour is handled by exactly one thread.
                    unsafe { element(m, s, row_base + i, &off, xm, yp.0) };
                }
            });
        }
    } else {
        let mut e = 0;
        for k in 0..grid.nz {
            for j in 0..grid.ny {
                let row_base = grid.node(0, j, k);
                for i in 0..grid.nx {
                    let (m, s) = get(e);
                    e += 1;
                    // SAFETY: every corner node index is < nnodes for i < nx, j < ny, k < nz.
                    unsafe { element(m, s, row_base + i, &off, &xm, yp.0) };
                }
            }
        }
    }
    for (d, &f) in fixed.iter().enumerate() {
        if f {
            y[d] = x[d];
        }
    }
}

/// The diagonal of K with supports (1 on fixed dofs).
pub fn diagonal(grid: &Grid, elems: &Elements, fixed: &[bool]) -> Vec<f64> {
    let mut diag = vec![0.0; grid.ndof()];
    for_each_element(grid, |e, dofs| {
        for a in 0..24 {
            diag[dofs[a]] += elems.diag(e, a);
        }
    });
    for (d, &f) in fixed.iter().enumerate() {
        if f || diag[d] == 0.0 {
            diag[d] = 1.0;
        }
    }
    diag
}

/// `out[e] = u_eᵀ · ke · u_e` for every element (the unscaled element strain energy, doubled).
pub fn element_energies(grid: &Grid, ke: &Ke, u: &[f64], out: &mut [f64]) {
    let plane = grid.nx * grid.ny;
    crate::par::chunks(out, plane, |k, slab| {
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let dofs = grid.elem_dofs(i, j, k);
                let mut ue = [0.0f64; 24];
                for a in 0..24 {
                    ue[a] = u[dofs[a]];
                }
                let mut sum = 0.0;
                for r in 0..24 {
                    let row = &ke[r * 24..r * 24 + 24];
                    let mut acc = 0.0;
                    for c in 0..24 {
                        acc += row[c] * ue[c];
                    }
                    sum += ue[r] * acc;
                }
                slab[i + grid.nx * j] = sum;
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ke::hex8_ke;

    /// The (coloured, possibly parallel) product equals a plain element-by-element sum.
    #[test]
    fn coloured_product_matches_reference() {
        let g = Grid::new(9, 6, 5, 1.0);
        let ke = hex8_ke(0.3, g.h);
        let emod: Vec<f64> = (0..g.nel()).map(|e| 1e-3 + ((e * 7919) % 101) as f64 / 101.0).collect();
        let mut fixed = vec![false; g.ndof()];
        fixed[5] = true;
        let x: Vec<f64> = (0..g.ndof()).map(|d| ((d * 31) % 17) as f64 - 8.0).collect();
        let mut y = vec![0.0; g.ndof()];
        apply(&g, &Elements::Scaled { ke: &ke, emod: &emod }, &fixed, &x, &mut y);
        let mut r = vec![0.0; g.ndof()];
        for_each_element(&g, |e, dofs| {
            for a in 0..24 {
                for b in 0..24 {
                    let xb = if fixed[dofs[b]] { 0.0 } else { x[dofs[b]] };
                    r[dofs[a]] += emod[e] * ke[a * 24 + b] * xb;
                }
            }
        });
        r[5] = x[5];
        for d in 0..g.ndof() {
            assert!((y[d] - r[d]).abs() < 1e-9 * (1.0 + r[d].abs()), "dof {d}: {} vs {}", y[d], r[d]);
        }
    }

    /// The matrix-free product equals an explicitly assembled matrix, and is symmetric.
    #[test]
    fn matches_assembled_and_symmetric() {
        let g = Grid::new(3, 2, 2, 1.5);
        let ke = hex8_ke(0.3, g.h);
        let emod: Vec<f64> = (0..g.nel()).map(|e| 0.1 + e as f64 * 0.05).collect();
        let mut fixed = vec![false; g.ndof()];
        for d in [0, 1, 2, 7, 11] {
            fixed[d] = true;
        }
        let n = g.ndof();
        let elems = Elements::Scaled { ke: &ke, emod: &emod };
        // Column j of the operator.
        let mut cols = vec![vec![0.0; n]; n];
        for j in 0..n {
            let mut x = vec![0.0; n];
            x[j] = 1.0;
            apply(&g, &elems, &fixed, &x, &mut cols[j]);
        }
        for i in 0..n {
            for j in 0..n {
                assert!((cols[j][i] - cols[i][j]).abs() < 1e-10, "not symmetric at {i},{j}");
            }
        }
        let diag = diagonal(&g, &elems, &fixed);
        for i in 0..n {
            assert!((diag[i] - cols[i][i]).abs() < 1e-10);
        }
    }
}
