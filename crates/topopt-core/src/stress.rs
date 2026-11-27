//! Stresses and displacements of a solved design: what a result is checked against the
//! material with.
//!
//! The strain of a Hex8 element is evaluated at its centre, where the trilinear shape functions'
//! derivatives are ±1/(4h) per corner; the stress is the element's modulus times the isotropic
//! constitutive matrix times that strain. On a voxel model this is an estimate: stair-stepped
//! surfaces, corners and point-like supports raise single elements' stresses, so look at a high
//! percentile rather than the single largest value.

use crate::fem::for_each_element;
use crate::grid::{CORNERS, Grid};

/// Von Mises stress at each element's centre, MPa. `e`: the solid's Young's modulus (MPa);
/// `emod`: each element's modulus relative to it (as `Optimizer::moduli` gives); `u`: the
/// displacements, mm.
pub fn von_mises(grid: &Grid, e: f64, nu: f64, emod: &[f64], u: &[f64]) -> Vec<f64> {
    let lam = nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = 1.0 / (2.0 * (1.0 + nu));
    let s = 1.0 / (4.0 * grid.h);
    let sign: Vec<[f64; 3]> = CORNERS.iter().map(|c| c.map(|v| if v == 1 { s } else { -s })).collect();
    let mut out = vec![0.0; grid.nel()];
    for_each_element(grid, |el, dofs| {
        // Engineering strains: xx, yy, zz, xy, yz, zx.
        let mut eps = [0.0f64; 6];
        for (a, d) in sign.iter().enumerate() {
            let (ux, uy, uz) = (u[dofs[3 * a]], u[dofs[3 * a + 1]], u[dofs[3 * a + 2]]);
            eps[0] += d[0] * ux;
            eps[1] += d[1] * uy;
            eps[2] += d[2] * uz;
            eps[3] += d[1] * ux + d[0] * uy;
            eps[4] += d[2] * uy + d[1] * uz;
            eps[5] += d[0] * uz + d[2] * ux;
        }
        let k = e * emod[el];
        let tr = eps[0] + eps[1] + eps[2];
        let sx = k * (lam * tr + 2.0 * mu * eps[0]);
        let sy = k * (lam * tr + 2.0 * mu * eps[1]);
        let sz = k * (lam * tr + 2.0 * mu * eps[2]);
        let (txy, tyz, tzx) = (k * mu * eps[3], k * mu * eps[4], k * mu * eps[5]);
        out[el] = (0.5 * ((sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2)) + 3.0 * (txy * txy + tyz * tyz + tzx * tzx)).sqrt();
    });
    out
}

/// The largest nodal displacement, mm.
pub fn max_displacement(u: &[f64]) -> f64 {
    u.chunks_exact(3).map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()).fold(0.0, f64::max)
}

/// The `q`-quantile (0…1) of `values` (nearest rank); 0 when empty.
pub fn quantile(values: &[f64], q: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut v = values.to_vec();
    v.sort_by(f64::total_cmp);
    v[((q.clamp(0.0, 1.0) * v.len() as f64).ceil() as usize).clamp(1, v.len()) - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniaxial_stress_is_exact() {
        // u_x = ε·x, u_y = −νε·y, u_z = −νε·z is uniaxial tension σ_xx = E·ε everywhere.
        let grid = Grid::new(3, 2, 2, 0.5);
        let (e, nu, eps) = (2100.0, 0.38, 1e-3);
        let mut u = vec![0.0; grid.ndof()];
        for k in 0..=grid.nz {
            for j in 0..=grid.ny {
                for i in 0..=grid.nx {
                    let n = grid.node(i, j, k);
                    let p = grid.node_pos(n);
                    u[3 * n] = eps * p[0];
                    u[3 * n + 1] = -nu * eps * p[1];
                    u[3 * n + 2] = -nu * eps * p[2];
                }
            }
        }
        let vm = von_mises(&grid, e, nu, &vec![1.0; grid.nel()], &u);
        for s in vm {
            assert!((s - e * eps).abs() < 1e-9, "{s}");
        }
    }

    #[test]
    fn pure_shear_gives_root_three_tau() {
        let grid = Grid::new(2, 2, 2, 1.0);
        let (e, nu, gamma) = (1000.0, 0.3, 2e-3);
        let mut u = vec![0.0; grid.ndof()];
        for k in 0..=grid.nz {
            for j in 0..=grid.ny {
                for i in 0..=grid.nx {
                    u[3 * grid.node(i, j, k)] = gamma * grid.node_pos(grid.node(i, j, k))[1];
                }
            }
        }
        let tau = e / (2.0 * (1.0 + nu)) * gamma;
        for s in von_mises(&grid, e, nu, &vec![1.0; grid.nel()], &u) {
            assert!((s - 3f64.sqrt() * tau).abs() < 1e-9, "{s}");
        }
    }

    #[test]
    fn quantiles() {
        let v: Vec<f64> = (1..=100).map(|i| i as f64).collect();
        assert_eq!(quantile(&v, 0.99), 99.0);
        assert_eq!(quantile(&v, 1.0), 100.0);
        assert_eq!(quantile(&[], 0.5), 0.0);
    }
}
