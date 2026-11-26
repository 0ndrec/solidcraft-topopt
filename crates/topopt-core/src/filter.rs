//! Density filter (linear "hat" weights, Bruns & Tortorelli) and smoothed Heaviside projection
//! (Wang, Lazarov & Sigmund 2011).
//!
//! `x̃_j = Σ_i w_ij·x_i / Σ_i w_ij` with `w_ij = max(0, r − ‖c_i − c_j‖)` over grid elements, the
//! radius `r` in element lengths. The filter sets a minimum feature size of about `2r` elements
//! and removes checkerboards.

use crate::grid::Grid;

pub struct DensityFilter {
    grid: Grid,
    /// Neighbour offsets and weights.
    stencil: Vec<([i32; 3], f64)>,
    /// Σ_i w_ij for each element j (boundary elements have fewer neighbours).
    sums: Vec<f64>,
}

impl DensityFilter {
    /// `radius` in element lengths; values ≤ 1 leave the field unchanged.
    pub fn new(grid: Grid, radius: f64) -> Self {
        let reach = (radius.ceil() as i32 - 1).max(0);
        let mut stencil = Vec::new();
        for dk in -reach..=reach {
            for dj in -reach..=reach {
                for di in -reach..=reach {
                    let dist = ((di * di + dj * dj + dk * dk) as f64).sqrt();
                    let w = radius - dist;
                    if w > 0.0 {
                        stencil.push(([di, dj, dk], w));
                    }
                }
            }
        }
        if stencil.is_empty() {
            stencil.push(([0, 0, 0], 1.0));
        }
        let mut f = DensityFilter { grid, stencil, sums: vec![0.0; grid.nel()] };
        let ones = vec![1.0; grid.nel()];
        let mut sums = vec![0.0; grid.nel()];
        f.weighted_sum(&ones, &mut sums);
        f.sums = sums;
        f
    }

    /// `out_j = Σ_i w_ij·v_i` over in-grid neighbours.
    fn weighted_sum(&self, v: &[f64], out: &mut [f64]) {
        let g = &self.grid;
        let (nx, ny, nz) = (g.nx as i32, g.ny as i32, g.nz as i32);
        crate::par::chunks(out, g.nx * g.ny, |k, slab| {
            let k = k as i32;
            for j in 0..ny {
                for i in 0..nx {
                    let mut acc = 0.0;
                    for &([di, dj, dk], w) in &self.stencil {
                        let (a, b, c) = (i + di, j + dj, k + dk);
                        if a >= 0 && a < nx && b >= 0 && b < ny && c >= 0 && c < nz {
                            acc += w * v[g.elem(a as usize, b as usize, c as usize)];
                        }
                    }
                    slab[(i + nx * j) as usize] = acc;
                }
            }
        });
    }

    /// `x̃ = filter(x)`.
    pub fn apply(&self, x: &[f64], out: &mut [f64]) {
        self.weighted_sum(x, out);
        for (o, s) in out.iter_mut().zip(&self.sums) {
            *o /= s;
        }
    }

    /// The chain rule through the filter: given `g = ∂f/∂x̃`, returns `∂f/∂x`.
    pub fn apply_transpose(&self, g: &[f64], out: &mut [f64]) {
        let scaled: Vec<f64> = g.iter().zip(&self.sums).map(|(v, s)| v / s).collect();
        // The weights are symmetric, so the transpose is the same weighted sum.
        self.weighted_sum(&scaled, out);
    }
}

/// Smoothed Heaviside: 0 below `eta`, 1 above, sharper as `beta` grows.
#[inline]
pub fn project(x: f64, beta: f64, eta: f64) -> f64 {
    let (a, b) = ((beta * eta).tanh(), (beta * (1.0 - eta)).tanh());
    (a + (beta * (x - eta)).tanh()) / (a + b)
}

/// d project / dx.
#[inline]
pub fn project_derivative(x: f64, beta: f64, eta: f64) -> f64 {
    let (a, b) = ((beta * eta).tanh(), (beta * (1.0 - eta)).tanh());
    let t = (beta * (x - eta)).tanh();
    beta * (1.0 - t * t) / (a + b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constant_field_is_kept_and_transpose_is_adjoint() {
        let g = Grid::new(5, 4, 3, 1.0);
        let f = DensityFilter::new(g, 2.2);
        let ones = vec![0.7; g.nel()];
        let mut out = vec![0.0; g.nel()];
        f.apply(&ones, &mut out);
        assert!(out.iter().all(|v| (v - 0.7).abs() < 1e-14));
        // <F x, y> = <x, Fᵀ y>
        let x: Vec<f64> = (0..g.nel()).map(|e| ((e * 37) % 11) as f64 / 11.0).collect();
        let y: Vec<f64> = (0..g.nel()).map(|e| ((e * 13) % 7) as f64 / 7.0).collect();
        let (mut fx, mut fty) = (vec![0.0; g.nel()], vec![0.0; g.nel()]);
        f.apply(&x, &mut fx);
        f.apply_transpose(&y, &mut fty);
        let lhs: f64 = fx.iter().zip(&y).map(|(a, b)| a * b).sum();
        let rhs: f64 = x.iter().zip(&fty).map(|(a, b)| a * b).sum();
        assert!((lhs - rhs).abs() < 1e-12);
    }

    #[test]
    fn projection_ends_and_derivative() {
        for beta in [1.0, 8.0, 64.0] {
            assert!(project(0.0, beta, 0.5).abs() < 1e-12);
            assert!((project(1.0, beta, 0.5) - 1.0).abs() < 1e-12);
            let (x, h) = (0.37, 1e-6);
            let fd = (project(x + h, beta, 0.5) - project(x - h, beta, 0.5)) / (2.0 * h);
            assert!((fd - project_derivative(x, beta, 0.5)).abs() < 1e-5 * beta);
        }
    }
}
