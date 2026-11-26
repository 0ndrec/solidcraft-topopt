//! Preconditioned conjugate gradients.

/// Outcome of a solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PcgStats {
    pub iterations: usize,
    /// ‖b − A·x‖ / ‖b‖ at exit.
    pub rel_residual: f64,
    pub converged: bool,
}

pub(crate) fn dot(a: &[f64], b: &[f64]) -> f64 {
    crate::par::dot(a, b)
}

/// Solves `A·x = b` from the initial `x` until `‖r‖ ≤ tol·‖b‖`. `apply(x, y)` sets `y = A·x`;
/// `precond(r, z)` sets `z ≈ A⁻¹·r` (symmetric positive definite).
pub fn pcg(
    mut apply: impl FnMut(&[f64], &mut [f64]),
    mut precond: impl FnMut(&[f64], &mut [f64]),
    b: &[f64],
    x: &mut [f64],
    tol: f64,
    max_iter: usize,
) -> PcgStats {
    let n = b.len();
    let bnorm = dot(b, b).sqrt();
    if bnorm == 0.0 {
        x.iter_mut().for_each(|v| *v = 0.0);
        return PcgStats { iterations: 0, rel_residual: 0.0, converged: true };
    }
    let mut ax = vec![0.0; n];
    apply(x, &mut ax);
    let mut r = vec![0.0; n];
    crate::par::fill(&mut r, |i| b[i] - ax[i]);
    let mut rel = dot(&r, &r).sqrt() / bnorm;
    if rel <= tol {
        return PcgStats { iterations: 0, rel_residual: rel, converged: true };
    }
    let mut z = vec![0.0; n];
    precond(&r, &mut z);
    let mut p = z.clone();
    let mut q = vec![0.0; n];
    let mut rz = dot(&r, &z);
    for it in 1..=max_iter {
        apply(&p, &mut q);
        let pq = dot(&p, &q);
        if pq <= 0.0 {
            return PcgStats { iterations: it, rel_residual: rel, converged: false };
        }
        let alpha = rz / pq;
        crate::par::update(x, |i, v| v + alpha * p[i]);
        crate::par::update(&mut r, |i, v| v - alpha * q[i]);
        rel = dot(&r, &r).sqrt() / bnorm;
        if rel <= tol {
            return PcgStats { iterations: it, rel_residual: rel, converged: true };
        }
        precond(&r, &mut z);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        crate::par::update(&mut p, |i, v| z[i] + beta * v);
    }
    PcgStats { iterations: max_iter, rel_residual: rel, converged: false }
}
