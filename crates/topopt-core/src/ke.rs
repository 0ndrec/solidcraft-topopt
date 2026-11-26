//! Hex8 (trilinear) element stiffness for isotropic linear elasticity, by 2×2×2 Gauss quadrature.
//!
//! The matrix is for Young's modulus 1 and a cube of side `h`; the stiffness of a cube scales
//! linearly with its side, so `ke(nu, h) = h · ke(nu, 1)`. Local nodes follow
//! [`crate::grid::CORNERS`]; dofs are node-major (x, y, z per node).

use crate::grid::CORNERS;

/// Row-major 24×24.
pub type Ke = [f64; 576];

/// Element stiffness of a cube of side `h` with E = 1 and Poisson ratio `nu`.
pub fn hex8_ke(nu: f64, h: f64) -> Ke {
    assert!((-1.0..0.5).contains(&nu), "Poisson ratio must be in (-1, 0.5)");
    let lam = nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = 1.0 / (2.0 * (1.0 + nu));
    // Constitutive matrix, Voigt order xx, yy, zz, xy, yz, zx (engineering shear strains).
    let mut d = [[0.0f64; 6]; 6];
    for a in 0..3 {
        for b in 0..3 {
            d[a][b] = lam;
        }
        d[a][a] = lam + 2.0 * mu;
        d[a + 3][a + 3] = mu;
    }

    let g = 1.0 / 3f64.sqrt();
    let mut ke = [0.0f64; 576];
    // Natural coordinates of the corners, ±1.
    let xi: Vec<[f64; 3]> = CORNERS.iter().map(|c| [2.0 * c[0] as f64 - 1.0, 2.0 * c[1] as f64 - 1.0, 2.0 * c[2] as f64 - 1.0]).collect();
    // On a cube of side h the Jacobian is (h/2)·I: dN/dx = (2/h)·dN/dξ, det J = (h/2)³.
    let jac = 2.0 / h;
    let det = (h / 2.0).powi(3);
    for &gx in &[-g, g] {
        for &gy in &[-g, g] {
            for &gz in &[-g, g] {
                // Shape-function derivatives in x, y, z.
                let mut dn = [[0.0f64; 3]; 8];
                for (a, c) in xi.iter().enumerate() {
                    let (fx, fy, fz) = (1.0 + c[0] * gx, 1.0 + c[1] * gy, 1.0 + c[2] * gz);
                    dn[a] = [c[0] * fy * fz / 8.0 * jac, fx * c[1] * fz / 8.0 * jac, fx * fy * c[2] / 8.0 * jac];
                }
                // Strain-displacement matrix B (6×24).
                let mut bm = [[0.0f64; 24]; 6];
                for a in 0..8 {
                    let [dx, dy, dz] = dn[a];
                    let c = 3 * a;
                    bm[0][c] = dx;
                    bm[1][c + 1] = dy;
                    bm[2][c + 2] = dz;
                    bm[3][c] = dy;
                    bm[3][c + 1] = dx;
                    bm[4][c + 1] = dz;
                    bm[4][c + 2] = dy;
                    bm[5][c] = dz;
                    bm[5][c + 2] = dx;
                }
                // ke += Bᵀ D B · det (Gauss weights are 1).
                let mut db = [[0.0f64; 24]; 6];
                for r in 0..6 {
                    for c in 0..24 {
                        db[r][c] = (0..6).map(|s| d[r][s] * bm[s][c]).sum();
                    }
                }
                for r in 0..24 {
                    for c in 0..24 {
                        ke[r * 24 + c] += det * (0..6).map(|s| bm[s][r] * db[s][c]).sum::<f64>();
                    }
                }
            }
        }
    }
    ke
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetric_and_scales_with_h() {
        let k1 = hex8_ke(0.3, 1.0);
        let k2 = hex8_ke(0.3, 2.5);
        for r in 0..24 {
            for c in 0..24 {
                assert!((k1[r * 24 + c] - k1[c * 24 + r]).abs() < 1e-14);
                assert!((k2[r * 24 + c] - 2.5 * k1[r * 24 + c]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn rigid_modes_have_no_energy() {
        let k = hex8_ke(0.3, 1.0);
        let pos: Vec<[f64; 3]> = CORNERS.iter().map(|c| [c[0] as f64, c[1] as f64, c[2] as f64]).collect();
        // Three translations and three infinitesimal rotations.
        let modes: Vec<Box<dyn Fn([f64; 3]) -> [f64; 3]>> = vec![
            Box::new(|_| [1.0, 0.0, 0.0]),
            Box::new(|_| [0.0, 1.0, 0.0]),
            Box::new(|_| [0.0, 0.0, 1.0]),
            Box::new(|p| [-p[1], p[0], 0.0]),
            Box::new(|p| [0.0, -p[2], p[1]]),
            Box::new(|p| [p[2], 0.0, -p[0]]),
        ];
        for m in &modes {
            let u: Vec<f64> = pos.iter().flat_map(|&p| m(p)).collect();
            for r in 0..24 {
                let f: f64 = (0..24).map(|c| k[r * 24 + c] * u[c]).sum();
                assert!(f.abs() < 1e-12, "rigid mode produced force {f}");
            }
        }
    }

    #[test]
    fn uniaxial_stretch_energy_is_exact() {
        // u_x = ε·x with lateral contraction u_y = -νε·y, u_z = -νε·z is a uniform uniaxial stress
        // σ_xx = E·ε; strain energy of the unit cube = ½·E·ε².
        let nu = 0.3;
        let k = hex8_ke(nu, 1.0);
        let eps = 1e-3;
        let u: Vec<f64> = CORNERS
            .iter()
            .flat_map(|c| [eps * c[0] as f64, -nu * eps * c[1] as f64, -nu * eps * c[2] as f64])
            .collect();
        let energy: f64 = 0.5 * (0..24).map(|r| u[r] * (0..24).map(|c| k[r * 24 + c] * u[c]).sum::<f64>()).sum::<f64>();
        assert!((energy - 0.5 * eps * eps).abs() < 1e-15, "energy {energy}");
    }
}
