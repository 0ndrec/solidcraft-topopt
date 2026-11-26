//! Solver checks: multigrid-preconditioned CG against plain CG, and a cantilever against beam theory.

use topopt_core::fem::{self, Elements};
use topopt_core::ke::hex8_ke;
use topopt_core::mg::{Hierarchy, MgSettings};
use topopt_core::pcg::pcg;
use topopt_core::{Grid, Problem};

/// Cantilever: x = 0 face clamped, total force `p` (−z) spread over the x = L face.
fn cantilever(nx: usize, ny: usize, nz: usize, h: f64, p: f64) -> Problem {
    let g = Grid::new(nx, ny, nz, h);
    let mut pr = Problem::new(g);
    let (l, w, t) = (nx as f64 * h, ny as f64 * h, nz as f64 * h);
    pr.fix_box([0.0, 0.0, 0.0], [0.0, w, t], [true; 3]);
    pr.load_box([l, 0.0, 0.0], [l, w, t], [0.0, 0.0, -p]);
    pr
}

/// Heterogeneous moduli with SIMP-like contrast, deterministic.
fn rough_moduli(n: usize) -> Vec<f64> {
    (0..n)
        .map(|e| {
            let x = ((e * 2654435761usize) % 1000) as f64 / 1000.0;
            1e-3 + x.powi(3)
        })
        .collect()
}

#[test]
fn multigrid_cg_matches_plain_cg_and_converges_fast() {
    let pr = cantilever(32, 16, 16, 1.0, 100.0);
    let g = pr.grid;
    let ke = hex8_ke(0.3, g.h);
    let fixed = pr.fixed_dofs();
    let f = &pr.force_vectors(&fixed)[0];
    for emod in [vec![1.0; g.nel()], rough_moduli(g.nel())] {
        let elems = Elements::Scaled { ke: &ke, emod: &emod };
        // A small coarsest level so three levels (two Galerkin products) are exercised.
        let mut h = Hierarchy::new(g, fixed.clone(), MgSettings { max_coarse_dofs: 1500, direct_dofs: 0, ..Default::default() });
        h.update(&ke, &emod).unwrap();
        assert!(h.depth() >= 3, "expected several levels, got {}", h.depth());

        let mut u_mg = vec![0.0; g.ndof()];
        let s_mg = pcg(|x, y| fem::apply(&g, &elems, &fixed, x, y), |r, z| h.vcycle(&ke, &emod, r, z), f, &mut u_mg, 1e-9, 200);
        println!("MG stats {s_mg:?}, depth {}", h.depth());
        assert!(s_mg.converged, "MG-CG did not converge: {s_mg:?}");

        let diag = fem::diagonal(&g, &elems, &fixed);
        let mut u_cg = vec![0.0; g.ndof()];
        let s_cg = pcg(
            |x, y| fem::apply(&g, &elems, &fixed, x, y),
            |r, z| z.iter_mut().zip(r.iter().zip(&diag)).for_each(|(zi, (ri, di))| *zi = ri / di),
            f,
            &mut u_cg,
            1e-11,
            50_000,
        );
        assert!(s_cg.converged, "Jacobi-CG did not converge: {s_cg:?}");
        let err = u_mg.iter().zip(&u_cg).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt()
            / u_cg.iter().map(|a| a * a).sum::<f64>().sqrt();
        println!("MG-CG {} iterations, Jacobi-CG {} iterations, relative difference {err:.2e}", s_mg.iterations, s_cg.iterations);
        assert!(err < 1e-6, "solutions differ by {err}");
        assert!(s_mg.iterations < 60, "MG-CG took {} iterations", s_mg.iterations);
    }
}

#[test]
fn cantilever_tip_deflection_matches_timoshenko() {
    // Steel-like beam 60 × 10 × 10 mm, 4 elements per mm of section height is too fine for a test;
    // h = 1 mm gives 10 elements through the height.
    let (e, nu, p) = (200_000.0, 0.3, 1000.0);
    let (l, b, t) = (60.0, 10.0, 10.0);
    let mut pr = cantilever(60, 10, 10, 1.0, p);
    pr.material.e = e;
    pr.material.nu = nu;
    let g = pr.grid;
    let mut ke = hex8_ke(nu, g.h);
    ke.iter_mut().for_each(|v| *v *= e);
    let fixed = pr.fixed_dofs();
    let f = &pr.force_vectors(&fixed)[0];
    let emod = vec![1.0; g.nel()];
    let elems = Elements::Scaled { ke: &ke, emod: &emod };
    let mut h = Hierarchy::new(g, fixed.clone(), MgSettings::default());
    h.update(&ke, &emod).unwrap();
    let mut u = vec![0.0; g.ndof()];
    let s = pcg(|x, y| fem::apply(&g, &elems, &fixed, x, y), |r, z| h.vcycle(&ke, &emod, r, z), f, &mut u, 1e-10, 300);
    assert!(s.converged);
    // Mean vertical displacement of the loaded face.
    let tip = pr.nodes_in_box([l, 0.0, 0.0], [l, b, t]);
    let uz = tip.iter().map(|&n| u[3 * n + 2]).sum::<f64>() / tip.len() as f64;
    let i = b * t.powi(3) / 12.0;
    let gmod = e / (2.0 * (1.0 + nu));
    let timoshenko = p * l.powi(3) / (3.0 * e * i) + p * l / (5.0 / 6.0 * gmod * b * t);
    let ratio = -uz / timoshenko;
    println!("tip deflection {:.5} mm, Timoshenko {:.5} mm, ratio {ratio:.4}, {} CG iterations", -uz, timoshenko, s.iterations);
    // Fully integrated Hex8 is slightly stiff in bending; the clamp adds a little stiffness too.
    assert!(ratio > 0.93 && ratio < 1.02, "ratio {ratio}");
}
