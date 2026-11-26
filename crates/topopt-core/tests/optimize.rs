//! Optimizer checks: gradient against finite differences, and a full run on a small cantilever.

use topopt_core::{ElementKind, Grid, Optimizer, Outcome, Problem, Settings};

/// top3d's cantilever: x = 0 face clamped, a line load along y at the free end's bottom edge.
fn cantilever(nx: usize, ny: usize, nz: usize) -> Problem {
    let g = Grid::new(nx, ny, nz, 1.0);
    let mut pr = Problem::new(g);
    let (l, w, t) = (nx as f64, ny as f64, nz as f64);
    pr.fix_box([0.0, 0.0, 0.0], [0.0, w, t], [true; 3]);
    pr.load_box([l, 0.0, 0.0], [l, w, 0.0], [0.0, 0.0, -1.0]);
    pr.material.e = 1.0;
    pr.material.nu = 0.3;
    pr
}

#[test]
fn gradient_matches_finite_differences() {
    let mut pr = cantilever(8, 4, 4);
    // A passive solid pad and a void pocket, to cover both.
    pr.set_kind_box([0.0, 0.0, 3.0], [2.0, 4.0, 4.0], ElementKind::Solid);
    pr.set_kind_box([5.0, 1.0, 1.0], [7.0, 3.0, 3.0], ElementKind::Void);
    for projection in [false, true] {
        let settings = Settings { projection, beta_start: 4.0, pcg_tol: 1e-13, pcg_max_iter: 500, ..Default::default() };
        let mut opt = Optimizer::new(pr.clone(), settings).unwrap();
        let nel = pr.grid.nel();
        let x0: Vec<f64> = (0..nel).map(|e| 0.2 + 0.6 * ((e * 2654435761usize) % 997) as f64 / 997.0).collect();
        opt.set_design(&x0);
        let ev = opt.evaluate().unwrap();
        let mut checked = 0;
        for e in (0..nel).step_by(7) {
            if pr.kind(e) != ElementKind::Design {
                assert_eq!(ev.dc[e], 0.0);
                continue;
            }
            let h = 1e-6;
            let mut xp = x0.clone();
            xp[e] += h;
            opt.set_design(&xp);
            let cp = opt.evaluate().unwrap().compliance;
            let mut xm = x0.clone();
            xm[e] -= h;
            opt.set_design(&xm);
            let cm = opt.evaluate().unwrap().compliance;
            let fd = (cp - cm) / (2.0 * h);
            let rel = (fd - ev.dc[e]).abs() / fd.abs().max(1e-12);
            assert!(rel < 1e-4, "projection {projection}, element {e}: analytic {} vs finite difference {fd} (rel {rel:e})", ev.dc[e]);
            checked += 1;
        }
        assert!(checked > 5);
    }
}

#[test]
fn small_cantilever_converges_to_a_crisp_design() {
    let pr = cantilever(24, 12, 4);
    let settings = Settings { volfrac: 0.3, max_iter: 200, beta_every: 25, ..Default::default() };
    let mut opt = Optimizer::new(pr, settings).unwrap();
    let mut log = Vec::new();
    let outcome = opt.run(|it| {
        log.push(*it);
        true
    });
    let outcome = outcome.unwrap();
    let first = log.first().unwrap();
    let last = log.last().unwrap();
    println!(
        "{outcome:?} after {} iterations: compliance {:.3} → {:.3}, volume {:.4}, grayness {:.3}, beta {}, CG/iter ≈ {}",
        log.len(),
        first.compliance,
        last.compliance,
        last.volume,
        last.grayness,
        last.beta,
        log.iter().map(|i| i.cg_iterations).sum::<usize>() / log.len()
    );
    assert_ne!(outcome, Outcome::Cancelled);
    assert!((last.volume - 0.3).abs() < 2e-3, "volume {}", last.volume);
    assert!(last.compliance < 0.6 * first.compliance, "compliance did not drop enough");
    assert!(last.grayness < 0.15, "grayness {}", last.grayness);
    assert!(log.iter().all(|i| i.cg_worst_residual <= 1e-6), "a solve did not converge");
}
