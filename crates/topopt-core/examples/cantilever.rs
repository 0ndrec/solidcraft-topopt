//! top3d's cantilever benchmark: x = 0 face clamped, unit line load (−z) along y at the free end's
//! bottom edge, E = 1, ν = 0.3, h = 1. Prints one JSON line per iteration and a summary.
//!
//!   cargo run --release -p topopt-core --example cantilever -- NX NY NZ [ITERS] [projection 0|1] [volfrac] [rmin] [smooth] [omega] [power_steps]

use std::time::Instant;
use topopt_core::{Grid, Optimizer, Problem, Settings};

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let num = |i: usize, d: f64| a.get(i).and_then(|s| s.parse::<f64>().ok()).unwrap_or(d);
    let (nx, ny, nz) = (num(0, 60.0) as usize, num(1, 20.0) as usize, num(2, 4.0) as usize);
    let iters = num(3, 300.0) as usize;
    let projection = num(4, 1.0) != 0.0;
    let volfrac = num(5, 0.3);
    let rmin = num(6, 1.5);

    let g = Grid::new(nx, ny, nz, 1.0);
    let mut pr = Problem::new(g);
    let (l, w, t) = (nx as f64, ny as f64, nz as f64);
    pr.fix_box([0.0, 0.0, 0.0], [0.0, w, t], [true; 3]);
    pr.load_box([l, 0.0, 0.0], [l, w, 0.0], [0.0, 0.0, -1.0]);
    pr.material.e = 1.0;
    pr.material.nu = 0.3;
    let mut settings = Settings { volfrac, rmin, projection, max_iter: iters, ..Default::default() };
    // Optional multigrid overrides: smooth steps, omega, warm power steps.
    settings.mg.smooth_steps = num(7, settings.mg.smooth_steps as f64) as usize;
    settings.mg.omega = num(8, settings.mg.omega);
    settings.mg.power_steps = num(9, settings.mg.power_steps as f64) as usize;
    settings.mg.max_coarse_dofs = num(10, settings.mg.max_coarse_dofs as f64) as usize;

    let t0 = Instant::now();
    let mut opt = Optimizer::new(pr, settings).expect("problem");
    eprintln!("grid {nx}x{ny}x{nz}: {} elements, {} dofs, MG levels {:?}", g.nel(), g.ndof(),
        opt.mg_levels().iter().map(|l| [l.nx, l.ny, l.nz]).collect::<Vec<_>>());
    let mut last = Instant::now();
    let outcome = opt
        .run(|it| {
            let ms = last.elapsed().as_secs_f64() * 1e3;
            last = Instant::now();
            println!(
                "{{\"iter\":{},\"c\":{:.9},\"vol\":{:.6},\"change\":{:.5},\"beta\":{},\"gray\":{:.4},\"cg\":{},\"ms\":{:.1}}}",
                it.iter, it.compliance, it.volume, it.change, it.beta, it.grayness, it.cg_iterations, ms
            );
            true
        })
        .expect("run");
    let h = opt.history.last().unwrap();
    eprintln!(
        "{outcome:?}: {} iterations in {:.2} s ({:.0} ms/iter), compliance {:.4}, volume {:.4}, grayness {:.4}",
        opt.history.len(),
        t0.elapsed().as_secs_f64(),
        t0.elapsed().as_secs_f64() * 1e3 / opt.history.len() as f64,
        h.compliance,
        h.volume,
        h.grayness
    );
}
