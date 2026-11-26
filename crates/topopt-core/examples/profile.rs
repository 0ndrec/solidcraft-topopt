//! Where an iteration's time goes, on the cantilever benchmark.
//!
//!   cargo run --release -p topopt-core --example profile -- NX NY NZ [SMOOTH_STEPS [DIRECT_DOFS [STEPS [MAX_COARSE_DOFS]]]]

use std::time::Instant;
use topopt_core::fem::{self, Elements};
use topopt_core::filter::DensityFilter;
use topopt_core::ke::hex8_ke;
use topopt_core::mg::{Hierarchy, MgSettings};
use topopt_core::pcg::pcg;
use topopt_core::{Grid, Optimizer, Problem, Settings};

fn main() {
    let a: Vec<usize> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let (nx, ny, nz) = (a.first().copied().unwrap_or(64), a.get(1).copied().unwrap_or(32), a.get(2).copied().unwrap_or(32));
    let mut mg = MgSettings { smooth_steps: a.get(3).copied().unwrap_or(2), ..Default::default() };
    if let Some(&c) = a.get(6) {
        mg.max_coarse_dofs = c;
    }
    if let Some(&d) = a.get(4) {
        mg.direct_dofs = d;
    }
    let g = Grid::new(nx, ny, nz, 1.0);
    let mut pr = Problem::new(g);
    pr.fix_box([0.0, 0.0, 0.0], [0.0, ny as f64, nz as f64], [true; 3]);
    pr.load_box([nx as f64, 0.0, 0.0], [nx as f64, ny as f64, 0.0], [0.0, 0.0, -1.0]);
    pr.material = topopt_core::Material { e: 1.0, nu: 0.3 };
    let ke = hex8_ke(0.3, 1.0);
    let fixed = pr.fixed_dofs();
    let f = pr.force_vectors(&fixed).remove(0);
    let emod: Vec<f64> = (0..g.nel()).map(|e| 1e-9 + (0.3f64 + 0.4 * ((e * 2654435761usize) % 1000) as f64 / 1000.0).powi(3)).collect();
    println!("grid {nx}x{ny}x{nz}: {} elements, {} dofs", g.nel(), g.ndof());

    let t = Instant::now();
    let mut h = Hierarchy::new(g, fixed.clone(), mg);
    println!("hierarchy new        {:8.1} ms  levels {:?}", ms(t), h.level_grids().iter().map(|l| [l.nx, l.ny, l.nz]).collect::<Vec<_>>());
    let t = Instant::now();
    h.update(&ke, &emod).unwrap();
    println!("update (first)       {:8.1} ms", ms(t));
    let t = Instant::now();
    h.update(&ke, &emod).unwrap();
    println!("update               {:8.1} ms", ms(t));

    for (what, ms) in h.timing_breakdown(&ke, &emod) {
        println!("  {what:<32} {ms:8.2} ms");
    }
    let elems = Elements::Scaled { ke: &ke, emod: &emod };
    let mut y = vec![0.0; g.ndof()];
    let t = Instant::now();
    for _ in 0..10 {
        fem::apply(&g, &elems, &fixed, &f, &mut y);
    }
    let matvec = ms(t) / 10.0;
    println!("fine matvec          {:8.2} ms", matvec);
    let mut z = vec![0.0; g.ndof()];
    let t = Instant::now();
    for _ in 0..5 {
        h.vcycle(&ke, &emod, &f, &mut z);
    }
    println!("V-cycle              {:8.2} ms  ({:.1} matvecs)", ms(t) / 5.0, ms(t) / 5.0 / matvec);
    let mut u = vec![0.0; g.ndof()];
    let t = Instant::now();
    let s = pcg(|x, y| fem::apply(&g, &elems, &fixed, x, y), |r, z| h.vcycle(&ke, &emod, r, z), &f, &mut u, 1e-6, 200);
    println!("PCG cold (1e-6)      {:8.1} ms  {} iterations", ms(t), s.iterations);

    let filt = DensityFilter::new(g, 1.5);
    let x = vec![0.3; g.nel()];
    let mut out = vec![0.0; g.nel()];
    let t = Instant::now();
    for _ in 0..10 {
        filt.apply(&x, &mut out);
    }
    println!("filter apply         {:8.2} ms", ms(t) / 10.0);

    let steps = a.get(5).copied().unwrap_or(6);
    let settings = Settings { max_iter: steps, mg, ..Default::default() };
    let mut opt = Optimizer::new(pr, settings).unwrap();
    let total = Instant::now();
    for _ in 0..steps {
        let t = Instant::now();
        let it = opt.step().unwrap();
        if steps <= 12 || it.iter % 10 == 0 {
            println!("optimizer step {:3}   {:8.1} ms  CG {:3}  beta {}", it.iter, ms(t), it.cg_iterations, it.beta);
        }
    }
    println!("{steps} steps          {:8.1} ms, compliance {:.6e}", ms(total), opt.history.last().unwrap().compliance);
}

fn ms(t: Instant) -> f64 {
    t.elapsed().as_secs_f64() * 1e3
}
