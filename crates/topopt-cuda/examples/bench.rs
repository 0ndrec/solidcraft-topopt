//! CPU against GPU on the cantilever: one optimization iteration's update and solve, per size.
//!
//!   cargo run --release -p topopt-cuda --example bench -- [NX NY NZ ...]   (default: a size ladder)

use std::time::Instant;
use topopt_core::backend::LinearSolver;
use topopt_core::ke::hex8_ke;
use topopt_core::mg::MgSettings;
use topopt_core::{CpuSolver, Grid, Material, Problem};
use topopt_cuda::CudaSolver;

fn main() {
    let a: Vec<usize> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let sizes: Vec<[usize; 3]> = if a.len() >= 3 {
        a.chunks(3).filter(|c| c.len() == 3).map(|c| [c[0], c[1], c[2]]).collect()
    } else {
        vec![[40, 20, 20], [64, 32, 32], [96, 48, 48], [128, 64, 64]]
    };
    let gpu = topopt_cuda::open(0).expect("GPU");
    println!("{} — CPU: {} threads", gpu.info.name, topopt_core::par::threads());
    println!("{:>14} {:>9} | {:>9} {:>9} {:>5} | {:>9} {:>9} {:>5} | {:>6}", "grid", "elements", "CPU upd", "CPU solve", "its", "GPU upd", "GPU solve", "its", "speed");
    for [nx, ny, nz] in sizes {
        let g = Grid::new(nx, ny, nz, 1.0);
        let mut pr = Problem::new(g);
        pr.fix_box([0.0, 0.0, 0.0], [0.0, ny as f64, nz as f64], [true; 3]);
        pr.load_box([nx as f64, 0.0, 0.0], [nx as f64, ny as f64, 0.0], [0.0, 0.0, -1.0]);
        pr.material = Material { e: 1.0, nu: 0.3 };
        let fixed = pr.fixed_dofs();
        let f = pr.force_vectors(&fixed).remove(0);
        let ke = hex8_ke(0.3, 1.0);
        // A converged-looking design: solid and void (1e-9) in slabs, some grey.
        let emod: Vec<f64> = (0..g.nel())
            .map(|e| {
                let [i, j, k] = g.elem_ijk(e);
                let solid = k < nz / 4 || k >= nz - nz / 4 || (i + j) % 8 < 3 || j < 2 || j >= ny - 2;
                if solid { 1.0 } else if (i * 7 + k * 3) % 11 == 0 { 0.1 } else { 1e-9 }
            })
            .collect();
        let mg = MgSettings::default();
        let mut row = Vec::new();
        let solvers: [Box<dyn Fn() -> Result<Box<dyn LinearSolver>, String>>; 2] = [
            Box::new(|| Ok(Box::new(CpuSolver::new(g, fixed.clone(), mg)) as Box<dyn LinearSolver>)),
            Box::new(|| Ok(Box::new(CudaSolver::new(gpu.clone(), g, fixed.clone(), mg)?) as Box<dyn LinearSolver>)),
        ];
        for make in &solvers {
            let mut s = match make() {
                Ok(s) => s,
                Err(e) => {
                    println!("{nx}x{ny}x{nz}: {e}");
                    row.clear();
                    break;
                }
            };
            // Warm up (first update: JIT, allocation, power iteration from scratch), then time a
            // typical iteration: update and a solve from the previous answer after a small change.
            s.update(&ke, &emod).unwrap();
            let mut u = vec![0.0; g.ndof()];
            s.solve(&ke, &emod, &f, &mut u, 1e-6, 500).unwrap();
            let emod2: Vec<f64> = emod.iter().enumerate().map(|(e, &m)| if e % 13 == 0 { (m * 1.5).min(1.0) } else { m }).collect();
            let t = Instant::now();
            s.update(&ke, &emod2).unwrap();
            let upd = t.elapsed().as_secs_f64() * 1e3;
            let t = Instant::now();
            let st = s.solve(&ke, &emod2, &f, &mut u, 1e-6, 500).unwrap();
            let sol = t.elapsed().as_secs_f64() * 1e3;
            assert!(st.converged, "{st:?}");
            row.push((upd, sol, st.iterations));
        }
        if row.len() == 2 {
            let (c, g2) = (row[0], row[1]);
            println!(
                "{:>14} {:>9} | {:>7.0}ms {:>7.0}ms {:>5} | {:>7.0}ms {:>7.0}ms {:>5} | {:>5.1}x",
                format!("{nx}x{ny}x{nz}"),
                g.nel(),
                c.0,
                c.1,
                c.2,
                g2.0,
                g2.1,
                g2.2,
                (c.0 + c.1) / (g2.0 + g2.1)
            );
        }
    }
}
