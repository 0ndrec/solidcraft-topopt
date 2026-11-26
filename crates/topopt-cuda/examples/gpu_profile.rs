//! Milliseconds per GPU operation on a cantilever grid.
//!
//!   cargo run --release -p topopt-cuda --example gpu_profile -- NX NY NZ

use topopt_core::backend::LinearSolver;
use topopt_core::ke::hex8_ke;
use topopt_core::mg::MgSettings;
use topopt_core::{Grid, Problem};

fn main() {
    let a: Vec<usize> = std::env::args().skip(1).filter_map(|s| s.parse().ok()).collect();
    let (nx, ny, nz) = (a.first().copied().unwrap_or(96), a.get(1).copied().unwrap_or(48), a.get(2).copied().unwrap_or(48));
    let g = Grid::new(nx, ny, nz, 1.0);
    let mut pr = Problem::new(g);
    pr.fix_box([0.0, 0.0, 0.0], [0.0, ny as f64, nz as f64], [true; 3]);
    let fixed = pr.fixed_dofs();
    let ke = hex8_ke(0.3, 1.0);
    let emod: Vec<f64> = (0..g.nel()).map(|e| if e % 3 == 0 { 1e-9 } else { 1.0 }).collect();
    let gpu = topopt_cuda::open(0).expect("GPU");
    let mut s = topopt_cuda::CudaSolver::new(gpu, g, fixed, MgSettings::default()).unwrap();
    s.update(&ke, &emod).unwrap();
    s.update(&ke, &emod).unwrap();
    println!("update: coarse operators {:.0} ms, upload {:.0} ms, smoothers {:.0} ms", s.update_ms[0], s.update_ms[1], s.update_ms[2]);
    println!("{} — {nx}x{ny}x{nz}, {} elements, levels {:?}", s.name(), g.nel(), s.levels().iter().map(|l| [l.nx, l.ny, l.nz]).collect::<Vec<_>>());
    for (what, ms) in s.profile(20).unwrap() {
        println!("{what:>28} {ms:9.3} ms");
    }
}
