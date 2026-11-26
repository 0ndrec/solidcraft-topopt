//! The GPU solver against the CPU one. Tests that need a GPU print "skipped" and pass on machines
//! without one (or without an NVIDIA driver).

use std::sync::Arc;
use topopt_core::backend::LinearSolver;
use topopt_core::fem::{self, Elements};
use topopt_core::ke::hex8_ke;
use topopt_core::mg::{Hierarchy, MgSettings};
use topopt_core::{CpuSolver, Grid, Material, Optimizer, Problem, Settings};
use topopt_cuda::{CudaOptions, CudaSolver, factory};

fn gpu() -> Option<Arc<topopt_cuda::driver::Gpu>> {
    match topopt_cuda::open(0) {
        Ok(g) => Some(g),
        Err(e) => {
            eprintln!("skipped: {e}");
            None
        }
    }
}

/// A cantilever: clamped at x = 0, loaded down along the free end's bottom edge.
fn cantilever(nx: usize, ny: usize, nz: usize) -> Problem {
    let mut pr = Problem::new(Grid::new(nx, ny, nz, 1.0));
    pr.fix_box([0.0, 0.0, 0.0], [0.0, ny as f64, nz as f64], [true; 3]);
    pr.load_box([nx as f64, 0.0, 0.0], [nx as f64, ny as f64, 0.0], [0.0, 0.0, -1.0]);
    pr.material = Material { e: 1.0, nu: 0.3 };
    pr
}

/// SIMP-like moduli: mostly solid or void (1e-9), some grey.
fn moduli(n: usize) -> Vec<f64> {
    (0..n)
        .map(|e| {
            let h = (e.wrapping_mul(2654435761) >> 7) % 100;
            if h < 45 { 1.0 } else if h < 85 { 1e-9 } else { 1e-9 + (h as f64 / 100.0).powi(3) }
        })
        .collect()
}

fn rel_diff(a: &[f64], b: &[f64]) -> f64 {
    let d: f64 = a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum::<f64>().sqrt();
    d / b.iter().map(|y| y * y).sum::<f64>().sqrt()
}

fn mg_small() -> MgSettings {
    MgSettings { max_coarse_dofs: 200, direct_dofs: 0, ..Default::default() }
}

#[test]
fn ptx_matches_the_kernel_source() {
    let src = std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/kernels/lib.rs")).unwrap();
    let src: Vec<u8> = String::from_utf8(src).unwrap().replace("\r\n", "\n").into_bytes();
    let mut h: u64 = 0xcbf29ce484222325;
    for b in src {
        h = (h ^ b as u64).wrapping_mul(0x100000001b3);
    }
    let line = topopt_cuda::PTX.lines().nth(1).unwrap();
    assert_eq!(line, format!("// source fnv1a64: {h:016x}"), "kernels/lib.rs changed: run scripts/build-ptx.sh");
}

#[test]
fn the_operator_matches_the_cpu() {
    let Some(gpu) = gpu() else { return };
    let pr = cantilever(9, 6, 5);
    let g = pr.grid;
    let fixed = pr.fixed_dofs();
    let ke = hex8_ke(0.3, g.h);
    let emod = moduli(g.nel());
    let mut s = CudaSolver::new(gpu, g, fixed.clone(), mg_small()).unwrap();
    s.update(&ke, &emod).unwrap();
    let x: Vec<f64> = (0..g.ndof()).map(|d| ((d * 31) % 17) as f64 - 8.0).collect();
    let mut y = vec![0.0; g.ndof()];
    fem::apply(&g, &Elements::Scaled { ke: &ke, emod: &emod }, &fixed, &x, &mut y);
    let d64 = rel_diff(&s.apply_fine(&x, true).unwrap(), &y);
    let d32 = rel_diff(&s.apply_fine(&x, false).unwrap(), &y);
    println!("K·x: double {d64:.1e}, single {d32:.1e}");
    assert!(d64 < 1e-13, "double: {d64}");
    assert!(d32 < 1e-6, "single: {d32}");
}

#[test]
fn the_vcycle_matches_the_cpu() {
    let Some(gpu) = gpu() else { return };
    let pr = cantilever(16, 8, 8);
    let g = pr.grid;
    let fixed = pr.fixed_dofs();
    let ke = hex8_ke(0.3, g.h);
    let emod = moduli(g.nel());
    let mut cpu = Hierarchy::new(g, fixed.clone(), mg_small());
    cpu.update(&ke, &emod).unwrap();
    let mut s = CudaSolver::new(gpu, g, fixed.clone(), mg_small()).unwrap();
    s.update(&ke, &emod).unwrap();
    assert_eq!(s.levels(), cpu.level_grids());
    let r: Vec<f64> = (0..g.ndof()).map(|d| if fixed[d] { 0.0 } else { (((d * 7919) % 101) as f64 - 50.0) / 50.0 }).collect();
    let mut z = vec![0.0; g.ndof()];
    cpu.vcycle(&ke, &emod, &r, &mut z);
    let d = rel_diff(&s.vcycle(&r).unwrap(), &z);
    println!("V-cycle: {d:.1e} (levels {:?}, damping {:?})", s.levels().iter().map(|l| [l.nx, l.ny, l.nz]).collect::<Vec<_>>(), s.omegas());
    assert!(d < 1e-4, "{d}");
}

#[test]
fn solves_match_the_cpu() {
    let Some(gpu) = gpu() else { return };
    let pr = cantilever(32, 16, 16);
    let g = pr.grid;
    let fixed = pr.fixed_dofs();
    let f = pr.force_vectors(&fixed).remove(0);
    let ke = hex8_ke(0.3, g.h);
    let emod: Vec<f64> = moduli(g.nel()).iter().map(|&m| m.max(1e-3)).collect();
    let mg = MgSettings::default();
    let tol = 1e-8;
    let mut cpu = CpuSolver::new(g, fixed.clone(), mg);
    cpu.update(&ke, &emod).unwrap();
    let mut u_cpu = vec![0.0; g.ndof()];
    let sc = cpu.solve(&ke, &emod, &f, &mut u_cpu, tol, 500).unwrap();
    let mut s = CudaSolver::new(gpu, g, fixed.clone(), mg).unwrap();
    s.update(&ke, &emod).unwrap();
    let mut u = vec![0.0; g.ndof()];
    let sg = s.solve(&ke, &emod, &f, &mut u, tol, 500).unwrap();
    // The GPU's answer, checked in double precision on the CPU.
    let mut ku = vec![0.0; g.ndof()];
    fem::apply(&g, &Elements::Scaled { ke: &ke, emod: &emod }, &fixed, &u, &mut ku);
    let res = rel_diff(&ku, &f);
    let d = rel_diff(&u, &u_cpu);
    println!("CPU {} iterations; GPU {} iterations, residual {:.1e} (checked {res:.1e}); difference {d:.1e}", sc.iterations, sg.iterations, sg.rel_residual);
    assert!(sg.converged && res < 2.0 * tol, "{sg:?}, checked residual {res}");
    assert!(d < 1e-6, "{d}");
    // A warm start from the answer needs no iteration.
    let again = s.solve(&ke, &emod, &f, &mut u, tol, 500).unwrap();
    assert_eq!(again.iterations, 0);
}

#[test]
fn optimizations_match_the_cpu() {
    if gpu().is_none() {
        return;
    }
    let pr = cantilever(24, 12, 12);
    let settings = Settings { max_iter: 25, ..Default::default() };
    let gpu_factory = factory(CudaOptions { min_elements: 0, fallback: false, ..Default::default() });
    let mut a = Optimizer::new(pr.clone(), settings).unwrap();
    let mut b = Optimizer::with_solver(pr, settings, Some(&gpu_factory)).unwrap();
    println!("{} against {}", a.solver_name(), b.solver_name());
    for _ in 0..settings.max_iter {
        let (ia, ib) = (a.step().unwrap(), b.step().unwrap());
        let dc = (ia.compliance - ib.compliance).abs() / ia.compliance;
        assert!(dc < 1e-4, "iteration {}: compliance {} against {}", ia.iter, ia.compliance, ib.compliance);
        assert!(ib.cg_worst_residual <= settings.pcg_tol, "iteration {}: residual {}", ib.iter, ib.cg_worst_residual);
    }
    let dx = a.x_phys.iter().zip(&b.x_phys).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
    println!("largest density difference after {} iterations: {dx:.1e}", settings.max_iter);
    assert!(dx < 1e-3, "{dx}");
}

#[test]
fn falls_back_to_the_cpu_without_a_usable_gpu() {
    let make = factory(CudaOptions { device: 99, min_elements: 0, fallback: true });
    let pr = cantilever(4, 2, 2);
    let s = make(&pr.grid, &pr.fixed_dofs(), &MgSettings::default()).unwrap();
    assert!(s.name().starts_with("CPU"), "{}", s.name());
    let note = s.note().unwrap();
    assert!(note.contains("GPU"), "{note}");
    let strict = factory(CudaOptions { device: 99, min_elements: 0, fallback: false });
    assert!(strict(&pr.grid, &pr.fixed_dofs(), &MgSettings::default()).is_err());
}

#[test]
fn small_problems_stay_on_the_cpu() {
    let make = factory(CudaOptions::default());
    let pr = cantilever(10, 5, 5);
    let s = make(&pr.grid, &pr.fixed_dofs(), &MgSettings::default()).unwrap();
    assert!(s.name().starts_with("CPU") && s.note().is_none(), "{}", s.name());
}
