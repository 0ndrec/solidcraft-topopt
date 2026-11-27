//! Deflection and stress of a TopOpt job's design body, all of it solid: the baseline a result is
//! compared with.
//!
//!     cargo run --release -p topopt-study --example check_job -- JOB.json
//!
//! JOB.json is a job the plugin handed to TopOpt Solver (`results/topopt/solver/jobs/*.part0`
//! when it is one part) or one `topopt-solver` wrote.

use topopt_core::{ElementKind, Optimizer, stress};

fn main() -> Result<(), String> {
    let path = std::env::args().nth(1).ok_or("usage: check_job JOB.json")?;
    let job: topopt_study::Job = serde_json::from_slice(&std::fs::read(&path).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
    let problem = job.problem.clone();
    let solid: Vec<f64> = (0..problem.grid.nel()).map(|e| if problem.kind(e) == ElementKind::Void { 0.0 } else { 1.0 }).collect();
    let volume = solid.iter().sum::<f64>() * problem.grid.h.powi(3);
    let mut settings = job.settings;
    settings.mg = topopt_core::mg::MgSettings::default();
    settings.pcg_tol = 1e-8;
    let mut opt = Optimizer::new(problem, settings)?;
    println!("solid design body: {volume:.0} mm³ of voxels, {:.1} g", volume * job.study.material.density / 1000.0);
    for (k, c) in opt.check(&solid)?.iter().enumerate() {
        let material: Vec<f64> = c.von_mises.iter().zip(&solid).filter(|(_, s)| **s > 0.5).map(|(v, _)| *v).collect();
        println!(
            "case {}: compliance {:.4} N·mm, max deflection {:.4} mm, von Mises 99 % below {:.2} MPa, peak {:.2} MPa",
            k + 1,
            c.compliance,
            c.max_displacement,
            stress::quantile(&material, 0.99),
            material.iter().copied().fold(0.0, f64::max)
        );
    }
    Ok(())
}
