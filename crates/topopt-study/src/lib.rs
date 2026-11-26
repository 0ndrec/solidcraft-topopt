//! A topology-optimization study as the SolidCraft plugin runs it, without SolidCraft: the host
//! supplies ray casting ([`voxel::Caster`]), the original surface for pinning, and progress,
//! cancel and clock ([`Hooks`]). Native tests drive it with parry3d meshes; the plugin drives it
//! with `view.pick` and `model.mesh`.
//!
//! A run has two halves so that the heavy one can move out of the plugin sandbox:
//! - [`prepare`]: voxelize the bodies and build the problem (needs the CAD app's ray casting);
//! - [`solve`]: optimize and make the mesh (pure computation: the plugin's single thread, or the
//!   multi-threaded TopOpt Solver through the files of [`solver`]).

pub mod build;
pub mod solver;
pub mod spec;
pub mod status;
pub mod study;
pub mod voxel;

pub use study::{Condition, FaceRef, MATERIALS, Material, MeshSettings, RunResult, Settings, Study};
pub use topopt_mesh::Mesh;
pub use spec::{MeshCaster, Spec, prepare_spec};
pub use voxel::{Caster, Hit};

use serde::{Deserialize, Serialize};
use topopt_core::{ElementKind, Grid, Optimizer, Problem};
use topopt_mesh::{Field, Options, Projector, extract};

pub trait Hooks {
    fn progress(&mut self, fraction: f32, message: &str);
    fn cancelled(&mut self) -> bool;
    /// Seconds since the run started (0 when the host has no clock).
    fn elapsed(&mut self) -> f64;
}

/// The optimized densities, kept so the mesh can be remade with other settings.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Density {
    pub grid: Grid,
    pub values: Vec<f32>,
    pub pads: Vec<bool>,
    /// Cells inside the design body (and outside keep-outs).
    pub inside: Vec<bool>,
}

pub struct Outcome {
    pub mesh: Mesh,
    pub density: Density,
    pub result: RunResult,
}

/// The design body as `inspect.measure` gives it.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct BodyInfo {
    pub lo: [f64; 3],
    pub hi: [f64; 3],
    /// Exact volume, mm³ (0 when unknown: the voxel volume is used).
    pub volume: f64,
}

/// The original surface for pinning: positions and triangles (world mm).
pub type Surface<'a> = Option<(&'a [[f64; 3]], &'a [[u32; 3]])>;

/// Everything [`solve`] needs, independent of the CAD app: what the plugin hands to the solver.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Job {
    pub version: u32,
    pub id: String,
    pub study: Study,
    pub body: BodyInfo,
    /// Element kinds encode the voxels: void outside, solid pads, design elsewhere.
    pub problem: Problem,
    pub settings: topopt_core::Settings,
    pub warnings: Vec<String>,
    pub surface: Option<(Vec<[f64; 3]>, Vec<[u32; 3]>)>,
    /// Where the mesh goes, inside the `results` folder.
    pub mesh_path: String,
}

pub const JOB_VERSION: u32 = 1;

/// Voxelize the bodies and build the problem.
pub fn prepare(
    study: &Study,
    body: BodyInfo,
    design: &mut dyn Caster,
    keep_out: &mut [&mut dyn Caster],
    surface: Surface,
    hooks: &mut dyn Hooks,
) -> Result<Job, String> {
    let problems = study.problems();
    if !problems.is_empty() {
        return Err(format!("the study is not ready: {}", problems.join("; ")));
    }
    let s = &study.settings;
    let (grid, mut warnings) = voxel::plan_grid(body.lo, body.hi, s.element_size, s.max_elements);
    hooks.progress(0.0, &format!("voxelizing ({}×{}×{}, {:.2} mm)", grid.nx, grid.ny, grid.nz, grid.h));
    let design_vox = voxel::voxelize(design, &grid, |f| hooks.progress(0.04 * f, "voxelizing"))?;
    let mut blocked = vec![false; grid.nel()];
    for caster in keep_out.iter_mut() {
        let v = voxel::voxelize(*caster, &grid, |_| {})?;
        for (b, i) in blocked.iter_mut().zip(&v.inside) {
            *b |= *i;
        }
    }
    let built = build::build(study, grid, &design_vox, &blocked)?;
    warnings.extend(built.warnings.iter().cloned());
    finish_job(study, body, built.problem, warnings, surface)
}

/// The core settings for a problem (the volume target is a fraction of the whole body; pads are
/// solid on top of it), and the job that carries everything [`solve`] needs.
pub(crate) fn finish_job(study: &Study, body: BodyInfo, problem: Problem, mut warnings: Vec<String>, surface: Surface) -> Result<Job, String> {
    let s = &study.settings;
    let grid = problem.grid;
    let n_inside = problem.kinds.iter().filter(|k| **k != ElementKind::Void).count();
    let n_pads = problem.kinds.iter().filter(|k| **k == ElementKind::Solid).count();
    let n_design = n_inside - n_pads;
    let mut volfrac = (s.volfrac * n_inside as f64 - n_pads as f64) / n_design as f64;
    if !(0.05..=0.95).contains(&volfrac) {
        warnings.push(format!("the volume target leaves {:.0} % of the free space for material; clamped", volfrac * 100.0));
        volfrac = volfrac.clamp(0.05, 0.95);
    }
    let rmin = if s.min_member > 0.0 { (s.min_member / (2.0 * grid.h)).max(1.5) } else { 1.5 };
    let settings = topopt_core::Settings {
        volfrac,
        penal: s.penal,
        rmin,
        beta_max: s.beta_max,
        move_limit: s.move_limit,
        max_iter: s.max_iterations,
        ..Default::default()
    };
    Ok(Job {
        version: JOB_VERSION,
        id: String::new(),
        study: study.clone(),
        body,
        problem,
        settings,
        warnings,
        surface: surface.map(|(p, t)| (p.to_vec(), t.to_vec())),
        mesh_path: String::new(),
    })
}

/// Optimize and make the mesh.
pub fn solve(job: &Job, hooks: &mut dyn Hooks) -> Result<Outcome, String> {
    let study = &job.study;
    let s = &study.settings;
    let problem = &job.problem;
    let grid = problem.grid;
    let mut warnings = job.warnings.clone();
    // Defaults of the solver settings depend on where it runs (direct-solve thresholds differ
    // between the sandbox and native code), so they are not taken from the job.
    let mut settings = job.settings;
    settings.mg = topopt_core::mg::MgSettings::default();
    let mut opt = Optimizer::new(problem.clone(), settings)?;
    let max_iter = s.max_iterations.max(1);
    let mut stopped_by_time = false;
    let outcome = opt.run(|it| {
        hooks.progress(
            0.05 + 0.85 * it.iter as f32 / max_iter as f32,
            &format!("iteration {} — compliance {:.5}, β {}", it.iter, it.compliance, it.beta),
        );
        if hooks.elapsed() > s.time_limit {
            stopped_by_time = true;
            return false;
        }
        !hooks.cancelled()
    })?;
    if hooks.cancelled() {
        return Err("cancelled".into());
    }
    if stopped_by_time {
        warnings.push(format!("stopped at the {} s time limit before converging", s.time_limit));
    }
    let last = *opt.history.last().ok_or("no iteration ran")?;
    let inside: Vec<bool> = problem.kinds.iter().map(|k| *k != ElementKind::Void).collect();
    let pads: Vec<bool> = problem.kinds.iter().map(|k| *k == ElementKind::Solid).collect();
    let n_inside = inside.iter().filter(|&&i| i).count();
    let density = Density { grid, values: opt.x_phys.iter().map(|&v| v as f32).collect(), pads, inside };

    hooks.progress(0.9, "making the mesh");
    let surface = job.surface.as_ref().map(|(p, t)| (p.as_slice(), t.as_slice()));
    let (mesh, mesh_report) = make_mesh(study, &density, surface)?;
    let outcome_name = if stopped_by_time {
        "time limit"
    } else {
        match outcome {
            topopt_core::Outcome::Converged => "converged",
            topopt_core::Outcome::MaxIterations => "iteration limit",
            topopt_core::Outcome::Cancelled => "stopped",
        }
    };
    let body = job.body;
    let result = RunResult {
        outcome: outcome_name.into(),
        iterations: last.iter,
        compliance: last.compliance,
        volume_fraction: mesh_report.volume / if body.volume > 0.0 { body.volume } else { body_volume(&density) },
        mass_g: mesh_report.volume * study.material.density / 1000.0,
        mesh_volume_mm3: mesh_report.volume,
        triangles: mesh_report.triangles,
        genus: mesh_report.genus,
        elements: n_inside,
        grid: [grid.nx, grid.ny, grid.nz],
        h: grid.h,
        seconds: hooks.elapsed(),
        stl: job.mesh_path.clone(),
        inserted: false,
        warnings,
    };
    hooks.progress(1.0, "done");
    Ok(Outcome { mesh, density, result })
}

/// Voxelize, optimize, mesh — in one place.
pub fn run(
    study: &Study,
    body: BodyInfo,
    design: &mut dyn Caster,
    keep_out: &mut [&mut dyn Caster],
    surface: Surface,
    hooks: &mut dyn Hooks,
) -> Result<Outcome, String> {
    let job = prepare(study, body, design, keep_out, surface, hooks)?;
    solve(&job, hooks)
}

/// Volume of the design body as voxelized (mm³).
fn body_volume(d: &Density) -> f64 {
    d.inside.iter().filter(|&&i| i).count() as f64 * d.grid.h.powi(3)
}

/// The mesh of a density field with the study's mesh settings.
pub fn make_mesh(study: &Study, density: &Density, surface: Surface) -> Result<(Mesh, topopt_mesh::MeshReport), String> {
    let g = density.grid;
    let dims = [g.nx, g.ny, g.nz];
    let values: Vec<f64> = density.values.iter().map(|&v| v as f64).collect();
    let field = Field::from_cells(dims, g.h, g.origin, &values, 1);
    let anchors = field.cell_mask(dims, 1, &density.pads);
    let projector = match (study.mesh.keep_faces, surface) {
        (true, Some((p, t))) if !t.is_empty() => Some(Projector::new(p, t)?),
        _ => None,
    };
    let m = &study.mesh;
    let triangles = m.triangles.min(study::MAX_TRIANGLES);
    let opts = Options {
        iso: m.iso,
        smooth_iterations: m.smoothing,
        remove_islands: m.remove_islands,
        fill_cavities: m.fill_cavities,
        target_triangles: Some(triangles),
        ..Default::default()
    };
    let (mesh, report) = extract(field, Some(&anchors), projector.as_ref(), &opts)?;
    Ok((mesh, report.mesh))
}
