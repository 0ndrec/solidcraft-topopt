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
pub mod trim;
pub mod voxel;

pub use study::{Condition, FaceRef, MATERIALS, Material, MeshSettings, RunResult, SMOOTHERS, Settings, Study};
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
    /// The filtered densities, before the Heaviside projection (pads 1, outside 0). The surface
    /// is made from these: the projection makes `values` nearly 0 or 1 per cell, so their iso
    /// surface follows the voxels, while the filtered field varies smoothly and places the same
    /// surface (the projection's threshold is the iso value) between the cell centres. Empty in
    /// densities saved before it existed: then `values` are used.
    #[serde(default)]
    pub filtered: Vec<f32>,
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
    let filtered = (0..grid.nel())
        .map(|e| match problem.kind(e) {
            ElementKind::Solid => 1.0,
            ElementKind::Void => 0.0,
            ElementKind::Design => opt.x_filtered[e] as f32,
        })
        .collect();
    let mut density = Density { grid, values: opt.x_phys.iter().map(|&v| v as f32).collect(), pads, inside, filtered };

    // Thin pieces that add almost no stiffness, when asked for.
    let mut trimmed = 0.0;
    if study.mesh.trim_percent > 0.0 {
        hooks.progress(0.88, "trimming dead weight");
        let max_loss = study.mesh.trim_percent / 100.0;
        let trim = trim::TrimSettings { max_loss, piece_loss: max_loss / 2.0, ..Default::default() };
        let x: Vec<f64> = density.values.iter().map(|&v| v as f64).collect();
        let (out, rep) = trim::trim_dead_ends(problem, &settings, &x, study.mesh.iso, &trim)?;
        trim::apply(&mut density, &out);
        trimmed = rep.removed_fraction;
        if rep.removed_pieces > 0 {
            warnings.push(format!(
                "trimmed {} thin pieces ({:.1} % of the material) for {:.2} % of the stiffness",
                rep.removed_pieces,
                100.0 * rep.removed_fraction,
                100.0 * (rep.compliance_after / rep.compliance_before - 1.0)
            ));
        }
    }

    hooks.progress(0.9, "making the mesh");
    let surface = job.surface.as_ref().map(|(p, t)| (p.as_slice(), t.as_slice()));
    let (mesh, full_report) = make_mesh(study, &density, surface)?;
    warnings.extend(mesh_warnings(study, &full_report, surface.is_some()));
    let mesh_report = full_report.mesh;
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
        trimmed,
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

/// The mesh of a density field with the study's mesh settings, and the mesher's report.
pub fn make_mesh(study: &Study, density: &Density, surface: Surface) -> Result<(Mesh, topopt_mesh::Report), String> {
    let g = density.grid;
    let dims = [g.nx, g.ny, g.nz];
    // The filtered field when there is one (see `Density::filtered`), sampled `upsample` times
    // finer by cubic interpolation: round members instead of faceted ones.
    let source = if density.filtered.len() == density.values.len() { &density.filtered } else { &density.values };
    let values: Vec<f64> = source.iter().map(|&v| v as f64).collect();
    let k = study.mesh.upsample.clamp(1, 3);
    let (fine, mut cells) = topopt_mesh::upsample_cells(dims, &values, k);
    let inside = topopt_mesh::upsample_mask(dims, &density.inside, k);
    let pads = topopt_mesh::upsample_mask(dims, &density.pads, k);
    for (c, (&i, &p)) in cells.iter_mut().zip(inside.iter().zip(&pads)) {
        if p {
            *c = 1.0;
        } else if !i {
            *c = 0.0;
        }
    }
    // Fillets where the optimized material meets the kept parts (pads, keep-in regions).
    if study.mesh.fillets {
        let r = if study.mesh.fillet_radius > 0.0 { study.mesh.fillet_radius } else { 1.5 * g.h };
        let material: Vec<bool> = (0..cells.len()).map(|c| cells[c] >= study.mesh.iso && !pads[c] && inside[c]).collect();
        topopt_mesh::round_joints(fine, g.h / k as f64, &mut cells, study.mesh.iso, &material, &pads, &inside, r);
    }
    let field = Field::from_cells(fine, g.h / k as f64, g.origin, &cells, 1);
    let anchors = field.cell_mask(fine, 1, &pads);
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
        smoother: topopt_mesh::Smoother::parse(&m.smoother).unwrap_or_default(),
        preserve_volume: m.preserve_volume,
        sharp_edges: m.sharp_edges,
        // The pinning band is in samples: keep it one *optimization* cell wide on a finer field.
        pin_band: Options::default().pin_band * k as f64,
        ..Default::default()
    };
    extract(field, Some(&anchors), projector.as_ref(), &opts)
}

/// Warnings about a mesh made with the study's settings.
pub fn mesh_warnings(study: &Study, report: &topopt_mesh::Report, surface_given: bool) -> Vec<String> {
    let mut w = Vec::new();
    if study.mesh.keep_faces && study.mesh.sharp_edges && surface_given && !report.sharp_edges {
        w.push("the sharp edges could not be kept without the surface crossing itself; they are rounded".into());
    }
    w
}
