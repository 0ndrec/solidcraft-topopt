//! TopOpt Solver as a library: topology optimization on all CPU cores, for any pipeline.
//!
//! Two ways in:
//! - [`run`]: one input ([`Input`]: a [`Spec`] described without a CAD app, or a plugin [`Job`])
//!   → optimized densities, a closed mesh and a summary, with progress [`Event`]s and an optional
//!   cancel flag. Nothing touches the disk unless [`RunOptions`] names output files.
//! - [`serve`]: the folder protocol of the SolidCraft plugin (`topopt_study::solver`): watch a
//!   folder, solve the jobs the plugin writes there, write the results back.
//!
//! The `topopt-solver` binary is a thin command line over both (`solve`, `watch`): [`cli`],
//! which `topopt-cuda-solver` reuses with its linear systems solved on a GPU
//! ([`RunOptions::solver`]).

pub mod cli;

use serde::Serialize;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};
use topopt_study::solver::{self, Alive, JobHeader, PROTOCOL, ResultFile};
use topopt_core::SolverFactory;
use topopt_study::{Density, Hooks, Job, Mesh, RunResult, Spec};

pub use topopt_study::spec;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub fn unix_now() -> u64 {
    SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0)
}

/// Worker threads in use (set `RAYON_NUM_THREADS` before the first run to change it).
pub fn threads() -> usize {
    topopt_core::par::threads()
}

// ---------------------------------------------------------------------------------------------
// One run.

/// What to solve.
#[derive(Clone, Debug)]
pub enum Input {
    /// A study described by regions and shapes; `base` resolves its relative file paths.
    Spec { spec: Spec, base: PathBuf },
    /// A job as the SolidCraft plugin writes it (already voxelized).
    Job(Box<Job>),
}

impl Input {
    /// Reads a JSON file: a plugin job when it has `problem` and `study`, else a spec.
    pub fn from_file(path: &Path) -> Result<Input, String> {
        let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let base = path.parent().map(Path::to_path_buf).unwrap_or_default();
        Input::from_json(&bytes, base).map_err(|e| format!("{}: {e}", path.display()))
    }

    pub fn from_json(bytes: &[u8], base: PathBuf) -> Result<Input, String> {
        let value: serde_json::Value = serde_json::from_slice(bytes).map_err(|e| format!("not JSON: {e}"))?;
        if value.get("problem").is_some() && value.get("study").is_some() {
            let job: Job = serde_json::from_value(value).map_err(|e| format!("not a TopOpt job: {e}"))?;
            if job.version != topopt_study::JOB_VERSION {
                return Err(format!("job version {} (this solver reads {})", job.version, topopt_study::JOB_VERSION));
            }
            Ok(Input::Job(Box::new(job)))
        } else {
            let spec: Spec = serde_json::from_value(value).map_err(|e| format!("not a TopOpt spec: {e}"))?;
            Ok(Input::Spec { spec, base })
        }
    }

    /// A name for reports and default file names.
    pub fn name(&self) -> String {
        match self {
            Input::Spec { spec, .. } => spec.name.clone().unwrap_or_else(|| "topopt".into()),
            Input::Job(job) => job.study.body.clone().unwrap_or_else(|| "topopt".into()),
        }
    }
}

/// Overrides and outputs of [`run`]. Every field is optional.
#[derive(Clone, Default)]
pub struct RunOptions {
    /// Write the mesh here (`.ply` or `.stl`, by extension).
    pub mesh_out: Option<PathBuf>,
    /// Write the densities here (JSON, [`Density`]).
    pub density_out: Option<PathBuf>,
    /// Write the summary here (JSON, [`RunResult`]).
    pub result_out: Option<PathBuf>,
    pub time_limit: Option<f64>,
    pub max_iterations: Option<usize>,
    /// Set to true from another thread to stop the run (it then returns `Err("cancelled")`).
    pub cancel: Option<Arc<AtomicBool>>,
    /// What solves the linear systems (the CPU multigrid-CG when `None`); `topopt-cuda`'s
    /// `factory` runs them on a GPU.
    pub solver: Option<SolverFactory>,
}

impl std::fmt::Debug for RunOptions {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RunOptions")
            .field("mesh_out", &self.mesh_out)
            .field("density_out", &self.density_out)
            .field("result_out", &self.result_out)
            .field("time_limit", &self.time_limit)
            .field("max_iterations", &self.max_iterations)
            .field("cancel", &self.cancel)
            .field("solver", &self.solver.as_ref().map(|_| "custom"))
            .finish()
    }
}

/// Progress of a run, in order: one `Started`, many `Progress`, then `Finished` or `Failed`.
/// Serialized (for `solve --progress json`) as `{"event": "progress", ...}`.
#[derive(Clone, Debug, PartialEq, Serialize)]
#[serde(tag = "event", rename_all = "snake_case")]
pub enum Event {
    Started { name: String, threads: usize, version: String },
    Progress { fraction: f32, message: String, elapsed: f64 },
    Finished { result: RunResult, files: Vec<String> },
    Failed { message: String },
}

/// What [`run`] produced.
pub struct RunOutput {
    pub result: RunResult,
    pub mesh: Mesh,
    pub density: Density,
    /// Files written (from [`RunOptions`]).
    pub files: Vec<PathBuf>,
}

struct EventHooks<'a> {
    start: Instant,
    cancel: Option<Arc<AtomicBool>>,
    solver: Option<SolverFactory>,
    on_event: &'a mut dyn FnMut(&Event),
}

impl Hooks for EventHooks<'_> {
    fn progress(&mut self, fraction: f32, message: &str) {
        (self.on_event)(&Event::Progress { fraction, message: message.into(), elapsed: self.start.elapsed().as_secs_f64() });
    }
    fn cancelled(&mut self) -> bool {
        self.cancel.as_ref().is_some_and(|c| c.load(Ordering::Relaxed))
    }
    fn elapsed(&mut self) -> f64 {
        self.start.elapsed().as_secs_f64()
    }
    fn solver(&self) -> Option<SolverFactory> {
        self.solver.clone()
    }
}

/// Writes through a temporary sibling and a rename, so a reader never sees half a file.
pub fn write_atomic(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(".tmp");
    std::fs::write(&tmp, bytes)?;
    std::fs::rename(&tmp, path)
}

/// The mesh as file bytes: STL for a `.stl` path, PLY otherwise.
pub fn mesh_bytes(mesh: &Mesh, path: &Path, comment: &str) -> Vec<u8> {
    match path.extension().and_then(|e| e.to_str()).map(str::to_ascii_lowercase).as_deref() {
        Some("stl") => mesh.to_stl(comment),
        _ => mesh.to_ply(comment),
    }
}

/// What [`check`] found: the problem the input becomes, without solving it.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct CheckReport {
    pub name: String,
    pub grid: [usize; 3],
    /// Element edge, mm.
    pub h: f64,
    /// Elements inside the domain (design + solid).
    pub elements: usize,
    pub design_elements: usize,
    pub solid_elements: usize,
    pub fixed_nodes: usize,
    pub load_cases: usize,
    pub loaded_nodes: usize,
    pub volfrac: f64,
    pub filter_radius: f64,
    /// Rough run time on this machine, s.
    pub estimate_seconds: f64,
    pub threads: usize,
    pub warnings: Vec<String>,
}

/// Turns the input into a job (voxelizes a spec) and checks it, without optimizing: a quick test
/// of a spec in a pipeline before the long run.
pub fn prepare(input: &Input, on_event: &mut dyn FnMut(&Event)) -> Result<Job, String> {
    match input {
        Input::Job(job) => Ok((**job).clone()),
        Input::Spec { spec, base } => {
            let mut hooks = EventHooks { start: Instant::now(), cancel: None, solver: None, on_event };
            let read = |p: &str| -> Result<Vec<u8>, String> {
                let path = base.join(p);
                std::fs::read(&path).map_err(|e| format!("{}: {e}", path.display()))
            };
            topopt_study::prepare_spec(spec, &read, &mut hooks)
        }
    }
}

/// What [`remesh`] made.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct RemeshReport {
    pub triangles: usize,
    pub vertices: usize,
    pub volume_mm3: f64,
    pub genus: i64,
    /// Vertices kept on the original surface, on its sharp edges and on its corners.
    pub pinned_vertices: usize,
    pub edge_vertices: usize,
    pub corner_vertices: usize,
    pub sharp_edges: bool,
    pub smooth_iterations: usize,
    /// Dead ends trimmed before meshing (with `trim`).
    pub trim: Option<topopt_study::trim::TrimReport>,
    pub warnings: Vec<String>,
}

/// Makes the mesh again from saved densities (a `density.json` of an earlier run of the same
/// input) with the input's mesh settings: tune the mesh without optimizing again. Dead weight is
/// trimmed first when the input's `mesh.trim_percent` (or `trim_percent`, which overrides it) is
/// above 0.
pub fn remesh(input: &Input, density: &Density, trim_percent: Option<f64>) -> Result<(Mesh, RemeshReport), String> {
    remesh_with(input, density, trim_percent, None)
}

/// [`remesh`] with the trimming's analyses solved by what `solver` makes.
pub fn remesh_with(input: &Input, density: &Density, trim_percent: Option<f64>, solver: Option<&SolverFactory>) -> Result<(Mesh, RemeshReport), String> {
    let job = prepare(input, &mut |_| {})?;
    let g = &job.problem.grid;
    if density.grid.nx != g.nx || density.grid.ny != g.ny || density.grid.nz != g.nz {
        return Err(format!(
            "the densities are on a {}×{}×{} grid, the input makes {}×{}×{}: they come from another study or settings",
            density.grid.nx, density.grid.ny, density.grid.nz, g.nx, g.ny, g.nz
        ));
    }
    let surface = job.surface.as_ref().map(|(p, t)| (p.as_slice(), t.as_slice()));
    let mut density = density.clone();
    let pct = trim_percent.unwrap_or(job.study.mesh.trim_percent);
    let trim_report = if pct > 0.0 {
        let mut settings = job.settings;
        settings.mg = topopt_core::mg::MgSettings::default();
        let max_loss = pct / 100.0;
        let t = topopt_study::trim::TrimSettings { max_loss, piece_loss: max_loss / 2.0, ..Default::default() };
        let x: Vec<f64> = density.values.iter().map(|&v| v as f64).collect();
        let (trimmed, rep) = topopt_study::trim::trim_dead_ends_with(&job.problem, &settings, &x, job.study.mesh.iso, &t, solver)?;
        topopt_study::trim::apply(&mut density, &trimmed);
        Some(rep)
    } else {
        None
    };
    let (mesh, r) = topopt_study::make_mesh(&job.study, &density, surface)?;
    let warnings = topopt_study::mesh_warnings(&job.study, &r, surface.is_some());
    let report = RemeshReport {
        triangles: r.mesh.triangles,
        vertices: r.mesh.vertices,
        volume_mm3: r.mesh.volume,
        genus: r.mesh.genus,
        pinned_vertices: r.pinned_vertices,
        edge_vertices: r.edge_vertices,
        corner_vertices: r.corner_vertices,
        sharp_edges: r.sharp_edges,
        smooth_iterations: r.smooth_iterations,
        trim: trim_report,
        warnings,
    };
    Ok((mesh, report))
}

pub fn check(input: &Input) -> Result<CheckReport, String> {
    use topopt_core::ElementKind;
    let job = prepare(input, &mut |_| {})?;
    job.problem.validate().map_err(|e| format!("the problem cannot be solved: {e}"))?;
    let p = &job.problem;
    let n = p.grid.nel();
    let count = |k: ElementKind| (0..n).filter(|&e| p.kind(e) == k).count();
    let (design, solid) = (count(ElementKind::Design), count(ElementKind::Solid));
    let elements = design + solid;
    Ok(CheckReport {
        name: input.name(),
        grid: [p.grid.nx, p.grid.ny, p.grid.nz],
        h: p.grid.h,
        elements,
        design_elements: design,
        solid_elements: solid,
        fixed_nodes: p.supports.len(),
        load_cases: p.load_cases.len(),
        loaded_nodes: p.load_cases.iter().map(|c| c.forces.len()).sum(),
        volfrac: job.settings.volfrac,
        filter_radius: job.settings.rmin * p.grid.h,
        estimate_seconds: topopt_study::status::estimate_seconds(&job.study, elements, threads()),
        threads: threads(),
        warnings: job.warnings.clone(),
    })
}

/// Solves one input: voxelize (for a spec), optimize, mesh; writes the requested files.
pub fn run(input: &Input, opts: &RunOptions, on_event: &mut dyn FnMut(&Event)) -> Result<RunOutput, String> {
    let name = input.name();
    on_event(&Event::Started { name: name.clone(), threads: threads(), version: VERSION.into() });
    let outcome = (|| -> Result<RunOutput, String> {
        let start = Instant::now();
        let mut job = prepare(input, &mut *on_event)?;
        let mut hooks = EventHooks { start, cancel: opts.cancel.clone(), solver: opts.solver.clone(), on_event: &mut *on_event };
        if let Some(t) = opts.time_limit {
            job.study.settings.time_limit = t;
        }
        if let Some(n) = opts.max_iterations {
            job.study.settings.max_iterations = n;
            job.settings.max_iter = n;
        }
        if let Some(p) = &opts.mesh_out {
            job.mesh_path = p.display().to_string();
        }
        let out = topopt_study::solve(&job, &mut hooks)?;
        let mut files = Vec::new();
        let mut write = |path: &Path, bytes: Vec<u8>| -> Result<(), String> {
            write_atomic(path, &bytes).map_err(|e| format!("{}: {e}", path.display()))?;
            files.push(path.to_path_buf());
            Ok(())
        };
        if let Some(p) = &opts.mesh_out {
            write(p, mesh_bytes(&out.mesh, p, &format!("TopOpt {name}")))?;
        }
        if let Some(p) = &opts.density_out {
            write(p, serde_json::to_vec(&out.density).map_err(|e| e.to_string())?)?;
        }
        if let Some(p) = &opts.result_out {
            write(p, serde_json::to_vec_pretty(&out.result).map_err(|e| e.to_string())?)?;
        }
        Ok(RunOutput { result: out.result, mesh: out.mesh, density: out.density, files })
    })();
    match &outcome {
        Ok(o) => on_event(&Event::Finished { result: o.result.clone(), files: o.files.iter().map(|f| f.display().to_string()).collect() }),
        Err(e) => on_event(&Event::Failed { message: e.clone() }),
    }
    outcome
}

// ---------------------------------------------------------------------------------------------
// The SolidCraft plugin's folder protocol.

/// The shared results folder.
#[derive(Clone)]
pub struct Folder {
    pub root: PathBuf,
}

impl Folder {
    pub fn path(&self, rel: &str) -> PathBuf {
        self.root.join(rel)
    }

    /// Ids of complete jobs without a result, oldest first.
    pub fn ready_jobs(&self) -> Vec<String> {
        let Ok(dir) = std::fs::read_dir(self.path(solver::JOBS)) else { return Vec::new() };
        let mut jobs: Vec<(SystemTime, String)> = Vec::new();
        for entry in dir.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();
            let Some(id) = name.strip_suffix(".json") else { continue };
            if self.path(&solver::result_file(id)).exists() {
                continue;
            }
            let Ok(text) = std::fs::read(entry.path()) else { continue };
            let Ok(h) = serde_json::from_slice::<JobHeader>(&text) else { continue };
            let total: u64 = (0..h.parts).map(|k| std::fs::metadata(self.path(&solver::job_part(id, k))).map(|m| m.len()).unwrap_or(0)).sum();
            if h.protocol == PROTOCOL && total == h.bytes as u64 {
                jobs.push((entry.metadata().and_then(|m| m.modified()).unwrap_or(UNIX_EPOCH), id.to_string()));
            }
        }
        jobs.sort();
        jobs.into_iter().map(|(_, id)| id).collect()
    }

    pub fn load_job(&self, id: &str) -> Result<Job, String> {
        let header: JobHeader = serde_json::from_slice(&std::fs::read(self.path(&solver::job_header(id))).map_err(|e| e.to_string())?)
            .map_err(|e| format!("job header: {e}"))?;
        let mut bytes = Vec::with_capacity(header.bytes);
        for k in 0..header.parts {
            bytes.extend(std::fs::read(self.path(&solver::job_part(id, k))).map_err(|e| format!("job part {k}: {e}"))?);
        }
        match Input::from_json(&bytes, self.root.clone())? {
            Input::Job(job) => Ok(*job),
            Input::Spec { .. } => Err("the job file holds a spec, not a plugin job".into()),
        }
    }

    fn remove_job_files(&self, id: &str) {
        let _ = std::fs::remove_file(self.path(&solver::job_header(id)));
        let _ = std::fs::remove_file(self.path(&solver::job_cancel(id)));
        for k in 0.. {
            if std::fs::remove_file(self.path(&solver::job_part(id, k))).is_err() {
                break;
            }
        }
    }
}

/// The heartbeat, rewritten every second by a background thread.
pub struct Heartbeat {
    pub state: Arc<Mutex<Alive>>,
}

impl Heartbeat {
    pub fn start(folder: &Folder) -> Heartbeat {
        let state = Arc::new(Mutex::new(Alive {
            protocol: PROTOCOL,
            version: VERSION.into(),
            time: unix_now(),
            threads: threads(),
            job: None,
            progress: 0.0,
            message: "idle".into(),
        }));
        let (s, path) = (state.clone(), folder.path(solver::ALIVE));
        std::thread::spawn(move || loop {
            let snapshot = {
                let mut a = s.lock().unwrap();
                a.time = unix_now();
                a.clone()
            };
            let _ = write_atomic(&path, &serde_json::to_vec(&snapshot).unwrap_or_default());
            std::thread::sleep(Duration::from_secs(1));
        });
        Heartbeat { state }
    }

    pub fn set(&self, job: Option<&str>, progress: f32, message: &str) {
        let mut a = self.state.lock().unwrap();
        a.job = job.map(str::to_owned);
        a.progress = progress;
        a.message = message.into();
    }
}

/// Solves one plugin job from the folder and writes its mesh, densities and result file (the
/// result file last). `on_event` sees the run's progress; `solver` as in [`RunOptions::solver`].
pub fn run_job(folder: &Folder, id: &str, beat: &Heartbeat, solver: Option<&SolverFactory>, on_event: &mut dyn FnMut(&Event)) -> ResultFile {
    let mut file = ResultFile { protocol: PROTOCOL, id: id.into(), threads: threads(), ..Default::default() };
    let cancel = Arc::new(AtomicBool::new(false));
    let cancel_path = folder.path(&solver::job_cancel(id));
    let outcome = folder.load_job(id).and_then(|job| {
        let mesh_path = folder.path(&job.mesh_path);
        let input = Input::Job(Box::new(job));
        let opts = RunOptions { mesh_out: Some(mesh_path), cancel: Some(cancel.clone()), solver: solver.cloned(), ..Default::default() };
        run(&input, &opts, &mut |e: &Event| {
            if let Event::Progress { fraction, message, .. } = e {
                beat.set(Some(id), *fraction, message);
                if cancel_path.exists() {
                    cancel.store(true, Ordering::Relaxed);
                }
            }
            on_event(e);
        })
    });
    match outcome {
        Ok(out) => {
            let density = serde_json::to_vec(&out.density).unwrap_or_default();
            let parts = solver::split(&density);
            let mut ok = true;
            for (k, part) in parts.iter().enumerate() {
                if let Err(e) = write_atomic(&folder.path(&solver::density_part(id, k)), part) {
                    file.error = Some(format!("writing densities: {e}"));
                    ok = false;
                    break;
                }
            }
            if ok {
                file.density_parts = parts.len();
                file.density_bytes = density.len();
                // The plugin expects the mesh path relative to its folder, as it sent it.
                let mut result = out.result;
                if let Ok(rel) = Path::new(&result.stl).strip_prefix(&folder.root) {
                    result.stl = rel.to_string_lossy().replace('\\', "/");
                }
                file.result = Some(result);
            }
        }
        Err(e) => file.error = Some(e),
    }
    let _ = write_atomic(&folder.path(&solver::result_file(id)), &serde_json::to_vec(&file).unwrap_or_default());
    folder.remove_job_files(id);
    file
}

/// Watches the folder and solves jobs as they come; with `once`, solves what is there and returns
/// the number solved.
pub fn serve(folder: &Folder, once: bool, on_event: &mut dyn FnMut(&str, &Event)) -> usize {
    serve_with(folder, once, None, on_event)
}

/// [`serve`] with the linear systems solved by what `solver` makes.
pub fn serve_with(folder: &Folder, once: bool, solver: Option<&SolverFactory>, on_event: &mut dyn FnMut(&str, &Event)) -> usize {
    let beat = Heartbeat::start(folder);
    let mut solved = 0;
    loop {
        for id in folder.ready_jobs() {
            run_job(folder, &id, &beat, solver, &mut |e| on_event(&id, e));
            solved += 1;
            beat.set(None, 0.0, "idle");
        }
        if once {
            return solved;
        }
        std::thread::sleep(Duration::from_millis(500));
    }
}
