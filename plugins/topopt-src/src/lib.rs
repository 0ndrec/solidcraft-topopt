//! TopOpt: the SolidCraft plugin.
//!
//! The study lives in the design (`storage.set` document); the panel shows and edits it; dialogs
//! (study, support, load) take faces and bodies picked in the viewport. `optimize` is a long job:
//! ray voxelization through `view.pick`, topology optimization, mesh, STL into the `results`
//! folder. `create-body` inserts that STL as a mesh body (a separate run: a file written by a run
//! is only published when the run ends) and can convert it to a solid.

wit_bindgen::generate!({ path: "../../wit", world: "panel-plugin" });

mod panel;

use craft::plugin::types::{FieldFailure, PluginError};
use craft::plugin::{clock, commands, control, files, forms, storage, ui};
use serde_json::{Value, json};
use topopt_study::solver::{self, Alive, JobHeader, Pending, ResultFile, SolverView};
use topopt_study::{BodyInfo, Caster, Condition, Density, FaceRef, Hit, Hooks, Material, Study};

struct Plugin;

pub(crate) fn invalid(field: &str, reason: impl ToString) -> PluginError {
    PluginError::InvalidInput(FieldFailure { field: field.into(), reason: reason.to_string() })
}

/// Runs an app command and parses its reply.
pub(crate) fn run(command: &str, args: Value) -> Result<Value, PluginError> {
    let reply = commands::run(command, &args.to_string())?;
    serde_json::from_str(&reply).map_err(|e| invalid(command, e))
}

fn load() -> Result<Study, PluginError> {
    let reply = run("storage.get", json!({}))?;
    let doc = &reply["document"];
    if doc.is_null() {
        return Ok(Study::default());
    }
    serde_json::from_value(doc.clone()).map_err(|e| invalid("study", format!("the stored study does not read: {e}")))
}

fn save(study: &Study) -> Result<(), PluginError> {
    run("storage.set", json!({"document": study}))?;
    Ok(())
}

/// Body names from a selection field: `[{"type": "body", "name": "Body1"}]`.
fn body_names(v: &Value) -> Vec<String> {
    v.as_array()
        .map(|a| a.iter().filter_map(|x| x["name"].as_str().or_else(|| x.as_str()).map(str::to_owned)).collect())
        .unwrap_or_default()
}

/// Faces from a selection field.
fn faces(v: &Value) -> Result<Vec<FaceRef>, PluginError> {
    let list = v.as_array().ok_or_else(|| invalid("faces", "pick at least one face"))?;
    list.iter()
        .map(|f| {
            let body = f["body"].as_str().ok_or_else(|| invalid("faces", "a pick without its body"))?;
            let index = f["index"].as_i64().ok_or_else(|| invalid("faces", "a pick without a face index"))?;
            let point = f["point"].as_array().and_then(|p| Some([p.first()?.as_f64()?, p.get(1)?.as_f64()?, p.get(2)?.as_f64()?]));
            Ok(FaceRef { body: body.into(), index, point, label: None })
        })
        .collect()
}

fn face_refs_json(faces: &[FaceRef]) -> Value {
    Value::Array(
        faces
            .iter()
            .map(|f| match f.point {
                Some(p) => json!({"type": "face", "body": f.body, "index": f.index, "point": p}),
                None => json!({"type": "face", "body": f.body, "index": f.index}),
            })
            .collect(),
    )
}

fn num(params: &Value, key: &str, default: f64) -> f64 {
    params[key].as_f64().unwrap_or(default)
}

impl exports::craft::plugin::entry::Guest for Plugin {
    fn execute(entry: String, params: String) -> Result<String, PluginError> {
        let params: Value = serde_json::from_str(&params).map_err(|e| invalid("params", e))?;
        // A panel button runs its entry with the click as params: open that entry's form.
        let pressed = params["panel"].as_str() == Some(panel::PANEL) && params["event"].as_str() == Some("press");
        let reply = match entry.as_str() {
            "design-body" if pressed => open_study()?,
            "add-support" | "add-load" | "delete-item" if pressed => {
                forms::open(&entry, "{}")?;
                json!({"opened": entry})
            }
            "design-body" => edit_study(&params)?,
            "add-support" => add_condition(&params, Kind::Support)?,
            "add-load" => add_condition(&params, Kind::Load)?,
            "keep-face" => add_condition(&params, Kind::Keep)?,
            "delete-item" => remove(&params)?,
            "clear-all" => {
                let mut s = load()?;
                s.conditions.clear();
                s.result = None;
                save(&s)?;
                json!({"conditions": 0})
            }
            "panel" => match params["event"].as_str() {
                Some("change") => set_option(&params)?,
                Some("edit") => edit_item(&params)?,
                // Run from the menu or the ribbon: show the panel (again).
                _ => {
                    ui::show_panel(panel::PANEL, true)?;
                    json!({"shown": panel::PANEL})
                }
            },
            "optimize" => optimize()?,
            "remesh" => remesh()?,
            "create-body" => create_body()?,
            "refresh" => json!({}),
            _ => return Err(invalid("entry", format!("no entry `{entry}`"))),
        };
        // Redraw the panel from the stored study; a panel failure never fails the entry.
        if let Ok(study) = load() {
            let info = study.body.as_deref().and_then(|b| body_info(b).ok());
            let view = solver_view(&study);
            let _ = panel::show(&study, info.as_ref(), &view);
        }
        Ok(reply.to_string())
    }
}

fn edit_study(p: &Value) -> Result<Value, PluginError> {
    let mut s = load()?;
    let body = body_names(&p["body"]).into_iter().next().ok_or_else(|| invalid("body", "pick the design body"))?;
    s.body = Some(body.clone());
    s.keep_out = body_names(&p["keep_out"]).into_iter().filter(|b| *b != body).collect();
    if let Some(m) = p["material"].as_str() {
        if let Some(mat) = Material::preset(m) {
            s.material = mat;
        } else {
            s.material.preset = "custom".into();
        }
    }
    let st = &mut s.settings;
    st.volfrac = match p["volume_percent"].as_f64() {
        Some(pc) => (pc / 100.0).clamp(0.02, 0.95),
        None => num(p, "volfrac", st.volfrac).clamp(0.02, 0.95),
    };
    st.element_size = num(p, "element_size", st.element_size).max(0.0);
    st.max_elements = num(p, "max_elements", st.max_elements as f64).clamp(1000.0, 400_000.0) as usize;
    st.min_member = num(p, "min_member", st.min_member).max(0.0);
    // A study on another body: its faces no longer apply.
    let before = s.conditions.len();
    s.conditions.retain(|c| c.faces().iter().all(|f| f.body == body));
    s.result = None;
    save(&s)?;
    Ok(json!({"body": body, "dropped_items": before - s.conditions.len()}))
}

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Support,
    Load,
    Keep,
}

fn add_condition(p: &Value, kind: Kind) -> Result<Value, PluginError> {
    let mut s = load()?;
    let replace = p["replace"].as_u64().filter(|&r| r > 0).map(|r| r as u32);
    if p["delete"].as_bool() == Some(true) {
        let id = replace.ok_or_else(|| invalid("replace", "give the item number to delete"))?;
        if !s.remove(id) {
            return Err(invalid("replace", format!("there is no item {id}")));
        }
        s.result = None;
        save(&s)?;
        return Ok(json!({"removed": id}));
    }
    let mut fs = faces(&p["faces"])?;
    if let Some(b) = &s.body {
        if let Some(f) = fs.iter().find(|f| &f.body != b) {
            return Err(invalid(
                "faces",
                format!("that face belongs to {}, but the design body is {b}: pick faces of {b}, or change the design body first", f.body),
            ));
        }
    } else {
        // No study yet: the faces' body becomes the design body.
        s.body = fs.first().map(|f| f.body.clone());
    }
    label_faces(&mut fs);
    let c = if kind == Kind::Keep {
        Condition::Keep { id: 0, faces: fs, thickness: num(p, "thickness", 0.0).max(0.0) }
    } else if kind == Kind::Support {
        let dirs = [p["fix_x"].as_bool().unwrap_or(true), p["fix_y"].as_bool().unwrap_or(true), p["fix_z"].as_bool().unwrap_or(true)];
        if !dirs.iter().any(|&d| d) {
            return Err(invalid("fix_x", "tick at least one of Hold in X, Y, Z"));
        }
        Condition::Support { id: 0, faces: fs, dirs }
    } else {
        let mode = p["mode"].as_str().unwrap_or("vector").to_owned();
        let force = [num(p, "fx", 0.0), num(p, "fy", 0.0), num(p, "fz", -1000.0)];
        let magnitude = num(p, "magnitude", 1000.0);
        if (mode == "vector" && force.iter().all(|f| *f == 0.0)) || (mode == "normal" && magnitude == 0.0) {
            return Err(invalid("fz", "the force is zero: enter Fx, Fy, Fz (or a Magnitude for 'Pushing into the faces')"));
        }
        Condition::Load { id: 0, faces: fs, mode, force, magnitude, case: num(p, "case", 1.0).clamp(1.0, 8.0) as u32 }
    };
    let id = s.add(c, replace);
    s.result = None;
    save(&s)?;
    Ok(json!({"item": id}))
}

fn remove(p: &Value) -> Result<Value, PluginError> {
    let mut s = load()?;
    let id = p["id"].as_u64().filter(|&i| i > 0).map(|i| i as u32).ok_or_else(|| invalid("id", "give the item number from the # column"))?;
    if !s.remove(id) {
        return Err(invalid("id", format!("there is no item {id} in the list")));
    }
    s.result = None;
    save(&s)?;
    Ok(json!({"removed": id}))
}

fn set_option(p: &Value) -> Result<Value, PluginError> {
    let path = p["path"].as_str().ok_or_else(|| invalid("path", "a panel field change"))?;
    let mut s = load()?;
    s.set_option(path, &p["value"]).map_err(|e| invalid("value", e))?;
    save(&s)?;
    Ok(json!({"set": path}))
}

fn open_study() -> Result<Value, PluginError> {
    let s = load()?;
    let st = &s.settings;
    let mut v = json!({
        "material": s.material.preset,
        "volume_percent": (st.volfrac * 1000.0).round() / 10.0,
        "element_size": st.element_size,
        "max_elements": st.max_elements,
        "min_member": st.min_member,
    });
    if let Some(b) = &s.body {
        v["body"] = json!([{"type": "body", "name": b}]);
    }
    if !s.keep_out.is_empty() {
        v["keep_out"] = Value::Array(s.keep_out.iter().map(|b| json!({"type": "body", "name": b})).collect());
    }
    forms::open("design-body", &v.to_string())?;
    Ok(json!({"opened": "design-body"}))
}

/// Double-click on a row: its dialog, filled in, replacing the item on Apply.
fn edit_item(p: &Value) -> Result<Value, PluginError> {
    let id: u32 = p["row"].as_str().and_then(|r| r.parse().ok()).ok_or_else(|| invalid("row", "no row"))?;
    let s = load()?;
    let c = s.conditions.iter().find(|c| c.id() == id).ok_or_else(|| invalid("row", format!("no item {id}")))?;
    match c {
        Condition::Support { faces, dirs, .. } => forms::open(
            "add-support",
            &json!({"faces": face_refs_json(faces), "fix_x": dirs[0], "fix_y": dirs[1], "fix_z": dirs[2], "replace": id}).to_string(),
        )?,
        Condition::Load { faces, mode, force, magnitude, case, .. } => forms::open(
            "add-load",
            &json!({"faces": face_refs_json(faces), "mode": mode, "fx": force[0], "fy": force[1], "fz": force[2],
                    "magnitude": magnitude, "case": case, "replace": id})
            .to_string(),
        )?,
        Condition::Keep { faces, thickness, .. } => {
            forms::open("keep-face", &json!({"faces": face_refs_json(faces), "thickness": thickness, "replace": id}).to_string())?
        }
    }
    Ok(json!({"editing": id}))
}

// ---------------------------------------------------------------------------------------------
// The job.

/// `view.pick` as a ray caster for one body.
struct PickCaster {
    body: String,
}

fn vec3(v: &Value) -> Option<[f64; 3]> {
    Some([v.get(0)?.as_f64()?, v.get(1)?.as_f64()?, v.get(2)?.as_f64()?])
}

impl Caster for PickCaster {
    fn cast(&mut self, rays: &[([f64; 3], [f64; 3])]) -> Result<Vec<Option<Hit>>, String> {
        let list: Vec<Value> = rays.iter().map(|(o, d)| json!({"origin": o, "dir": d})).collect();
        let reply = run("view.pick", json!({"rays": list, "all": true})).map_err(|e| format!("view.pick: {e:?}"))?;
        let hits = reply["hits"].as_array().ok_or("view.pick: no hits list")?;
        Ok(hits
            .iter()
            .map(|h| {
                let candidates: Vec<&Value> = match h {
                    Value::Array(a) => a.iter().collect(),
                    Value::Object(_) => vec![h],
                    _ => vec![],
                };
                let hit = candidates.into_iter().find(|c| c["body"].as_str() == Some(self.body.as_str()))?;
                Some(Hit { point: vec3(&hit["point"])?, normal: vec3(&hit["normal"])?, face: hit["face"].as_i64()? })
            })
            .collect())
    }
}

struct JobHooks {
    start: u64,
}

impl Hooks for JobHooks {
    fn progress(&mut self, fraction: f32, message: &str) {
        control::progress(fraction.clamp(0.0, 1.0), message);
    }
    fn cancelled(&mut self) -> bool {
        control::cancelled()
    }
    fn elapsed(&mut self) -> f64 {
        clock::now().saturating_sub(self.start) as f64
    }
}

fn body_info(body: &str) -> Result<BodyInfo, PluginError> {
    let m = run("inspect.measure", json!({"bodies": [body]}))?;
    let b = m["bodies"].get(0).ok_or_else(|| invalid("body", format!("no body named {body}")))?;
    if b["closed"].as_bool() == Some(false) {
        return Err(invalid("body", format!("{body} is not a closed solid")));
    }
    let lo = vec3(&b["bbox"]["min"]).ok_or_else(|| invalid("body", "no bounding box"))?;
    let hi = vec3(&b["bbox"]["max"]).ok_or_else(|| invalid("body", "no bounding box"))?;
    Ok(BodyInfo { lo, hi, volume: b["volume_mm3"].as_f64().unwrap_or(0.0) })
}

/// The body's surface: positions, triangles, and (when the study keeps surfaces) the face each
/// triangle lies on. Coarser when the body is too detailed for `model.mesh`.
struct BodySurface {
    positions: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    faces: Vec<i64>,
}

impl BodySurface {
    fn as_ref(&self) -> topopt_study::SurfaceRef<'_> {
        topopt_study::SurfaceRef::new(&self.positions, &self.triangles).with_faces(&self.faces)
    }
}

fn surface(body: &str, study: &Study) -> Option<BodySurface> {
    for tol in [0.05, 0.2, 0.5, 2.0] {
        if let Ok(m) = run("model.mesh", json!({"body": body, "tolerance": tol})) {
            let positions: Vec<[f64; 3]> = m["positions"].as_array()?.iter().filter_map(vec3).collect();
            let triangles: Vec<[u32; 3]> = m["triangles"]
                .as_array()?
                .iter()
                .filter_map(|t| Some([t.get(0)?.as_u64()? as u32, t.get(1)?.as_u64()? as u32, t.get(2)?.as_u64()? as u32]))
                .collect();
            if triangles.is_empty() {
                continue;
            }
            // Kept surfaces take these triangles as they are: which face is each on? A ray from
            // just outside (beyond the chord error) says.
            let faces = if study.conditions.iter().any(|c| matches!(c, Condition::Keep { .. })) {
                let mut caster = PickCaster { body: body.to_string() };
                topopt_study::label_faces(&mut caster, &positions, &triangles, 2.0 * tol + 0.01).unwrap_or_default()
            } else {
                Vec::new()
            };
            return Some(BodySurface { positions, triangles, faces });
        }
    }
    None
}

fn mesh_path(body: &str) -> String {
    let name: String = body.chars().map(|c| if c.is_ascii_alphanumeric() || c == '-' || c == '_' { c } else { '_' }).collect();
    format!("topopt/{}.ply", if name.is_empty() { "result".into() } else { name })
}

/// One plugin file write may hold at most 1 MiB, so the mesh goes out as indexed binary PLY.
fn write_mesh(path: &str, mesh: &topopt_study::Mesh) -> Result<(), PluginError> {
    let bytes = mesh.to_ply("SolidCraft TopOpt");
    if bytes.len() > 1 << 20 {
        return Err(invalid(
            "mesh",
            format!("the mesh is {} KB, over the 1 MB a plugin may write: lower Mesh › Triangles and Re-mesh", bytes.len() / 1024),
        ));
    }
    files::write_file("results", path, &bytes)
}

/// Readable names of picked faces, from `model.faces`.
fn label_faces(faces: &mut [FaceRef]) {
    let Some(body) = faces.first().map(|f| f.body.clone()) else { return };
    let Ok(m) = run("model.faces", json!({"body": body})) else { return };
    let list = m["faces"].as_array().cloned().unwrap_or_default();
    for f in faces.iter_mut() {
        if let Some(info) = list.iter().find(|x| x["index"].as_i64() == Some(f.index)) {
            f.label = Some(face_label(info));
        }
    }
}

fn face_label(f: &Value) -> String {
    if let Some(r) = f["cylinder"]["radius"].as_f64() {
        return format!("Ø{} cylinder", topopt_study::study::fmt_n((2.0 * r * 100.0).round() / 100.0));
    }
    if let Some(n) = vec3(&f["plane_normal"]) {
        for (axis, name) in ["X", "Y", "Z"].iter().enumerate() {
            if n[axis].abs() > 0.999 {
                return format!("{}{name} face", if n[axis] < 0.0 { "−" } else { "+" });
            }
        }
        return "flat face".into();
    }
    "curved face".into()
}

const DENSITY_BLOB: &str = "density.json";

fn save_density(d: &Density) -> Result<(), PluginError> {
    let bytes = serde_json::to_vec(d).map_err(|e| invalid("density", e))?;
    for (k, chunk) in bytes.chunks(1 << 20).enumerate() {
        storage::put_blob(DENSITY_BLOB, chunk, true, k > 0)?;
    }
    Ok(())
}

fn load_density() -> Result<Density, PluginError> {
    let mut bytes = Vec::new();
    loop {
        let part = storage::get_blob(None, DENSITY_BLOB, bytes.len() as u64, 1 << 20)?;
        if part.is_empty() {
            break;
        }
        bytes.extend_from_slice(&part);
    }
    if bytes.is_empty() {
        return Err(invalid("density", "no optimized result to mesh yet: run Optimize"));
    }
    serde_json::from_slice(&bytes).map_err(|e| invalid("density", e))
}

fn optimize() -> Result<Value, PluginError> {
    let mut study = load()?;
    let problems = study.problems();
    if !problems.is_empty() {
        return Err(invalid("study", format!("not ready yet: {}", problems.join("; "))));
    }
    let body = study.body.clone().unwrap();
    // A job still with the solver is replaced by this one.
    cancel_pending(&study);
    study.pending = None;
    study.last_error = None;
    let mut hooks = JobHooks { start: clock::now() };
    control::progress(0.0, "reading the design body");
    let info = body_info(&body)?;
    let surf = if study.pins_faces() { surface(&body, &study) } else { None };
    let mut design = PickCaster { body: body.clone() };
    let mut outs: Vec<PickCaster> = study.keep_out.iter().map(|b| PickCaster { body: b.clone() }).collect();
    let mut out_refs: Vec<&mut dyn Caster> = outs.iter_mut().map(|c| c as &mut dyn Caster).collect();
    let surface_ref = surf.as_ref().map(|s| s.as_ref());
    let mut job = topopt_study::prepare(&study, info, &mut design, &mut out_refs, surface_ref, &mut hooks).map_err(|e| invalid("optimize", e))?;
    let path = mesh_path(&body);
    let alive = if study.settings.use_solver { read_alive() } else { None };
    if let Some(alive) = alive {
        // TopOpt Solver runs: hand the job over (it starts when this run ends and its files publish).
        let now = clock::now();
        job.id = format!("j{now}-{}", study.next_id);
        job.mesh_path = path;
        let bytes = serde_json::to_vec(&job).map_err(|e| invalid("job", e))?;
        let parts = solver::split(&bytes);
        for (k, part) in parts.iter().enumerate() {
            files::write_file("results", &solver::job_part(&job.id, k), part)?;
        }
        let header = JobHeader { protocol: solver::PROTOCOL, id: job.id.clone(), parts: parts.len(), bytes: bytes.len() };
        files::write_file("results", &solver::job_header(&job.id), &serde_json::to_vec(&header).unwrap_or_default())?;
        study.pending = Some(Pending { id: job.id.clone(), sent: now });
        study.result = None;
        save(&study)?;
        control::progress(1.0, "sent to TopOpt Solver");
        return Ok(json!({"sent_to_solver": job.id, "solver_threads": alive.threads, "job_bytes": bytes.len()}));
    }
    let outcome = topopt_study::solve(&job, &mut hooks).map_err(|e| {
        if e == "cancelled" { PluginError::Cancelled } else { invalid("optimize", e) }
    })?;
    save_density(&outcome.density)?;
    let mut result = outcome.result;
    if outcome.mesh.triangles.is_empty() {
        result.stl = String::new();
    } else {
        write_mesh(&path, &outcome.mesh)?;
        result.stl = path;
    }
    result.inserted = false;
    if surf.is_none() && study.pins_faces() {
        result.warnings.push("the body's surface was not available: faces are not kept exactly".into());
    }
    study.result = Some(result.clone());
    save(&study)?;
    Ok(serde_json::to_value(&result).unwrap_or_default())
}

fn remesh() -> Result<Value, PluginError> {
    let mut study = load()?;
    take_solver_result(&mut study)?;
    let body = study.body.clone().ok_or_else(|| invalid("study", "no design body"))?;
    control::progress(0.0, "reading the optimized densities");
    let density = load_density()?;
    let info = body_info(&body)?;
    let surf = if study.pins_faces() { surface(&body, &study) } else { None };
    control::progress(0.3, "making the mesh");
    let surface_ref = surf.as_ref().map(|s| s.as_ref());
    let (mesh, full) = topopt_study::make_mesh(&study, &density, surface_ref).map_err(|e| invalid("mesh", e))?;
    let report = &full.mesh;
    let path = mesh_path(&body);
    write_mesh(&path, &mesh)?;
    let mut result = study.result.clone().unwrap_or_default();
    result.stl = path;
    result.inserted = false;
    result.warnings.retain(|w| !w.starts_with("the sharp edges"));
    result.warnings.extend(topopt_study::mesh_warnings(&study, &full, surface_ref.is_some()));
    result.triangles = report.triangles;
    result.genus = report.genus;
    result.mesh_volume_mm3 = report.volume;
    result.mass_g = report.volume * study.material.density / 1000.0;
    if info.volume > 0.0 {
        result.volume_fraction = report.volume / info.volume;
    }
    study.result = Some(result.clone());
    save(&study)?;
    control::progress(1.0, "done");
    Ok(serde_json::to_value(&result).unwrap_or_default())
}

fn create_body() -> Result<Value, PluginError> {
    let study = load()?;
    let mut study = study;
    take_solver_result(&mut study)?;
    if study.pending.is_some() {
        let view = solver_view(&study);
        return Err(invalid(
            "result",
            match view.alive {
                Some(a) => format!("TopOpt Solver is still working ({:.0} %): wait, click Refresh, then Create body", a.progress * 100.0),
                None => "the job is waiting for TopOpt Solver, which is not running: start it, or click Optimize to run here".into(),
            },
        ));
    }
    let result = study.result.clone().ok_or_else(|| invalid("result", "there is no result yet: click 3. Optimize first and wait until it finishes"))?;
    if result.stl.is_empty() {
        return Err(invalid("result", "the last run made no mesh: click Re-mesh"));
    }
    let name = format!("TopOpt {}", study.body.as_deref().unwrap_or("result"));
    let inserted = run("file.insert_mesh", json!({"path": format!("grant:results/{}", result.stl), "units": "mm", "name": name}))?;
    let body = inserted["bodies"][0].as_str().unwrap_or_default().to_owned();
    let mut converted = Value::Null;
    let mut convert_error = None;
    if study.mesh.convert_solid && !body.is_empty() {
        // The mesh body stays when the CAD kernel cannot rebuild it as a solid: say why, and how
        // to try again, rather than losing the insert too.
        match run("mesh.convert", json!({"body": body, "name": format!("{name} solid")})) {
            Ok(v) => converted = v,
            Err(e) => convert_error = Some(format!("{e:?}")),
        }
    }
    let measured = run("inspect.measure", json!({"bodies": [body]}))?;
    if let Some(r) = study.result.as_mut() {
        r.inserted = true;
        r.warnings.retain(|w| !w.starts_with("the body was inserted as a mesh"));
        if let Some(e) = &convert_error {
            r.warnings.push(format!(
                "the body was inserted as a mesh: SolidCraft could not make it a solid ({e}); Re-mesh with other settings (smoothing, triangles) and Create body again, or use the mesh body"
            ));
        }
    }
    save(&study)?;
    Ok(json!({"body": body, "converted": !converted.is_null(), "convert_error": convert_error, "measure": measured["bodies"][0]}))
}

// ---------------------------------------------------------------------------------------------
// TopOpt Solver (see topopt_study::solver for the files).

/// A file in the results folder, or None when it does not exist (or cannot be read).
fn read_result_file(path: &str) -> Option<Vec<u8>> {
    files::read_file("results", path).ok()
}

/// The solver's heartbeat, when it is fresh.
fn read_alive() -> Option<Alive> {
    let a: Alive = serde_json::from_slice(&read_result_file(solver::ALIVE)?).ok()?;
    a.online(clock::now()).then_some(a)
}

/// What the panel shows about the solver (reads only).
fn solver_view(study: &Study) -> SolverView {
    let ready = study
        .pending
        .as_ref()
        .and_then(|p| read_result_file(&solver::result_file(&p.id)))
        .and_then(|b| serde_json::from_slice::<ResultFile>(&b).ok());
    SolverView { alive: read_alive(), ready }
}

/// Asks the solver to drop the pending job (best effort).
fn cancel_pending(study: &Study) {
    if let Some(p) = &study.pending {
        let _ = files::write_file("results", &solver::job_cancel(&p.id), b"cancel");
    }
}

/// Moves a finished solver job into the study: its result, and its densities for Re-mesh.
fn take_solver_result(study: &mut Study) -> Result<(), PluginError> {
    let Some(p) = study.pending.clone() else { return Ok(()) };
    let Some(bytes) = read_result_file(&solver::result_file(&p.id)) else { return Ok(()) };
    let file: ResultFile = serde_json::from_slice(&bytes).map_err(|e| invalid("solver", format!("the solver's result does not read: {e}")))?;
    study.pending = None;
    match (file.result, file.error) {
        (_, Some(e)) => {
            study.last_error = Some(e);
            study.result = None;
        }
        (Some(mut r), None) => {
            let mut density = Vec::with_capacity(file.density_bytes);
            for k in 0..file.density_parts {
                density.extend(read_result_file(&solver::density_part(&p.id, k)).ok_or_else(|| invalid("solver", "a density file of the solver is missing"))?);
            }
            if density.len() == file.density_bytes && !density.is_empty() {
                for (k, chunk) in density.chunks(1 << 20).enumerate() {
                    storage::put_blob(DENSITY_BLOB, chunk, true, k > 0)?;
                }
            }
            r.inserted = false;
            // The solver names what solved the systems (a GPU with TopOpt CUDA Solver); older
            // solvers only report their threads.
            let on = if r.solver.is_empty() { format!("{} threads", file.threads) } else { r.solver.clone() };
            r.warnings.push(format!("solved by TopOpt Solver on {on}"));
            study.result = Some(r);
            study.last_error = None;
        }
        (None, None) => study.last_error = Some("the solver returned no result".into()),
    }
    save(study)
}

impl exports::craft::plugin::ui_provider::Guest for Plugin {
    fn describe() -> craft::plugin::ui::Declarations {
        panel::declarations()
    }
    fn on_event(_event: craft::plugin::ui::UiEvent) -> Result<(), PluginError> {
        Ok(())
    }
    fn on_enter(_workspace: String) -> Result<(), PluginError> {
        Ok(())
    }
    fn on_exit(_workspace: String) -> Result<(), PluginError> {
        Ok(())
    }
}

export!(Plugin);
