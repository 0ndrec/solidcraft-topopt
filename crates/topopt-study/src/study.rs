//! The study the plugin keeps in the design (`storage.set` document): what to optimize, how,
//! and the last result. Every user-facing setting lives here, with its default and range.

use serde::{Deserialize, Serialize};
use serde_json::Value;

/// A face as SolidCraft's selection reference gives it:
/// `{"type": "face", "body": "Block", "index": 3, "point": [x, y, z]}`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FaceRef {
    pub body: String,
    pub index: i64,
    #[serde(default)]
    pub point: Option<[f64; 3]>,
    /// A readable name for the panel ("−X face", "Ø10 cylinder"); not part of the geometry.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
}

impl FaceRef {
    pub fn display(&self) -> String {
        self.label.clone().unwrap_or_else(|| format!("face {}", self.index))
    }
}

/// The most faces a mesh may have: the plugin writes it as one binary PLY of at most 1 MiB
/// (the host's limit on a file write), about 19 bytes per triangle.
pub const MAX_TRIANGLES: usize = 50_000;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Material {
    /// Preset id from [`MATERIALS`], or `custom`.
    pub preset: String,
    /// Young's modulus, MPa.
    pub e: f64,
    pub nu: f64,
    /// g/cm³, for the reported mass.
    pub density: f64,
}

/// Presets: id, label, E (MPa), ν, density (g/cm³).
pub const MATERIALS: [(&str, &str, f64, f64, f64); 9] = [
    ("aluminium-6061", "Aluminium 6061-T6", 68_900.0, 0.33, 2.70),
    ("aluminium-7075", "Aluminium 7075-T6", 71_700.0, 0.33, 2.81),
    ("steel", "Steel (structural)", 200_000.0, 0.29, 7.85),
    ("stainless-316", "Stainless steel 316", 193_000.0, 0.30, 8.00),
    ("titanium", "Titanium Ti-6Al-4V", 113_800.0, 0.342, 4.43),
    ("pla", "PLA (3D print)", 3_500.0, 0.36, 1.24),
    ("petg", "PETG (3D print)", 2_100.0, 0.38, 1.27),
    ("abs", "ABS (3D print)", 2_200.0, 0.35, 1.05),
    ("nylon", "Nylon PA12", 1_700.0, 0.39, 1.01),
];

impl Material {
    pub fn preset(id: &str) -> Option<Material> {
        MATERIALS.iter().find(|m| m.0 == id).map(|m| Material { preset: m.0.into(), e: m.2, nu: m.3, density: m.4 })
    }
    pub fn label(&self) -> String {
        MATERIALS.iter().find(|m| m.0 == self.preset).map(|m| m.1.to_string()).unwrap_or_else(|| "Custom".into())
    }
}

impl Default for Material {
    fn default() -> Self {
        Material::preset("aluminium-6061").unwrap()
    }
}

/// Optimization settings.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Settings {
    /// Material kept, as a fraction of the design body's volume (0.05 … 0.9).
    pub volfrac: f64,
    /// Element edge length, mm; 0 = from `max_elements`.
    pub element_size: f64,
    /// Upper bound on elements when `element_size` is 0 (and a hard cap otherwise).
    pub max_elements: usize,
    /// Smallest member thickness, mm; 0 = about three elements.
    pub min_member: f64,
    /// Solid layer kept under supports and loads, mm; 0 = 1.5 elements.
    pub pad: f64,
    pub max_iterations: usize,
    /// Stop optimizing after this many seconds (the job then still makes the mesh).
    pub time_limit: f64,
    /// Hand Optimize to TopOpt Solver (all CPU cores) when it is running.
    pub use_solver: bool,
    /// SIMP penalty.
    pub penal: f64,
    /// Final sharpness of the projection.
    pub beta_max: f64,
    /// OC move limit.
    pub move_limit: f64,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            volfrac: 0.3,
            element_size: 0.0,
            max_elements: 40_000,
            min_member: 0.0,
            pad: 0.0,
            max_iterations: 150,
            time_limit: 480.0,
            use_solver: true,
            penal: 3.0,
            beta_max: 16.0,
            move_limit: 0.2,
        }
    }
}

/// Settings of the mesh made from the result.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct MeshSettings {
    pub triangles: usize,
    pub smoothing: usize,
    /// Keep vertices near the original body's faces on those faces.
    pub keep_faces: bool,
    pub remove_islands: bool,
    pub fill_cavities: bool,
    /// `create-body` also converts the mesh body to a solid.
    pub convert_solid: bool,
    /// Density threshold of the surface.
    pub iso: f64,
    /// Smoother: `taubin`, `hc` or `bilateral` (see [`SMOOTHERS`]).
    pub smoother: String,
    /// Restore the volume smoothing takes, by offsetting the free surface along its normals.
    pub preserve_volume: bool,
    /// Keep the sharp edges and corners of the design body sharp (with `keep_faces`).
    pub sharp_edges: bool,
    /// Sample the density field this many times finer (1…3) before making the surface: 2 gives
    /// round members instead of faceted ones, for about 8× the meshing work.
    pub upsample: usize,
    /// Trim thin pieces as long as the stiffness drops by at most this many percent (0 = off).
    pub trim_percent: f64,
    /// Round the joints where the optimized material meets the kept parts (the pads under
    /// supports and loads, keep-in regions) with fillets.
    pub fillets: bool,
    /// Fillet radius, mm; 0 = 1.5 elements.
    pub fillet_radius: f64,
}

/// Smoothers of the result mesh: id and label.
pub const SMOOTHERS: [(&str, &str); 3] = [("taubin", "Taubin"), ("hc", "HC (keeps volume)"), ("bilateral", "Bilateral (keeps creases)")];

impl Default for MeshSettings {
    fn default() -> Self {
        MeshSettings {
            triangles: 30_000,
            smoothing: 20,
            keep_faces: true,
            remove_islands: true,
            fill_cavities: true,
            convert_solid: true,
            iso: 0.5,
            smoother: "bilateral".into(),
            preserve_volume: true,
            sharp_edges: true,
            upsample: 2,
            trim_percent: 0.0,
            fillets: true,
            fillet_radius: 0.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Condition {
    /// Faces held in the chosen directions.
    Support { id: u32, faces: Vec<FaceRef>, dirs: [bool; 3] },
    /// A force on faces, spread evenly. `mode`: `vector` uses `force` (N, total); `normal` pushes
    /// into the faces with `magnitude` N.
    Load { id: u32, faces: Vec<FaceRef>, mode: String, force: [f64; 3], magnitude: f64, case: u32 },
    /// Faces that must not change: a solid layer of `thickness` mm (0 = the pad thickness) under
    /// them is kept out of the optimization, and the mesh keeps them exactly where they are.
    Keep {
        id: u32,
        faces: Vec<FaceRef>,
        #[serde(default)]
        thickness: f64,
    },
}

impl Condition {
    pub fn id(&self) -> u32 {
        match self {
            Condition::Support { id, .. } | Condition::Load { id, .. } | Condition::Keep { id, .. } => *id,
        }
    }
    pub fn faces(&self) -> &[FaceRef] {
        match self {
            Condition::Support { faces, .. } | Condition::Load { faces, .. } | Condition::Keep { faces, .. } => faces,
        }
    }
    /// Type and value columns of the panel.
    pub fn summary(&self) -> (String, String) {
        match self {
            Condition::Support { dirs, .. } => {
                let held: Vec<&str> = ["X", "Y", "Z"].iter().zip(dirs).filter(|(_, d)| **d).map(|(a, _)| *a).collect();
                ("Support".into(), if held.len() == 3 { "fixed".into() } else { format!("holds {}", held.join(", ")) })
            }
            Condition::Load { mode, force, magnitude, case, .. } => {
                let what = if mode == "normal" {
                    format!("{} N, into face", fmt_n(*magnitude))
                } else {
                    let total = (force[0].powi(2) + force[1].powi(2) + force[2].powi(2)).sqrt();
                    let axes: Vec<(usize, f64)> = force.iter().copied().enumerate().filter(|(_, f)| *f != 0.0).collect();
                    let sign = |f: f64| if f < 0.0 { "−" } else { "+" };
                    match axes.as_slice() {
                        [(a, f)] => format!("{} N, {}{}", fmt_n(f.abs()), sign(*f), ["X", "Y", "Z"][*a]),
                        _ => format!(
                            "{} N ({})",
                            fmt_n(total),
                            axes.iter().map(|(a, f)| format!("{}{} {}", ["X", "Y", "Z"][*a], sign(*f), fmt_n(f.abs()))).collect::<Vec<_>>().join(", ")
                        ),
                    }
                };
                (if *case > 1 { format!("Load, case {case}") } else { "Load".into() }, what)
            }
            Condition::Keep { thickness, .. } => (
                "Keep surface".into(),
                if *thickness > 0.0 { format!("unchanged, solid {} mm deep", fmt_n(*thickness)) } else { "unchanged, solid pad under it".into() },
            ),
        }
    }

    /// The faces, readable.
    pub fn where_(&self) -> String {
        self.faces().iter().map(FaceRef::display).collect::<Vec<_>>().join(", ")
    }

    /// Total force magnitude of a load (0 for a support).
    pub fn force_magnitude(&self) -> f64 {
        match self {
            Condition::Load { mode, force, magnitude, .. } => {
                if mode == "normal" { *magnitude } else { (force[0].powi(2) + force[1].powi(2) + force[2].powi(2)).sqrt() }
            }
            Condition::Support { .. } | Condition::Keep { .. } => 0.0,
        }
    }
}

/// A number for people: 2000 → "2 000", 0.5 → "0.5".
pub fn fmt_n(x: f64) -> String {
    if x.fract().abs() > 1e-9 && x.abs() < 1e4 {
        let s = format!("{x:.2}");
        return s.trim_end_matches('0').trim_end_matches('.').to_string();
    }
    let digits = format!("{:.0}", x.abs());
    let mut out = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push('\u{2009}');
        }
        out.push(c);
    }
    if x < 0.0 { format!("−{out}") } else { out }
}

/// Outcome of the last run, for the panel and `create-body`.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct RunResult {
    pub outcome: String,
    pub iterations: usize,
    pub compliance: f64,
    pub volume_fraction: f64,
    pub mass_g: f64,
    pub mesh_volume_mm3: f64,
    pub triangles: usize,
    pub genus: i64,
    pub elements: usize,
    pub grid: [usize; 3],
    pub h: f64,
    pub seconds: f64,
    /// Mesh file (PLY) inside the `results` folder grant.
    pub stl: String,
    /// `create-body` already inserted this result.
    pub inserted: bool,
    pub warnings: Vec<String>,
    /// Share of the material trimmed as dead weight (`MeshSettings::trim_percent`).
    pub trimmed: f64,
    /// What solved the linear systems ("CPU, 6 threads", or a GPU's name).
    pub solver: String,
    /// The result checked against its loads: the densities made solid (at or above the iso
    /// value) or empty, analysed once more. Worst over the load cases.
    pub max_displacement_mm: f64,
    /// Von Mises stress in that material, MPa: 99 % of it is below `stress_p99_mpa`;
    /// `stress_max_mpa` is the single highest element (the voxels exaggerate it at supports and
    /// corners).
    pub stress_p99_mpa: f64,
    pub stress_max_mpa: f64,
    /// The same per load case, in case-number order.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub cases: Vec<CaseResult>,
}

/// One load case of [`RunResult::cases`].
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CaseResult {
    pub case: u32,
    /// F·u, N·mm.
    pub compliance: f64,
    pub max_displacement_mm: f64,
    pub stress_p99_mpa: f64,
    pub stress_max_mpa: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Study {
    pub version: u32,
    /// The design body: the space material may occupy.
    pub body: Option<String>,
    /// Bodies whose space must stay empty.
    pub keep_out: Vec<String>,
    pub material: Material,
    pub settings: Settings,
    pub mesh: MeshSettings,
    pub conditions: Vec<Condition>,
    pub next_id: u32,
    pub result: Option<RunResult>,
    /// A job sent to TopOpt Solver whose result has not been taken yet.
    pub pending: Option<crate::solver::Pending>,
    /// Why the last solver job failed, until the next run.
    pub last_error: Option<String>,
}

impl Default for Study {
    fn default() -> Self {
        Study {
            version: 1,
            body: None,
            keep_out: Vec::new(),
            material: Material::default(),
            settings: Settings::default(),
            mesh: MeshSettings::default(),
            conditions: Vec::new(),
            next_id: 1,
            result: None,
            pending: None,
            last_error: None,
        }
    }
}

/// A panel value as the host sends it: a bare JSON scalar or a one-key object such as
/// `{"number": 0.3}` / `{"Number": 0.3}`.
fn scalar(v: &Value) -> &Value {
    match v {
        Value::Object(m) if m.len() == 1 => m.values().next().unwrap(),
        _ => v,
    }
}

impl Study {
    pub fn add(&mut self, mut c: Condition, replace: Option<u32>) -> u32 {
        let id = match replace.filter(|r| self.conditions.iter().any(|x| x.id() == *r)) {
            Some(r) => r,
            None => {
                self.next_id += 1;
                self.next_id - 1
            }
        };
        match &mut c {
            Condition::Support { id: i, .. } | Condition::Load { id: i, .. } | Condition::Keep { id: i, .. } => *i = id,
        }
        match self.conditions.iter().position(|x| x.id() == id) {
            Some(p) => self.conditions[p] = c,
            None => self.conditions.push(c),
        }
        id
    }

    pub fn remove(&mut self, id: u32) -> bool {
        let before = self.conditions.len();
        self.conditions.retain(|c| c.id() != id);
        before != self.conditions.len()
    }

    /// Applies one panel field change: `key` is the field's path tail (`settings/volfrac`,
    /// `material/preset`, `mesh/triangles`, …). Returns false for an unknown key.
    pub fn set_option(&mut self, key: &str, value: &Value) -> Result<(), String> {
        let v = scalar(value);
        let num = || v.as_f64().ok_or_else(|| format!("{key}: expected a number"));
        let flag = || v.as_bool().ok_or_else(|| format!("{key}: expected on/off"));
        let pos = |x: f64, lo: f64, hi: f64| -> Result<f64, String> {
            if x.is_finite() && x >= lo && x <= hi { Ok(x) } else { Err(format!("{key}: {x} is outside {lo} … {hi}")) }
        };
        let s = &mut self.settings;
        let m = &mut self.mesh;
        match key.trim_start_matches('/') {
            "settings/volfrac" => s.volfrac = pos(num()?, 0.02, 0.95)?,
            "settings/volume_percent" => s.volfrac = pos(num()?, 2.0, 95.0)? / 100.0,
            "settings/element_size" => s.element_size = pos(num()?, 0.0, 1e4)?,
            "settings/max_elements" => s.max_elements = pos(num()?, 1000.0, 400_000.0)? as usize,
            "settings/min_member" => s.min_member = pos(num()?, 0.0, 1e4)?,
            "settings/pad" => s.pad = pos(num()?, 0.0, 1e4)?,
            "settings/max_iterations" => s.max_iterations = pos(num()?, 1.0, 2000.0)? as usize,
            "settings/time_limit" => s.time_limit = pos(num()?, 10.0, 86_400.0)?,
            "settings/use_solver" => s.use_solver = flag()?,
            "settings/penal" => s.penal = pos(num()?, 1.0, 6.0)?,
            "settings/beta_max" => s.beta_max = pos(num()?, 1.0, 128.0)?,
            "settings/move_limit" => s.move_limit = pos(num()?, 0.01, 0.5)?,
            "material/preset" | "material/label" => {
                let text = v.as_str().ok_or("material: expected a name")?;
                let id = MATERIALS.iter().find(|m| m.0 == text || m.1 == text).map(|m| m.0).unwrap_or("custom");
                match Material::preset(id) {
                    Some(mat) => self.material = mat,
                    None => self.material.preset = "custom".into(),
                }
            }
            "material/e" => {
                self.material.e = pos(num()?, 1.0, 1e7)?;
                self.material.preset = "custom".into();
            }
            "material/nu" => {
                self.material.nu = pos(num()?, 0.0, 0.49)?;
                self.material.preset = "custom".into();
            }
            "material/density" => {
                self.material.density = pos(num()?, 0.01, 30.0)?;
                self.material.preset = "custom".into();
            }
            "mesh/triangles" => m.triangles = pos(num()?, 500.0, MAX_TRIANGLES as f64)? as usize,
            "mesh/smoothing" => m.smoothing = pos(num()?, 0.0, 100.0)? as usize,
            "mesh/keep_faces" => m.keep_faces = flag()?,
            "mesh/remove_islands" => m.remove_islands = flag()?,
            "mesh/fill_cavities" => m.fill_cavities = flag()?,
            "mesh/convert_solid" => m.convert_solid = flag()?,
            "mesh/iso" => m.iso = pos(num()?, 0.1, 0.9)?,
            "mesh/smoother" => {
                let text = v.as_str().ok_or("smoother: expected a name")?;
                let id = SMOOTHERS.iter().find(|x| x.0 == text || x.1 == text).ok_or(format!("unknown smoother `{text}`"))?.0;
                m.smoother = id.into();
            }
            // The panel's radio toggles: ticking one selects it; unticking the current one keeps
            // it (one smoother is always on).
            k if k.starts_with("mesh/smoother_") => {
                let id = &k["mesh/smoother_".len()..];
                let id = SMOOTHERS.iter().find(|x| x.0 == id).ok_or(format!("unknown smoother `{id}`"))?.0;
                if flag()? {
                    m.smoother = id.into();
                }
            }
            "mesh/preserve_volume" => m.preserve_volume = flag()?,
            "mesh/sharp_edges" => m.sharp_edges = flag()?,
            "mesh/upsample" => m.upsample = pos(num()?, 1.0, 3.0)?.round() as usize,
            "mesh/trim_percent" => m.trim_percent = pos(num()?, 0.0, 10.0)?,
            "mesh/fillets" => m.fillets = flag()?,
            "mesh/fillet_radius" => m.fillet_radius = pos(num()?, 0.0, 1e4)?,
            other => return Err(format!("unknown setting `{other}`")),
        }
        Ok(())
    }

    /// Whether the mesh is pinned to the design body's faces: asked for (`mesh.keep_faces`), or
    /// needed by a kept surface.
    pub fn pins_faces(&self) -> bool {
        self.mesh.keep_faces || self.conditions.iter().any(|c| matches!(c, Condition::Keep { .. }))
    }

    /// What is missing before the study can run.
    pub fn problems(&self) -> Vec<String> {
        let mut p = Vec::new();
        if self.body.is_none() {
            p.push("choose the design body".into());
        }
        if !self.conditions.iter().any(|c| matches!(c, Condition::Support { .. })) {
            p.push("add at least one support".into());
        }
        if !self.conditions.iter().any(|c| matches!(c, Condition::Load { .. })) {
            p.push("add at least one load".into());
        }
        if let Some(b) = &self.body {
            for c in &self.conditions {
                if c.faces().iter().any(|f| &f.body != b) {
                    p.push(format!("item {} uses a face of another body than {b}", c.id()));
                }
            }
        }
        p
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn options_and_round_trip() {
        let mut s = Study::default();
        s.set_option("/settings/volfrac", &json!({"number": 0.25})).unwrap();
        s.set_option("material/preset", &json!("Steel (structural)")).unwrap();
        assert_eq!(s.material.e, 200_000.0);
        s.set_option("material/e", &json!(210_000.0)).unwrap();
        assert_eq!(s.material.preset, "custom");
        s.set_option("mesh/convert_solid", &json!({"Flag": false})).unwrap();
        assert!(s.set_option("settings/volfrac", &json!(2.0)).is_err());
        assert!(s.set_option("nope", &json!(1)).is_err());
        // The smoother: radio toggles from the panel (as the host sends them), or a name.
        assert_eq!(s.mesh.smoother, "bilateral");
        s.set_option("/mesh/smoother_hc", &json!({"flag": true})).unwrap();
        assert_eq!(s.mesh.smoother, "hc");
        s.set_option("/mesh/smoother_hc", &json!({"flag": false})).unwrap();
        assert_eq!(s.mesh.smoother, "hc", "unticking the current smoother keeps it");
        s.set_option("/mesh/smoother_bilateral", &json!({"flag": false})).unwrap();
        assert_eq!(s.mesh.smoother, "hc");
        s.set_option("/mesh/smoother_taubin", &json!({"flag": true})).unwrap();
        assert_eq!(s.mesh.smoother, "taubin");
        s.set_option("mesh/smoother", &json!({"text": "HC (keeps volume)"})).unwrap();
        assert_eq!(s.mesh.smoother, "hc");
        assert!(s.set_option("/mesh/smoother_fast", &json!({"flag": true})).is_err());
        let face = FaceRef { body: "B".into(), index: 2, point: None, label: None };
        let a = s.add(Condition::Support { id: 0, faces: vec![face.clone()], dirs: [true; 3] }, None);
        let b = s.add(Condition::Load { id: 0, faces: vec![face], mode: "vector".into(), force: [0., 0., -1.], magnitude: 0., case: 1 }, None);
        assert_ne!(a, b);
        s.body = Some("B".into());
        assert!(s.problems().is_empty());
        let back: Study = serde_json::from_value(serde_json::to_value(&s).unwrap()).unwrap();
        assert_eq!(back, s);
        assert!(s.remove(a));
        assert_eq!(s.problems(), vec!["add at least one support".to_string()]);
        // Older or partial documents fill in defaults.
        let partial: Study = serde_json::from_value(json!({"body": "B"})).unwrap();
        assert_eq!(partial.settings, Settings::default());
    }
}
