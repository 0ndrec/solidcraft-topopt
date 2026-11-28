//! A study described without a CAD app: the input format of `topopt-solver solve` and of the
//! library's [`prepare_spec`]. The design domain is a box or a triangle mesh (STL, or positions
//! and triangles inline); supports and loads act on the grid nodes inside regions (boxes, spheres,
//! planes), so a pipeline needs no face ids.
//!
//! ```json
//! {
//!   "name": "bracket",
//!   "domain": {"stl": {"path": "bracket.stl"}},
//!   "keep_out": [{"box": {"min": [0, 0, 0], "max": [10, 10, 40]}}],
//!   "material": {"preset": "aluminium-6061"},
//!   "settings": {"volfrac": 0.3, "max_elements": 60000},
//!   "supports": [{"region": {"plane": {"point": [0, 0, 0], "normal": [1, 0, 0]}}}],
//!   "loads": [{"region": {"sphere": {"center": [100, 20, 0], "radius": 6}}, "force": [0, 0, -1000]}]
//! }
//! ```

use crate::study::{Material, MeshSettings, Settings, Study};
use crate::voxel::{Caster, Hit};
use crate::{BodyInfo, Hooks, Job, Mesh};
use parry3d_f64::math::Vector;
use parry3d_f64::query::{Ray, RayCast};
use parry3d_f64::shape::TriMesh;
use serde::{Deserialize, Serialize};
use topopt_core::{ElementKind, Grid, LoadCase, NodeFix, NodeForce, Problem};

/// A solid: an axis-aligned box, an STL file, or an inline triangle mesh (closed, mm).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum Shape {
    Box { min: [f64; 3], max: [f64; 3] },
    /// An STL file (binary or ASCII). `path` is relative to the spec file; `scale` multiplies its
    /// units into mm (1000 for metres, 25.4 for inches).
    Stl {
        path: String,
        #[serde(default = "one")]
        scale: f64,
    },
    Mesh { positions: Vec<[f64; 3]>, triangles: Vec<[u32; 3]> },
}

fn one() -> f64 {
    1.0
}

/// Where a support or load acts: the grid nodes of material inside it (preferring the nodes on
/// the part's surface). `tolerance` widens it; default 0.6 element so a region drawn exactly on a
/// face still catches the nodes near that face.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum Region {
    Box {
        min: [f64; 3],
        max: [f64; 3],
        #[serde(default)]
        tolerance: Option<f64>,
    },
    Sphere {
        center: [f64; 3],
        radius: f64,
        #[serde(default)]
        tolerance: Option<f64>,
    },
    /// The slab within `tolerance` of a plane, optionally limited to a box.
    Plane {
        point: [f64; 3],
        normal: [f64; 3],
        #[serde(default)]
        tolerance: Option<f64>,
        #[serde(default)]
        within: Option<([f64; 3], [f64; 3])>,
    },
}

impl Region {
    fn tolerance(&self) -> Option<f64> {
        match self {
            Region::Box { tolerance, .. } | Region::Sphere { tolerance, .. } | Region::Plane { tolerance, .. } => *tolerance,
        }
    }

    /// Whether point `p` is in the region widened by `tol`.
    pub fn contains(&self, p: [f64; 3], default_tol: f64) -> bool {
        let tol = self.tolerance().unwrap_or(default_tol);
        let in_box = |lo: &[f64; 3], hi: &[f64; 3]| (0..3).all(|d| p[d] >= lo[d] - tol && p[d] <= hi[d] + tol);
        match self {
            Region::Box { min, max, .. } => in_box(min, max),
            Region::Sphere { center, radius, .. } => {
                ((p[0] - center[0]).powi(2) + (p[1] - center[1]).powi(2) + (p[2] - center[2]).powi(2)).sqrt() <= radius + tol
            }
            Region::Plane { point, normal, within, .. } => {
                let len = (normal[0].powi(2) + normal[1].powi(2) + normal[2].powi(2)).sqrt().max(1e-300);
                let dist = ((p[0] - point[0]) * normal[0] + (p[1] - point[1]) * normal[1] + (p[2] - point[2]) * normal[2]) / len;
                dist.abs() <= tol && within.as_ref().is_none_or(|(lo, hi)| in_box(lo, hi))
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SupportSpec {
    pub region: Region,
    /// Held in x, y, z (default: all).
    #[serde(default = "all_dirs")]
    pub dirs: [bool; 3],
}

fn all_dirs() -> [bool; 3] {
    [true; 3]
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LoadSpec {
    pub region: Region,
    /// Total force, N, spread evenly over the region's nodes.
    pub force: [f64; 3],
    /// Load case (loads with the same number act together).
    #[serde(default = "case_one")]
    pub case: u32,
}

fn case_one() -> u32 {
    1
}

/// A material: a preset from [`crate::MATERIALS`] and/or explicit values (which win).
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MaterialSpec {
    #[serde(default)]
    pub preset: Option<String>,
    /// Young's modulus, MPa.
    #[serde(default)]
    pub e: Option<f64>,
    #[serde(default)]
    pub nu: Option<f64>,
    /// g/cm³.
    #[serde(default)]
    pub density: Option<f64>,
}

impl MaterialSpec {
    pub fn resolve(&self) -> Result<Material, String> {
        let mut m = match &self.preset {
            Some(p) => Material::preset(p).ok_or_else(|| {
                format!("unknown material preset `{p}`; known: {}", crate::MATERIALS.iter().map(|m| m.0).collect::<Vec<_>>().join(", "))
            })?,
            None => Material::default(),
        };
        if self.e.is_some() || self.nu.is_some() || self.density.is_some() {
            m.preset = "custom".into();
        }
        if let Some(e) = self.e {
            m.e = e;
        }
        if let Some(nu) = self.nu {
            m.nu = nu;
        }
        if let Some(d) = self.density {
            m.density = d;
        }
        if !(m.e > 0.0 && (0.0..0.5).contains(&m.nu) && m.density > 0.0) {
            return Err("material: need E > 0, 0 ≤ nu < 0.5 and density > 0".into());
        }
        Ok(m)
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Spec {
    /// A name for reports and the default output file names.
    #[serde(default)]
    pub name: Option<String>,
    /// The design domain: material may only be placed inside it.
    pub domain: Shape,
    /// Solids whose space must stay empty.
    #[serde(default)]
    pub keep_out: Vec<Shape>,
    /// Regions that must stay solid (besides the pads under supports and loads).
    #[serde(default)]
    pub keep_in: Vec<Region>,
    #[serde(default)]
    pub material: MaterialSpec,
    /// Optimization settings (same fields and defaults as the plugin's; see [`Settings`]).
    #[serde(default)]
    pub settings: Settings,
    /// Mesh settings (see [`MeshSettings`]).
    #[serde(default)]
    pub mesh: MeshSettings,
    pub supports: Vec<SupportSpec>,
    pub loads: Vec<LoadSpec>,
}

impl Spec {
    /// A small cantilever, as a starting point.
    pub fn example() -> Spec {
        Spec {
            name: Some("cantilever".into()),
            domain: Shape::Box { min: [0.0, 0.0, 0.0], max: [120.0, 40.0, 40.0] },
            keep_out: Vec::new(),
            keep_in: Vec::new(),
            material: MaterialSpec { preset: Some("aluminium-6061".into()), ..Default::default() },
            settings: Settings { max_elements: 30_000, ..Default::default() },
            mesh: MeshSettings::default(),
            supports: vec![SupportSpec { region: Region::Plane { point: [0.0; 3], normal: [1.0, 0.0, 0.0], tolerance: None, within: None }, dirs: [true; 3] }],
            loads: vec![LoadSpec {
                region: Region::Box { min: [120.0, 0.0, 0.0], max: [120.0, 40.0, 0.0], tolerance: None },
                force: [0.0, 0.0, -1000.0],
                case: 1,
            }],
        }
    }
}

/// The box as a closed outward mesh.
pub fn box_mesh(lo: [f64; 3], hi: [f64; 3]) -> Mesh {
    let c = |i: usize| [if i & 1 == 1 { hi[0] } else { lo[0] }, if i & 2 == 2 { hi[1] } else { lo[1] }, if i & 4 == 4 { hi[2] } else { lo[2] }];
    Mesh {
        positions: vec![c(0), c(1), c(3), c(2), c(4), c(5), c(7), c(6)],
        triangles: vec![[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]],
    }
}

/// A closed triangle mesh as a ray caster. Entry and exit come from the parity of the crossings
/// along the whole line, so the mesh's winding does not matter.
pub struct MeshCaster {
    mesh: TriMesh,
    reach: f64,
}

impl MeshCaster {
    pub fn new(m: &Mesh) -> Result<Self, String> {
        let verts = m.positions.iter().map(|p| Vector::new(p[0], p[1], p[2])).collect();
        let mesh = TriMesh::new(verts, m.triangles.clone()).map_err(|e| format!("mesh: {e:?}"))?;
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for p in &m.positions {
            for d in 0..3 {
                lo[d] = lo[d].min(p[d]);
                hi[d] = hi[d].max(p[d]);
            }
        }
        let diag = ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt();
        Ok(MeshCaster { mesh, reach: 2.0 * diag + 1.0 })
    }

    /// Distances of every crossing along the line through `o` with direction `d`, measured from
    /// `o - reach·d`.
    fn crossings(&self, o: [f64; 3], d: [f64; 3]) -> Vec<f64> {
        let start = [o[0] - self.reach * d[0], o[1] - self.reach * d[1], o[2] - self.reach * d[2]];
        let eps = 1e-9 * self.reach;
        let mut out = Vec::new();
        let mut t0 = 0.0;
        while out.len() < 1024 {
            let p = [start[0] + t0 * d[0], start[1] + t0 * d[1], start[2] + t0 * d[2]];
            let ray = Ray::new(Vector::new(p[0], p[1], p[2]), Vector::new(d[0], d[1], d[2]));
            match self.mesh.cast_local_ray(&ray, 3.0 * self.reach, false) {
                Some(t) => {
                    out.push(t0 + t);
                    t0 += t + eps;
                }
                None => break,
            }
        }
        out
    }
}

impl Caster for MeshCaster {
    fn cast(&mut self, rays: &[([f64; 3], [f64; 3])]) -> Result<Vec<Option<Hit>>, String> {
        Ok(rays
            .iter()
            .map(|&(o, d)| {
                let all = self.crossings(o, d);
                let k = all.iter().position(|&t| t >= self.reach)?;
                let t = all[k] - self.reach;
                // Even crossings enter the solid, odd ones leave it.
                let s = if k % 2 == 0 { -1.0 } else { 1.0 };
                Some(Hit { point: [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]], normal: [s * d[0], s * d[1], s * d[2]], face: -1 })
            })
            .collect())
    }
}

fn bounds(m: &Mesh) -> ([f64; 3], [f64; 3]) {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for p in &m.positions {
        for d in 0..3 {
            lo[d] = lo[d].min(p[d]);
            hi[d] = hi[d].max(p[d]);
        }
    }
    (lo, hi)
}

/// A shape as a mesh; `read` returns a file's bytes for [`Shape::Stl`].
pub fn shape_mesh(shape: &Shape, read: &dyn Fn(&str) -> Result<Vec<u8>, String>) -> Result<Mesh, String> {
    match shape {
        Shape::Box { min, max } => {
            if (0..3).any(|d| max[d] <= min[d]) {
                return Err("a box needs max > min on every axis".into());
            }
            Ok(box_mesh(*min, *max))
        }
        Shape::Stl { path, scale } => topopt_mesh::read_stl(&read(path)?, *scale).map_err(|e| format!("{path}: {e}")),
        Shape::Mesh { positions, triangles } => {
            if triangles.iter().flatten().any(|&i| i as usize >= positions.len()) {
                return Err("mesh: a triangle index is out of range".into());
            }
            Ok(Mesh { positions: positions.clone(), triangles: triangles.clone() })
        }
    }
}

/// Builds a [`Job`] from a spec: voxelize, select the regions' nodes, assemble the problem.
pub fn prepare_spec(spec: &Spec, read: &dyn Fn(&str) -> Result<Vec<u8>, String>, hooks: &mut dyn Hooks) -> Result<Job, String> {
    if spec.supports.is_empty() {
        return Err("the spec has no supports: the part would float".into());
    }
    if spec.loads.is_empty() {
        return Err("the spec has no loads".into());
    }
    let material = spec.material.resolve()?;
    let domain = shape_mesh(&spec.domain, read)?;
    let (lo, hi) = bounds(&domain);
    let report = domain.validate(false);
    let mut warnings = Vec::new();
    if report.boundary_edges > 0 || report.nonmanifold_edges > 0 {
        warnings.push(format!(
            "the domain mesh is not closed ({} open, {} non-manifold edges): inside/outside may be wrong near them",
            report.boundary_edges, report.nonmanifold_edges
        ));
    }
    let body = BodyInfo { lo, hi, volume: domain.volume().abs() };
    let s = &spec.settings;
    let (grid, plan_warnings) = crate::voxel::plan_grid(lo, hi, s.element_size, s.max_elements);
    warnings.extend(plan_warnings);
    hooks.progress(0.0, &format!("voxelizing ({}×{}×{}, {:.2} mm)", grid.nx, grid.ny, grid.nz, grid.h));
    let vox = crate::voxel::voxelize(&mut MeshCaster::new(&domain)?, &grid, |f| hooks.progress(0.04 * f, "voxelizing"))?;
    if vox.unbalanced_rays > 0 {
        warnings.push(format!("{} rays crossed the domain an odd number of times; check that it is closed", vox.unbalanced_rays));
    }
    let mut inside = vox.inside;
    for shape in &spec.keep_out {
        let m = shape_mesh(shape, read)?;
        let v = crate::voxel::voxelize(&mut MeshCaster::new(&m)?, &grid, |_| {})?;
        for (i, b) in inside.iter_mut().zip(&v.inside) {
            *i &= !b;
        }
    }
    if !inside.iter().any(|&i| i) {
        return Err("the domain produced no elements: is it a closed solid, in millimetres?".into());
    }
    let problem = assemble(spec, grid, &inside, material.e, material.nu, &mut warnings)?;
    let study = Study {
        body: Some(spec.name.clone().unwrap_or_else(|| "spec".into())),
        material,
        settings: spec.settings.clone(),
        mesh: spec.mesh.clone(),
        ..Default::default()
    };
    crate::finish_job(&study, body, problem, warnings, Some(crate::SurfaceRef::new(&domain.positions, &domain.triangles)))
}

/// Supports, loads and pads on the grid from the spec's regions.
fn assemble(spec: &Spec, grid: Grid, inside: &[bool], e: f64, nu: f64, warnings: &mut Vec<String>) -> Result<Problem, String> {
    let h = grid.h;
    let tol = 0.6 * h;
    // Nodes on material, and those on its surface (touching a cell that is not material).
    let mut solid_cells = vec![0u8; grid.nnodes()];
    let mut all_cells = vec![0u8; grid.nnodes()];
    for el in 0..grid.nel() {
        let [i, j, k] = grid.elem_ijk(el);
        for n in grid.elem_nodes(i, j, k) {
            all_cells[n] += 1;
            if inside[el] {
                solid_cells[n] += 1;
            }
        }
    }
    let on_solid = |n: usize| solid_cells[n] > 0;
    let on_surface = |n: usize| solid_cells[n] > 0 && (solid_cells[n] < all_cells[n] || all_cells[n] < 8);
    let select = |region: &Region, what: &str, idx: usize| -> Result<Vec<usize>, String> {
        let candidates: Vec<usize> = (0..grid.nnodes()).filter(|&n| on_solid(n) && region.contains(grid.node_pos(n), tol)).collect();
        let surface: Vec<usize> = candidates.iter().copied().filter(|&n| on_surface(n)).collect();
        let chosen = if surface.is_empty() { candidates } else { surface };
        if chosen.is_empty() {
            return Err(format!("{what} {} : its region holds no node of the part (grid step {h:.2} mm)", idx + 1));
        }
        Ok(chosen)
    };
    // Solid pads under supports and loads: cells within `pad` of a chosen node (the plugin's rule).
    let pad = if spec.settings.pad > 0.0 { spec.settings.pad } else { 1.5 * h };
    let mut pads = vec![false; grid.nel()];
    let mark_pads = |nodes: &[usize], pads: &mut Vec<bool>| {
        let pts: Vec<[f64; 3]> = nodes.iter().map(|&n| grid.node_pos(n)).collect();
        crate::build::cells_near(&grid, inside, &pts, pad, pads);
    };
    let mut problem = Problem::new(grid);
    problem.material = topopt_core::Material { e, nu };
    for (i, s) in spec.supports.iter().enumerate() {
        if !s.dirs.iter().any(|&d| d) {
            return Err(format!("support {}: holds no direction", i + 1));
        }
        let nodes = select(&s.region, "support", i)?;
        mark_pads(&nodes, &mut pads);
        problem.supports.extend(nodes.iter().map(|&node| NodeFix { node, dirs: s.dirs }));
    }
    let mut cases: Vec<(u32, Vec<NodeForce>)> = Vec::new();
    for (i, l) in spec.loads.iter().enumerate() {
        let nodes = select(&l.region, "load", i)?;
        mark_pads(&nodes, &mut pads);
        let k = nodes.len() as f64;
        let forces = nodes.iter().map(|&node| NodeForce { node, f: [l.force[0] / k, l.force[1] / k, l.force[2] / k] });
        match cases.iter_mut().find(|(c, _)| *c == l.case) {
            Some((_, list)) => list.extend(forces),
            None => cases.push((l.case, forces.collect())),
        }
    }
    cases.sort_by_key(|(c, _)| *c);
    problem.load_cases = cases.into_iter().map(|(_, forces)| LoadCase { weight: 1.0, forces, supports: None, relative: false }).collect();
    for (i, r) in spec.keep_in.iter().enumerate() {
        let before = pads.iter().filter(|&&p| p).count();
        for el in 0..grid.nel() {
            if inside[el] && r.contains(grid.elem_center(el), 0.0) {
                pads[el] = true;
            }
        }
        if pads.iter().filter(|&&p| p).count() == before {
            warnings.push(format!("keep_in {} holds no element of the part", i + 1));
        }
    }
    problem.kinds = (0..grid.nel())
        .map(|el| if !inside[el] { ElementKind::Void } else if pads[el] { ElementKind::Solid } else { ElementKind::Design })
        .collect();
    if !problem.kinds.contains(&ElementKind::Design) {
        return Err("supports, loads and keep_in regions cover the whole part: nothing left to optimize".into());
    }
    problem.validate()?;
    Ok(problem)
}
