//! From a study plus voxels to the topopt-core problem: element kinds, supports and loads on the
//! nodes next to the chosen faces, solid pads under them and under kept surfaces.

use crate::study::{Condition, Study};
use crate::voxel::{Hit, Voxels};
use std::collections::HashMap;
use topopt_core::{ElementKind, Grid, LoadCase, NodeFix, NodeForce, Problem};

pub struct Built {
    pub problem: Problem,
    /// Solid pad cells (anchors for the mesh clean-up).
    pub pads: Vec<bool>,
    /// The cells under kept surfaces (a subset of `pads`): the mesh is not simplified there.
    pub kept: Vec<bool>,
    pub warnings: Vec<String>,
    /// Nodes used per condition id, for the report.
    pub nodes_per_condition: Vec<(u32, usize)>,
}

/// Face samples by face index.
fn samples(hits: &[Hit]) -> HashMap<i64, Vec<Hit>> {
    let mut m: HashMap<i64, Vec<Hit>> = HashMap::new();
    for h in hits {
        m.entry(h.face).or_default().push(*h);
    }
    m
}

/// Nodes within `reach` of any sample, that belong to at least one inside cell.
fn nodes_near(grid: &Grid, inside: &[bool], pts: &[[f64; 3]], reach: f64) -> Vec<usize> {
    let h = grid.h;
    let mut on_solid = vec![false; grid.nnodes()];
    for e in 0..grid.nel() {
        if inside[e] {
            let [i, j, k] = grid.elem_ijk(e);
            for n in grid.elem_nodes(i, j, k) {
                on_solid[n] = true;
            }
        }
    }
    let mut picked = vec![false; grid.nnodes()];
    let r = (reach / h).ceil() as i64 + 1;
    for p in pts {
        let c: [i64; 3] = std::array::from_fn(|d| ((p[d] - grid.origin[d]) / h).round() as i64);
        for dk in -r..=r {
            for dj in -r..=r {
                for di in -r..=r {
                    let (i, j, k) = (c[0] + di, c[1] + dj, c[2] + dk);
                    if i < 0 || j < 0 || k < 0 || i > grid.nx as i64 || j > grid.ny as i64 || k > grid.nz as i64 {
                        continue;
                    }
                    let n = grid.node(i as usize, j as usize, k as usize);
                    let q = grid.node_pos(n);
                    let d2 = (q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2) + (q[2] - p[2]).powi(2);
                    if d2 <= reach * reach && on_solid[n] {
                        picked[n] = true;
                    }
                }
            }
        }
    }
    (0..grid.nnodes()).filter(|&n| picked[n]).collect()
}

/// Inside cells whose centre is within `depth` of a sample.
pub(crate) fn cells_near(grid: &Grid, inside: &[bool], pts: &[[f64; 3]], depth: f64, out: &mut [bool]) {
    let h = grid.h;
    let r = (depth / h).ceil() as i64 + 1;
    for p in pts {
        let c: [i64; 3] = std::array::from_fn(|d| ((p[d] - grid.origin[d]) / h - 0.5).round() as i64);
        for dk in -r..=r {
            for dj in -r..=r {
                for di in -r..=r {
                    let (i, j, k) = (c[0] + di, c[1] + dj, c[2] + dk);
                    if i < 0 || j < 0 || k < 0 || i >= grid.nx as i64 || j >= grid.ny as i64 || k >= grid.nz as i64 {
                        continue;
                    }
                    let e = grid.elem(i as usize, j as usize, k as usize);
                    let q = grid.elem_center(e);
                    let d2 = (q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2) + (q[2] - p[2]).powi(2);
                    if inside[e] && d2 <= depth * depth {
                        out[e] = true;
                    }
                }
            }
        }
    }
}

/// Builds the problem. `keep_out`: cells inside keep-out bodies.
pub fn build(study: &Study, grid: Grid, design: &Voxels, keep_out: &[bool]) -> Result<Built, String> {
    let mut warnings = Vec::new();
    let h = grid.h;
    let inside: Vec<bool> = design.inside.iter().zip(keep_out).map(|(&i, &k)| i && !k).collect();
    if !inside.iter().any(|&i| i) {
        return Err("the design body produced no elements: is it a closed solid?".into());
    }
    if design.unbalanced_rays > 0 {
        warnings.push(format!("{} rays crossed the body an odd number of times; check that it is closed", design.unbalanced_rays));
    }
    let by_face = samples(&design.hits);
    let reach = 1.0 * h;
    let pad = if study.settings.pad > 0.0 { study.settings.pad } else { 1.5 * h };
    let mut pads = vec![false; grid.nel()];
    let mut kept = vec![false; grid.nel()];
    let mut problem = Problem::new(grid);
    problem.material = topopt_core::Material { e: study.material.e, nu: study.material.nu };
    let mut cases: Vec<(u32, Vec<NodeForce>)> = Vec::new();
    let mut nodes_per_condition = Vec::new();
    for c in &study.conditions {
        let mut pts: Vec<[f64; 3]> = Vec::new();
        let mut normal = [0.0; 3];
        for f in c.faces() {
            if let Some(list) = by_face.get(&f.index) {
                for hit in list {
                    pts.push(hit.point);
                    for d in 0..3 {
                        normal[d] += hit.normal[d];
                    }
                }
            }
            if let Some(p) = f.point {
                pts.push(p);
            }
        }
        if pts.is_empty() {
            return Err(format!("item {}: its faces were not found on the body (try a smaller element size)", c.id()));
        }
        let nodes = nodes_near(&grid, &inside, &pts, reach);
        if nodes.is_empty() {
            return Err(format!("item {}: no element next to its faces; use a smaller element size", c.id()));
        }
        nodes_per_condition.push((c.id(), nodes.len()));
        let depth = match c {
            Condition::Keep { thickness, .. } if *thickness > 0.0 => thickness.max(0.5 * h),
            _ => pad,
        };
        cells_near(&grid, &inside, &pts, depth, &mut pads);
        if matches!(c, Condition::Keep { .. }) {
            cells_near(&grid, &inside, &pts, depth, &mut kept);
        }
        match c {
            // Only material: no support, no force.
            Condition::Keep { .. } => {}
            Condition::Support { dirs, .. } => {
                problem.supports.extend(nodes.iter().map(|&node| NodeFix { node, dirs: *dirs }));
            }
            Condition::Load { mode, force, magnitude, case, .. } => {
                let total = if mode == "normal" {
                    let len = (normal[0].powi(2) + normal[1].powi(2) + normal[2].powi(2)).sqrt();
                    if len < 1e-12 {
                        return Err(format!("item {}: the faces have no common normal; give the force as a vector", c.id()));
                    }
                    // Into the faces: against the outward normal.
                    [-magnitude * normal[0] / len, -magnitude * normal[1] / len, -magnitude * normal[2] / len]
                } else {
                    *force
                };
                let k = nodes.len() as f64;
                let forces = nodes.iter().map(|&node| NodeForce { node, f: [total[0] / k, total[1] / k, total[2] / k] });
                match cases.iter_mut().find(|(n, _)| n == case) {
                    Some((_, list)) => list.extend(forces),
                    None => cases.push((*case, forces.collect())),
                }
            }
        }
    }
    cases.sort_by_key(|(n, _)| *n);
    problem.load_cases = cases.into_iter().map(|(_, forces)| LoadCase { weight: 1.0, forces, supports: None, relative: false }).collect();
    problem.kinds = (0..grid.nel())
        .map(|e| {
            if !inside[e] {
                ElementKind::Void
            } else if pads[e] {
                ElementKind::Solid
            } else {
                ElementKind::Design
            }
        })
        .collect();
    let design_cells = problem.kinds.iter().filter(|k| **k == ElementKind::Design).count();
    let solid_cells = pads.iter().filter(|&&p| p).count();
    if design_cells == 0 {
        return Err("pads cover the whole body: lower the pad thickness".into());
    }
    let solid_share = solid_cells as f64 / (design_cells + solid_cells) as f64;
    if solid_share > 0.5 * study.settings.volfrac.max(0.1) {
        warnings.push(format!(
            "solid pads under supports, loads and kept surfaces are {:.0} % of the body; the volume target applies to the rest",
            solid_share * 100.0
        ));
    }
    problem.validate()?;
    Ok(Built { problem, pads, kept, warnings, nodes_per_condition })
}
