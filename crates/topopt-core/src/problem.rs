//! The optimization problem: grid, material, passive regions, supports and load cases.
//!
//! Units are millimetres, newtons and megapascals (N/mm²), so displacements come out in mm and
//! compliance in N·mm.

use crate::grid::Grid;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Material {
    /// Young's modulus of the solid, MPa.
    pub e: f64,
    /// Poisson ratio.
    pub nu: f64,
}

impl Default for Material {
    fn default() -> Self {
        // Aluminium-like.
        Material { e: 70_000.0, nu: 0.33 }
    }
}

/// What an element may be.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[repr(u8)]
pub enum ElementKind {
    /// The optimizer chooses its density.
    Design = 0,
    /// Always solid (keep-in, load and support pads).
    Solid = 1,
    /// Always empty (keep-out, outside the design domain).
    Void = 2,
}

/// A node held in some directions.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct NodeFix {
    pub node: usize,
    /// Fixed in x, y, z.
    pub dirs: [bool; 3],
}

/// A force on a node, N.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct NodeForce {
    pub node: usize,
    pub f: [f64; 3],
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct LoadCase {
    #[serde(default = "one")]
    pub weight: f64,
    pub forces: Vec<NodeForce>,
}

fn one() -> f64 {
    1.0
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Problem {
    pub grid: Grid,
    #[serde(default)]
    pub material: Material,
    /// One [`ElementKind`] per element; empty means every element is a design element.
    #[serde(default)]
    pub kinds: Vec<ElementKind>,
    pub supports: Vec<NodeFix>,
    pub load_cases: Vec<LoadCase>,
}

impl Problem {
    pub fn new(grid: Grid) -> Self {
        Problem { grid, material: Material::default(), kinds: Vec::new(), supports: Vec::new(), load_cases: Vec::new() }
    }

    pub fn kind(&self, e: usize) -> ElementKind {
        self.kinds.get(e).copied().unwrap_or(ElementKind::Design)
    }

    /// Checks indices and that the problem can be solved.
    pub fn validate(&self) -> Result<(), String> {
        let (nel, nn) = (self.grid.nel(), self.grid.nnodes());
        if !self.kinds.is_empty() && self.kinds.len() != nel {
            return Err(format!("kinds has {} entries for {nel} elements", self.kinds.len()));
        }
        if self.supports.is_empty() {
            return Err("no supports: the part could move freely".into());
        }
        if self.load_cases.is_empty() || self.load_cases.iter().all(|c| c.forces.is_empty()) {
            return Err("no loads".into());
        }
        for s in &self.supports {
            if s.node >= nn {
                return Err(format!("support node {} is outside the grid ({nn} nodes)", s.node));
            }
        }
        for (k, c) in self.load_cases.iter().enumerate() {
            if !(c.weight >= 0.0) {
                return Err(format!("load case {k} has a negative weight"));
            }
            for f in &c.forces {
                if f.node >= nn {
                    return Err(format!("load case {k}: node {} is outside the grid", f.node));
                }
            }
        }
        if !(self.material.e > 0.0) {
            return Err("Young's modulus must be positive".into());
        }
        if (0..nel).all(|e| self.kind(e) != ElementKind::Design) {
            return Err("no design elements".into());
        }
        Ok(())
    }

    /// One flag per dof: held by a support.
    pub fn fixed_dofs(&self) -> Vec<bool> {
        let mut fixed = vec![false; self.grid.ndof()];
        for s in &self.supports {
            for d in 0..3 {
                if s.dirs[d] {
                    fixed[3 * s.node + d] = true;
                }
            }
        }
        fixed
    }

    /// The force vector of each load case (zero on fixed dofs).
    pub fn force_vectors(&self, fixed: &[bool]) -> Vec<Vec<f64>> {
        self.load_cases
            .iter()
            .map(|c| {
                let mut f = vec![0.0; self.grid.ndof()];
                for nf in &c.forces {
                    for d in 0..3 {
                        f[3 * nf.node + d] += nf.f[d];
                    }
                }
                for (fi, &x) in f.iter_mut().zip(fixed) {
                    if x {
                        *fi = 0.0;
                    }
                }
                f
            })
            .collect()
    }

    /// Nodes whose position lies in the box `[lo, hi]` (mm, inclusive, with a small tolerance).
    pub fn nodes_in_box(&self, lo: [f64; 3], hi: [f64; 3]) -> Vec<usize> {
        let tol = 1e-9 * self.grid.h.max(1.0);
        (0..self.grid.nnodes())
            .filter(|&n| {
                let p = self.grid.node_pos(n);
                (0..3).all(|d| p[d] >= lo[d] - tol && p[d] <= hi[d] + tol)
            })
            .collect()
    }

    /// Elements whose centre lies in the box `[lo, hi]` (mm).
    pub fn elems_in_box(&self, lo: [f64; 3], hi: [f64; 3]) -> Vec<usize> {
        (0..self.grid.nel())
            .filter(|&e| {
                let p = self.grid.elem_center(e);
                (0..3).all(|d| p[d] >= lo[d] && p[d] <= hi[d])
            })
            .collect()
    }

    /// Holds every node in the box in the given directions.
    pub fn fix_box(&mut self, lo: [f64; 3], hi: [f64; 3], dirs: [bool; 3]) -> usize {
        let nodes = self.nodes_in_box(lo, hi);
        self.supports.extend(nodes.iter().map(|&node| NodeFix { node, dirs }));
        nodes.len()
    }

    /// Spreads a total force evenly over the nodes in the box, as a new load case.
    pub fn load_box(&mut self, lo: [f64; 3], hi: [f64; 3], total: [f64; 3]) -> usize {
        let nodes = self.nodes_in_box(lo, hi);
        let n = nodes.len().max(1) as f64;
        let forces = nodes.iter().map(|&node| NodeForce { node, f: [total[0] / n, total[1] / n, total[2] / n] }).collect();
        self.load_cases.push(LoadCase { weight: 1.0, forces });
        nodes.len()
    }

    /// Marks the elements in the box.
    pub fn set_kind_box(&mut self, lo: [f64; 3], hi: [f64; 3], kind: ElementKind) -> usize {
        if self.kinds.is_empty() {
            self.kinds = vec![ElementKind::Design; self.grid.nel()];
        }
        let elems = self.elems_in_box(lo, hi);
        for &e in &elems {
            self.kinds[e] = kind;
        }
        elems.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_round_trip_and_defaults() {
        let mut p = Problem::new(Grid::new(4, 2, 2, 2.5));
        p.fix_box([0.0, 0.0, 0.0], [0.0, 5.0, 5.0], [true, true, false]);
        p.load_box([10.0, 0.0, 0.0], [10.0, 5.0, 0.0], [0.0, 0.0, -100.0]);
        p.set_kind_box([0.0, 0.0, 0.0], [2.5, 5.0, 5.0], ElementKind::Solid);
        p.validate().unwrap();
        let text = serde_json::to_string(&p).unwrap();
        let back: Problem = serde_json::from_str(&text).unwrap();
        assert_eq!(back, p);
        // Minimal JSON: material, kinds, origin and weight fall back to their defaults.
        let minimal = r#"{"grid":{"nx":2,"ny":2,"nz":2,"h":1.0},"supports":[{"node":0,"dirs":[true,true,true]}],
                          "load_cases":[{"forces":[{"node":26,"f":[0,0,-1]}]}]}"#;
        let m: Problem = serde_json::from_str(minimal).unwrap();
        m.validate().unwrap();
        assert_eq!(m.load_cases[0].weight, 1.0);
        assert_eq!(m.material, Material::default());
        // Total load is spread evenly and supports zero their dofs.
        let f = p.force_vectors(&p.fixed_dofs());
        let total: f64 = f[0].iter().skip(2).step_by(3).sum();
        assert!((total + 100.0).abs() < 1e-9);
    }
}
