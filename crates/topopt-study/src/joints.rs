//! Every support joined to the part.
//!
//! Minimum compliance with several rigid supports leads the load into the supports nearest to
//! it: a support further away adds little stiffness, so the optimizer may leave it with no
//! material at all. For a part bolted down at several holes that is wrong: every hole must belong
//! to the part.
//!
//! The remedy is a helper load case per support that was left out: the loads, scaled down, with
//! the part held by that support alone. Its compliance is unbounded while nothing joins the
//! support to the loads, so the optimizer has to build that path; once it exists, the case only
//! weighs as much as its share of the load. [`crate::solve`] optimizes once without helper cases
//! and, when supports were left out, again from the start with them.

use topopt_core::{ElementKind, LoadCase, NodeFix, NodeForce, Problem};

/// The supports of a problem in pieces that do not touch: one per bolt hole, say.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportGroup {
    pub fixes: Vec<NodeFix>,
    /// Mean position of its nodes, mm.
    pub centre: [f64; 3],
}

impl SupportGroup {
    /// Whether the group alone holds the part: all three directions fixed at three or more
    /// nodes that are not on one line. A helper case needs that, or its system is singular.
    pub fn holds_alone(&self, problem: &Problem) -> bool {
        let g = problem.grid;
        let full: Vec<[f64; 3]> = self.fixes.iter().filter(|f| f.dirs == [true; 3]).map(|f| g.node_pos(f.node)).collect();
        if full.len() < 3 {
            return false;
        }
        let d = |a: [f64; 3], b: [f64; 3]| ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
        let a = full[0];
        let b = *full.iter().max_by(|p, q| d(a, **p).total_cmp(&d(a, **q))).unwrap();
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let l = d(a, b);
        if l < 0.5 * g.h {
            return false;
        }
        // Distance of each node from the line a–b.
        full.iter().any(|p| {
            let ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
            let c = [ap[1] * ab[2] - ap[2] * ab[1], ap[2] * ab[0] - ap[0] * ab[2], ap[0] * ab[1] - ap[1] * ab[0]];
            (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt() / l > 0.5 * g.h
        })
    }
}

/// The problem's supports grouped by contact: nodes next to each other (26 neighbours on the
/// node lattice) are one group. Groups come in the order of their first node.
pub fn support_groups(problem: &Problem) -> Vec<SupportGroup> {
    let g = problem.grid;
    // Every fix of each node (a node may be listed more than once).
    let mut at: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::with_capacity(problem.supports.len());
    for (k, f) in problem.supports.iter().enumerate() {
        at.entry(f.node).or_default().push(k);
    }
    let mut group = vec![usize::MAX; problem.supports.len()];
    let mut groups = Vec::new();
    for start in 0..problem.supports.len() {
        if group[start] != usize::MAX {
            continue;
        }
        let id = groups.len();
        let mut members = Vec::new();
        let mut stack = vec![start];
        group[start] = id;
        while let Some(k) = stack.pop() {
            members.push(k);
            let [i, j, l] = g.node_ijk(problem.supports[k].node);
            for dk in -1i64..=1 {
                for dj in -1i64..=1 {
                    for di in -1i64..=1 {
                        let (ni, nj, nl) = (i as i64 + di, j as i64 + dj, l as i64 + dk);
                        if ni < 0 || nj < 0 || nl < 0 || ni > g.nx as i64 || nj > g.ny as i64 || nl > g.nz as i64 {
                            continue;
                        }
                        let n = g.node(ni as usize, nj as usize, nl as usize);
                        for &q in at.get(&n).into_iter().flatten() {
                            if group[q] == usize::MAX {
                                group[q] = id;
                                stack.push(q);
                            }
                        }
                    }
                }
            }
        }
        members.sort_unstable();
        let fixes: Vec<NodeFix> = members.iter().map(|&k| problem.supports[k]).collect();
        let mut c = [0.0; 3];
        for f in &fixes {
            let p = g.node_pos(f.node);
            for d in 0..3 {
                c[d] += p[d] / fixes.len() as f64;
            }
        }
        groups.push(SupportGroup { fixes, centre: c });
    }
    groups
}

/// For each group: whether the material (elements at or above `iso`, and solid ones) joins it to
/// an element holding a load of the first `cases` load cases. Elements that share a node count as
/// joined, as they do in the analysis (on a coarse grid a diagonal strut is a chain of elements
/// touching at their edges).
pub fn joined(problem: &Problem, phys: &[f64], iso: f64, groups: &[SupportGroup], cases: usize) -> Vec<bool> {
    let g = problem.grid;
    let solid: Vec<bool> = (0..g.nel())
        .map(|e| match problem.kind(e) {
            ElementKind::Solid => true,
            ElementKind::Void => false,
            ElementKind::Design => phys[e] >= iso,
        })
        .collect();
    let mut loaded = vec![false; g.nnodes()];
    for c in problem.load_cases.iter().take(cases) {
        for f in &c.forces {
            loaded[f.node] = true;
        }
    }
    let touches = |e: usize, nodes: &[bool]| {
        let [i, j, k] = g.elem_ijk(e);
        g.elem_nodes(i, j, k).iter().any(|&n| nodes[n])
    };
    let mut reached = vec![false; g.nel()];
    let mut stack: Vec<usize> = (0..g.nel()).filter(|&e| solid[e] && touches(e, &loaded)).collect();
    for &e in &stack {
        reached[e] = true;
    }
    while let Some(e) = stack.pop() {
        let [i, j, k] = g.elem_ijk(e);
        let mut visit = |i: i64, j: i64, k: i64| {
            if i < 0 || j < 0 || k < 0 || i >= g.nx as i64 || j >= g.ny as i64 || k >= g.nz as i64 {
                return;
            }
            let f = g.elem(i as usize, j as usize, k as usize);
            if solid[f] && !reached[f] {
                reached[f] = true;
                stack.push(f);
            }
        };
        let (i, j, k) = (i as i64, j as i64, k as i64);
        for dk in -1..=1 {
            for dj in -1..=1 {
                for di in -1..=1 {
                    visit(i + di, j + dj, k + dk);
                }
            }
        }
    }
    groups
        .iter()
        .map(|grp| {
            let mut nodes = vec![false; g.nnodes()];
            for f in &grp.fixes {
                nodes[f.node] = true;
            }
            (0..g.nel()).any(|e| reached[e] && touches(e, &nodes))
        })
        .collect()
}

/// Weight of each helper case, relative to the real ones on the starting design
/// ([`topopt_core::LoadCase::relative`]): enough to keep a path to its support, little enough
/// not to reshape the part for it.
pub const HELPER_WEIGHT: f64 = 0.1;

/// `problem` with a helper load case for each of `helped` (indices into `groups`): the forces
/// of all its load cases together, divided by the number of groups, held by that group alone.
pub fn with_helper_cases(problem: &Problem, groups: &[SupportGroup], helped: &[usize]) -> Problem {
    let mut out = problem.clone();
    let share = 1.0 / groups.len().max(1) as f64;
    let mut total: std::collections::BTreeMap<usize, [f64; 3]> = std::collections::BTreeMap::new();
    for c in &problem.load_cases {
        for f in &c.forces {
            let t = total.entry(f.node).or_insert([0.0; 3]);
            for d in 0..3 {
                t[d] += f.f[d] * share;
            }
        }
    }
    let forces: Vec<NodeForce> = total.into_iter().map(|(node, f)| NodeForce { node, f }).collect();
    for &k in helped {
        out.load_cases.push(LoadCase { weight: HELPER_WEIGHT, forces: forces.clone(), supports: Some(groups[k].fixes.clone()), relative: true });
    }
    out
}

/// Where a group is, for messages: "(x, y, z) mm".
pub fn place(group: &SupportGroup) -> String {
    let c = group.centre;
    format!("({:.1}, {:.1}, {:.1}) mm", c[0], c[1], c[2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use topopt_core::Grid;

    #[test]
    fn groups_split_where_supports_do_not_touch() {
        let mut p = Problem::new(Grid::new(10, 4, 4, 1.0));
        p.fix_box([0.0, 0.0, 0.0], [0.0, 4.0, 4.0], [true; 3]);
        p.fix_box([10.0, 0.0, 0.0], [10.0, 4.0, 4.0], [true; 3]);
        p.fix_box([5.0, 0.0, 0.0], [5.0, 0.0, 0.0], [false, false, true]);
        let groups = support_groups(&p);
        assert_eq!(groups.len(), 3);
        assert_eq!(groups.iter().map(|g| g.fixes.len()).sum::<usize>(), p.supports.len());
        assert!(groups[0].holds_alone(&p));
        assert!(!groups.iter().find(|g| g.fixes.len() == 1).unwrap().holds_alone(&p));
    }
}
