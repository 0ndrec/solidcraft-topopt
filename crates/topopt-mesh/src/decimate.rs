//! Quadric-error edge-collapse decimation (Garland & Heckbert 1997) that keeps a closed
//! 2-manifold closed and manifold.
//!
//! A collapse `u ← v` is refused when it would:
//! - break the link condition (the edge's endpoints must share exactly the two opposite
//!   vertices of its two triangles), which is what keeps the surface a manifold;
//! - leave a vertex with fewer than three neighbours;
//! - turn an adjacent triangle's normal by more than `max_normal_turn_deg`, or make one that is
//!   needle-like and worse than before;
//! - move a pinned vertex: a pinned vertex stays where it is and absorbs its neighbour; two
//!   pinned vertices merge at one of their two positions.
//!
//! Quadrics are area-weighted; the error reported and limited is the mean squared distance to
//! the planes a vertex has absorbed.

use crate::mesh::{Mesh, cross, dot, norm, sub};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};

#[derive(Clone, Copy, Debug)]
pub struct DecimateOptions {
    /// Stop at this many triangles.
    pub target_triangles: usize,
    /// Stop when the cheapest collapse would move the surface by more than this (RMS distance).
    pub max_error: f64,
    pub max_normal_turn_deg: f64,
    /// Triangle quality 4√3·A / Σℓ² below which a collapse may not make a triangle worse.
    pub min_quality: f64,
}

impl Default for DecimateOptions {
    fn default() -> Self {
        DecimateOptions { target_triangles: 0, max_error: f64::INFINITY, max_normal_turn_deg: 60.0, min_quality: 0.15 }
    }
}

/// Symmetric 4×4 as a00 a01 a02 a03 a11 a12 a13 a22 a23 a33, plus the absorbed area.
type Quadric = [f64; 11];

fn plane_quadric(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Quadric {
    let n2 = cross(sub(b, a), sub(c, a));
    let len = norm(n2);
    if len == 0.0 {
        return [0.0; 11];
    }
    let n = [n2[0] / len, n2[1] / len, n2[2] / len];
    let d = -dot(n, a);
    let w = len / 2.0;
    [
        w * n[0] * n[0], w * n[0] * n[1], w * n[0] * n[2], w * n[0] * d,
        w * n[1] * n[1], w * n[1] * n[2], w * n[1] * d,
        w * n[2] * n[2], w * n[2] * d,
        w * d * d, w,
    ]
}

fn eval(q: &Quadric, p: [f64; 3]) -> f64 {
    let [x, y, z] = p;
    q[0] * x * x + 2.0 * q[1] * x * y + 2.0 * q[2] * x * z + 2.0 * q[3] * x
        + q[4] * y * y + 2.0 * q[5] * y * z + 2.0 * q[6] * y
        + q[7] * z * z + 2.0 * q[8] * z
        + q[9]
}

#[derive(PartialEq)]
struct Candidate {
    cost: f64,
    u: u32,
    v: u32,
    stamp_u: u32,
    stamp_v: u32,
    target: [f64; 3],
}
impl Eq for Candidate {}
impl Ord for Candidate {
    fn cmp(&self, o: &Self) -> Ordering {
        // Min-heap on cost.
        o.cost.partial_cmp(&self.cost).unwrap_or(Ordering::Equal)
    }
}
impl PartialOrd for Candidate {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

fn quality(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let area2 = norm(cross(sub(b, a), sub(c, a)));
    let l2 = dot(sub(b, a), sub(b, a)) + dot(sub(c, b), sub(c, b)) + dot(sub(a, c), sub(a, c));
    if l2 == 0.0 { 0.0 } else { 2.0 * 3f64.sqrt() * area2 / l2 }
}

struct State<'a> {
    pos: Vec<[f64; 3]>,
    tris: Vec<[u32; 3]>,
    alive: Vec<bool>,
    vtris: Vec<Vec<u32>>,
    q: Vec<Quadric>,
    stamp: Vec<u32>,
    dead: Vec<bool>,
    pinned: &'a [bool],
}

impl State<'_> {
    fn ring(&self, v: u32) -> HashSet<u32> {
        let mut r = HashSet::new();
        for &t in &self.vtris[v as usize] {
            if self.alive[t as usize] {
                for &w in &self.tris[t as usize] {
                    if w != v {
                        r.insert(w);
                    }
                }
            }
        }
        r
    }

    fn candidate(&self, u: u32, v: u32) -> Option<Candidate> {
        let (pu, pv) = (self.pinned[u as usize], self.pinned[v as usize]);
        // The kept vertex is the pinned one, if any. Two pinned vertices may merge into either
        // one's position: both lie on the original surface, so the result does too.
        let (keep, gone) = if pu && pv {
            let mut q = self.q[u as usize];
            for (a, b) in q.iter_mut().zip(&self.q[v as usize]) {
                *a += b;
            }
            if eval(&q, self.pos[u as usize]) <= eval(&q, self.pos[v as usize]) { (u, v) } else { (v, u) }
        } else if pv {
            (v, u)
        } else {
            (u, v)
        };
        let mut q = self.q[keep as usize];
        for (a, b) in q.iter_mut().zip(&self.q[gone as usize]) {
            *a += b;
        }
        let (a, b) = (self.pos[keep as usize], self.pos[gone as usize]);
        let mid = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0];
        let target = if pu || pv {
            a
        } else {
            // Minimize the quadric: solve A·p = −b.
            let m = [[q[0], q[1], q[2]], [q[1], q[4], q[5]], [q[2], q[5], q[7]]];
            let rhs = [-q[3], -q[6], -q[8]];
            let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
            let scale = (q[0] + q[4] + q[7]).powi(3);
            let edge = norm(sub(a, b));
            let solved = if det.abs() > 1e-9 * scale.max(1e-300) {
                let col = |c: usize| -> f64 {
                    let mut mm = m;
                    for r in 0..3 {
                        mm[r][c] = rhs[r];
                    }
                    mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1]) - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
                        + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0])
                };
                let p = [col(0) / det, col(1) / det, col(2) / det];
                if norm(sub(p, mid)) <= edge { Some(p) } else { None }
            } else {
                None
            };
            solved.unwrap_or_else(|| {
                *[a, b, mid].iter().min_by(|x, y| eval(&q, **x).partial_cmp(&eval(&q, **y)).unwrap()).unwrap()
            })
        };
        let cost = eval(&q, target).max(0.0) / q[10].max(1e-300);
        Some(Candidate { cost, u: keep, v: gone, stamp_u: self.stamp[keep as usize], stamp_v: self.stamp[gone as usize], target })
    }

    /// Whether collapsing `v` into `u` at `p` keeps the mesh valid (see module docs).
    fn allowed(&self, u: u32, v: u32, p: [f64; 3], cos_limit: f64, min_quality: f64) -> bool {
        let (ru, rv) = (self.ring(u), self.ring(v));
        let common: Vec<u32> = ru.intersection(&rv).copied().collect();
        if common.len() != 2 {
            return false;
        }
        for &w in &common {
            if self.ring(w).len() <= 3 {
                return false;
            }
        }
        if ru.len() + rv.len() - 2 < 3 + 2 {
            return false;
        }
        for (moving, other) in [(u, v), (v, u)] {
            for &t in &self.vtris[moving as usize] {
                if !self.alive[t as usize] {
                    continue;
                }
                let tri = self.tris[t as usize];
                if tri.contains(&other) {
                    continue; // removed by the collapse
                }
                let old = tri.map(|i| self.pos[i as usize]);
                let new = tri.map(|i| if i == moving { p } else { self.pos[i as usize] });
                let (n0, n1) = (cross(sub(old[1], old[0]), sub(old[2], old[0])), cross(sub(new[1], new[0]), sub(new[2], new[0])));
                let (l0, l1) = (norm(n0), norm(n1));
                if l1 <= 1e-14 * (l0 + 1e-300) || dot(n0, n1) < cos_limit * l0 * l1 {
                    return false;
                }
                let (q0, q1) = (quality(old[0], old[1], old[2]), quality(new[0], new[1], new[2]));
                if q1 < min_quality && q1 < q0 {
                    return false;
                }
            }
        }
        true
    }

    fn collapse(&mut self, u: u32, v: u32, p: [f64; 3]) -> usize {
        let mut removed = 0;
        let vt = std::mem::take(&mut self.vtris[v as usize]);
        for &t in &vt {
            if !self.alive[t as usize] {
                continue;
            }
            let tri = &mut self.tris[t as usize];
            if tri.contains(&u) {
                self.alive[t as usize] = false;
                removed += 1;
            } else {
                for w in tri.iter_mut() {
                    if *w == v {
                        *w = u;
                    }
                }
                self.vtris[u as usize].push(t);
            }
        }
        self.vtris[u as usize].retain(|&t| self.alive[t as usize]);
        self.pos[u as usize] = p;
        let qv = self.q[v as usize];
        for (a, b) in self.q[u as usize].iter_mut().zip(&qv) {
            *a += b;
        }
        self.dead[v as usize] = true;
        self.stamp[u as usize] += 1;
        self.stamp[v as usize] += 1;
        removed
    }
}

/// Collapses edges until the mesh has `target_triangles` or the next collapse would exceed
/// `max_error`. Returns the number of collapses. `pinned` has one flag per vertex.
pub fn decimate(mesh: &mut Mesh, pinned: &[bool], opts: &DecimateOptions) -> usize {
    assert_eq!(pinned.len(), mesh.positions.len());
    let nv = mesh.positions.len();
    let mut s = State {
        pos: mesh.positions.clone(),
        tris: mesh.triangles.clone(),
        alive: vec![true; mesh.triangles.len()],
        vtris: vec![Vec::new(); nv],
        q: vec![[0.0; 11]; nv],
        stamp: vec![0; nv],
        dead: vec![false; nv],
        pinned,
    };
    for (t, tri) in s.tris.iter().enumerate() {
        let pq = plane_quadric(s.pos[tri[0] as usize], s.pos[tri[1] as usize], s.pos[tri[2] as usize]);
        for &v in tri {
            s.vtris[v as usize].push(t as u32);
            for (a, b) in s.q[v as usize].iter_mut().zip(&pq) {
                *a += b;
            }
        }
    }
    let mut heap = BinaryHeap::new();
    let mut seen = HashSet::new();
    for tri in &s.tris {
        for e in 0..3 {
            let (a, b) = (tri[e].min(tri[(e + 1) % 3]), tri[e].max(tri[(e + 1) % 3]));
            if seen.insert((a, b)) {
                if let Some(c) = s.candidate(a, b) {
                    heap.push(c);
                }
            }
        }
    }
    let cos_limit = opts.max_normal_turn_deg.to_radians().cos();
    let max_cost = opts.max_error * opts.max_error;
    let mut alive_tris = s.tris.len();
    let mut collapses = 0;
    while alive_tris > opts.target_triangles.max(4) {
        let Some(c) = heap.pop() else { break };
        let (u, v) = (c.u as usize, c.v as usize);
        if s.dead[u] || s.dead[v] || s.stamp[u] != c.stamp_u || s.stamp[v] != c.stamp_v {
            continue;
        }
        if c.cost > max_cost {
            break;
        }
        if !s.allowed(c.u, c.v, c.target, cos_limit, opts.min_quality) {
            continue;
        }
        alive_tris -= s.collapse(c.u, c.v, c.target);
        collapses += 1;
        for w in s.ring(c.u) {
            if let Some(nc) = s.candidate(c.u.min(w), c.u.max(w)) {
                heap.push(nc);
            }
        }
    }
    mesh.triangles = s.tris.iter().zip(&s.alive).filter(|(_, a)| **a).map(|(t, _)| *t).collect();
    mesh.positions = s.pos;
    mesh.compact();
    collapses
}
