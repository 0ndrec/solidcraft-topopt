//! Self-intersection test: spatial hashing of triangle boxes, then for every candidate pair that
//! shares no vertex, segment/triangle tests of each triangle's edges against the other. Two
//! triangles that cross (not coplanar) always have an edge of one passing through the other.
//! Coplanar overlaps are not detected; they do not arise from contouring or smoothing in practice.

use crate::mesh::{Mesh, cross, dot, sub};
use std::collections::HashMap;

fn aabb(m: &Mesh, t: usize) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for &v in &m.triangles[t] {
        let p = m.positions[v as usize];
        for d in 0..3 {
            lo[d] = lo[d].min(p[d]);
            hi[d] = hi[d].max(p[d]);
        }
    }
    (lo, hi)
}

/// Segment `p→q` crosses the interior of triangle `abc` (Möller–Trumbore, strict).
fn segment_hits_triangle(p: [f64; 3], q: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> bool {
    let dir = sub(q, p);
    let (e1, e2) = (sub(b, a), sub(c, a));
    let h = cross(dir, e2);
    let det = dot(e1, h);
    let scale = (dot(e1, e1) * dot(e2, e2)).sqrt() * dot(dir, dir).sqrt();
    if det.abs() <= 1e-12 * scale {
        return false; // parallel or coplanar
    }
    let inv = 1.0 / det;
    let s = sub(p, a);
    let u = inv * dot(s, h);
    let qv = cross(s, e1);
    let v = inv * dot(dir, qv);
    let t = inv * dot(e2, qv);
    let eps = 1e-9;
    u > eps && v > eps && u + v < 1.0 - eps && t > eps && t < 1.0 - eps
}

fn triangles_cross(m: &Mesh, s: usize, t: usize) -> bool {
    let ps = m.triangles[s].map(|i| m.positions[i as usize]);
    let pt = m.triangles[t].map(|i| m.positions[i as usize]);
    (0..3).any(|e| segment_hits_triangle(ps[e], ps[(e + 1) % 3], pt[0], pt[1], pt[2]))
        || (0..3).any(|e| segment_hits_triangle(pt[e], pt[(e + 1) % 3], ps[0], ps[1], ps[2]))
}

/// Number of crossing triangle pairs (pairs sharing a vertex are skipped).
pub fn count_self_intersections(m: &Mesh) -> usize {
    find_self_intersections(m).len()
}

/// The crossing triangle pairs.
pub fn find_self_intersections(m: &Mesh) -> Vec<(u32, u32)> {
    if m.triangles.is_empty() {
        return Vec::new();
    }
    let cell = (2.0 * m.mean_edge()).max(1e-12);
    let key = |p: [f64; 3]| [(p[0] / cell).floor() as i64, (p[1] / cell).floor() as i64, (p[2] / cell).floor() as i64];
    let boxes: Vec<([f64; 3], [f64; 3])> = (0..m.triangles.len()).map(|t| aabb(m, t)).collect();
    let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    for (t, (lo, hi)) in boxes.iter().enumerate() {
        let (a, b) = (key(*lo), key(*hi));
        for i in a[0]..=b[0] {
            for j in a[1]..=b[1] {
                for k in a[2]..=b[2] {
                    grid.entry([i, j, k]).or_default().push(t as u32);
                }
            }
        }
    }
    let mut found = Vec::new();
    for (cell_key, list) in &grid {
        for x in 0..list.len() {
            for y in x + 1..list.len() {
                let (s, t) = (list[x] as usize, list[y] as usize);
                let (bs, bt) = (&boxes[s], &boxes[t]);
                if (0..3).any(|d| bs.0[d] > bt.1[d] || bt.0[d] > bs.1[d]) {
                    continue;
                }
                // Test each pair once: in the cell holding the lower corner of the boxes' overlap.
                let lo = [bs.0[0].max(bt.0[0]), bs.0[1].max(bt.0[1]), bs.0[2].max(bt.0[2])];
                if key(lo) != *cell_key {
                    continue;
                }
                let (ts, tt) = (m.triangles[s], m.triangles[t]);
                if ts.iter().any(|v| tt.contains(v)) {
                    continue;
                }
                if triangles_cross(m, s, t) {
                    found.push((s as u32, t as u32));
                }
            }
        }
    }
    found
}
