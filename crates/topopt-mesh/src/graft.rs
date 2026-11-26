//! Kept faces exactly as the CAD model has them: the result's triangles on chosen faces of the
//! original surface are replaced by the original surface's own triangles there (its tessellation,
//! from `model.mesh`), stitched watertight to the rest of the result.
//!
//! For each group of connected kept faces (a **patch**):
//! 1. **cut**: the result's vertices that lie on the patch (within a quarter of a sample spacing,
//!    on its outer side, not on its outline) are covered; every triangle with a covered vertex or
//!    centroid is removed. What is left has one hole per loop of the patch's outline;
//! 2. **insert** the patch's triangles;
//! 3. **weld**: each hole loop is matched to the nearest outline loop. Hole vertices on the outline
//!    (the crease vertices the pinning put on the sharp edges) become outline vertices: an outline
//!    segment is split where one lands, a vertex within a hair is merged into the outline's;
//! 4. **bridge**: between two welded vertices, the hole runs along the outline (its triangle's edge
//!    is split at the outline vertices in between) or away from it, where a ladder of triangles
//!    joins the two chains: each hole vertex gets a foot on the outline (the outline is split
//!    there), and neighbouring rungs are joined by two triangles. A gap that is nearly flat (most are:
//!    it lies in the neighbouring face) is ear-clipped in its plane instead, which cannot lay
//!    triangles over each other where the outline turns a corner;
//! 5. **check**: the whole mesh is validated (closed, manifold, oriented, no crossings or folds).
//!
//! [`graft_kept_faces`] grafts all patches at once and, if that fails, one at a time, keeping
//! those that pass; a patch that fails stays as the pinned extraction made it.

use crate::mesh::{Mesh, cross, dot, norm, sub, vertex_normals};
use parry3d_f64::math::Vector;
use parry3d_f64::query::{PointQuery, PointQueryWithLocation};
use parry3d_f64::shape::{Polyline, TriMesh};
use std::collections::HashMap;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GraftReport {
    /// Groups of connected kept faces asked for, and how many were grafted.
    pub patches: usize,
    pub grafted: usize,
    /// Result triangles removed from the kept faces, and original triangles put in their place.
    pub removed_triangles: usize,
    pub patch_triangles: usize,
    /// Result vertices welded onto the kept faces' outlines.
    pub welded: usize,
    /// Triangles added between the result and an outline where they do not meet.
    pub bridged: usize,
    /// Why patches were not grafted (one line each).
    pub failures: Vec<String>,
}

/// The original surface, welded, outward, with what the graft asks of it.
struct Original {
    positions: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    normals: Vec<[f64; 3]>,
    shape: TriMesh,
}

impl Original {
    fn new(positions: &[[f64; 3]], triangles: &[[u32; 3]]) -> Result<Self, String> {
        // Weld by position: CAD meshes may duplicate vertices along face boundaries.
        let mut weld: HashMap<[u64; 3], u32> = HashMap::new();
        let mut welded = Vec::new();
        let ids: Vec<u32> = positions
            .iter()
            .map(|p| {
                *weld.entry(p.map(|c| (c + 0.0).to_bits())).or_insert_with(|| {
                    welded.push(*p);
                    welded.len() as u32 - 1
                })
            })
            .collect();
        let mut tris: Vec<[u32; 3]> = triangles.iter().map(|t| t.map(|v| ids[v as usize])).collect();
        let volume: f64 = tris.iter().map(|t| dot(welded[t[0] as usize], cross(welded[t[1] as usize], welded[t[2] as usize])) / 6.0).sum();
        if volume < 0.0 {
            tris.iter_mut().for_each(|t| t.swap(1, 2));
        }
        let normals = tris
            .iter()
            .map(|t| {
                let n = cross(sub(welded[t[1] as usize], welded[t[0] as usize]), sub(welded[t[2] as usize], welded[t[0] as usize]));
                let l = norm(n);
                if l > 0.0 { n.map(|c| c / l) } else { [0.0; 3] }
            })
            .collect();
        let shape = TriMesh::new(welded.iter().map(|&p| v3(p)).collect(), tris.clone()).map_err(|e| format!("surface mesh: {e:?}"))?;
        Ok(Original { positions: welded, triangles: tris, normals, shape })
    }

    /// Patch id of every kept triangle (kept triangles joined by edges), `u32::MAX` for the rest.
    fn patches(&self, keep: &[bool]) -> (Vec<u32>, usize) {
        let mut parent: Vec<u32> = (0..self.triangles.len() as u32).collect();
        fn find(p: &mut [u32], mut x: u32) -> u32 {
            while p[x as usize] != x {
                p[x as usize] = p[p[x as usize] as usize];
                x = p[x as usize];
            }
            x
        }
        let mut first: HashMap<(u32, u32), u32> = HashMap::new();
        for (t, tri) in self.triangles.iter().enumerate() {
            if !keep[t] {
                continue;
            }
            for e in 0..3 {
                let (a, b) = (tri[e], tri[(e + 1) % 3]);
                match first.get(&(a.min(b), a.max(b))) {
                    Some(&o) => {
                        let (ra, rb) = (find(&mut parent, o), find(&mut parent, t as u32));
                        parent[ra as usize] = rb;
                    }
                    None => {
                        first.insert((a.min(b), a.max(b)), t as u32);
                    }
                }
            }
        }
        let mut ids: HashMap<u32, u32> = HashMap::new();
        let mut patch = vec![u32::MAX; self.triangles.len()];
        for t in 0..self.triangles.len() {
            if keep[t] {
                let r = find(&mut parent, t as u32);
                let n = ids.len() as u32;
                patch[t] = *ids.entry(r).or_insert(n);
            }
        }
        (patch, ids.len())
    }
}

/// Replaces the result's triangles on the kept triangles of the original surface (`keep`: one flag
/// per triangle of `surface`) by those triangles. `h`: the sample spacing the result was made at.
/// Never fails: a patch that cannot be grafted is left as it was, and the report says why.
pub fn graft_kept_faces(mesh: &Mesh, surface: (&[[f64; 3]], &[[u32; 3]]), keep: &[bool], h: f64, check_self_intersections: bool) -> (Mesh, GraftReport) {
    let mut report = GraftReport::default();
    if keep.len() != surface.1.len() || !keep.iter().any(|&k| k) {
        return (mesh.clone(), report);
    }
    let original = match Original::new(surface.0, surface.1) {
        Ok(o) => o,
        Err(e) => {
            report.failures.push(e);
            return (mesh.clone(), report);
        }
    };
    let (patch, count) = original.patches(keep);
    report.patches = count;
    let all: Vec<bool> = patch.iter().map(|&p| p != u32::MAX).collect();
    match graft_either(mesh, &original, &all, h, check_self_intersections) {
        Ok((m, s)) => {
            report.grafted = count;
            s.add_to(&mut report);
            return (m, report);
        }
        Err(e) if count == 1 => {
            report.failures.push(e);
            return (mesh.clone(), report);
        }
        Err(_) => {}
    }
    let mut current = mesh.clone();
    for p in 0..count as u32 {
        let one: Vec<bool> = patch.iter().map(|&q| q == p).collect();
        match graft_either(&current, &original, &one, h, check_self_intersections) {
            Ok((m, s)) => {
                current = m;
                report.grafted += 1;
                s.add_to(&mut report);
            }
            Err(e) => report.failures.push(format!("patch {}: {e}", p + 1)),
        }
    }
    (current, report)
}

/// Flat gaps by ear clipping (ladders elsewhere); if that fails validation, ladders everywhere.
fn graft_either(mesh: &Mesh, original: &Original, keep: &[bool], h: f64, check: bool) -> Result<(Mesh, Stats), String> {
    graft(mesh, original, keep, h, check, true).or_else(|e| graft(mesh, original, keep, h, check, false).map_err(|_| e))
}

#[derive(Default)]
struct Stats {
    removed: usize,
    inserted: usize,
    welded: usize,
    bridged: usize,
}

impl Stats {
    fn add_to(&self, r: &mut GraftReport) {
        r.removed_triangles += self.removed;
        r.patch_triangles += self.inserted;
        r.welded += self.welded;
        r.bridged += self.bridged;
    }
}

/// Boundary loops of a triangle set: each loop follows the boundary half-edges (those without a
/// twin) in their direction. Fails where the boundary touches itself at a vertex.
fn boundary_loops(triangles: &[[u32; 3]]) -> Result<Vec<Vec<u32>>, String> {
    let mut half: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    for t in triangles {
        for e in 0..3 {
            half.insert((t[e], t[(e + 1) % 3]));
        }
    }
    let mut next: HashMap<u32, u32> = HashMap::new();
    for &(a, b) in &half {
        if !half.contains(&(b, a)) && next.insert(a, b).is_some() {
            return Err("the outline touches itself at a vertex".into());
        }
    }
    let mut starts: Vec<u32> = next.keys().copied().collect();
    starts.sort_unstable();
    let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let mut loops = Vec::new();
    for s in starts {
        if seen.contains(&s) {
            continue;
        }
        let mut lp = vec![s];
        seen.insert(s);
        let mut v = next[&s];
        while v != s {
            if !seen.insert(v) {
                return Err("an outline does not close".into());
            }
            lp.push(v);
            v = *next.get(&v).ok_or("an outline does not close")?;
        }
        loops.push(lp);
    }
    Ok(loops)
}

/// Closest point of segment `a`–`b` to `p`: (parameter in [0, 1], distance).
fn closest_on_segment(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> (f64, f64) {
    let ab = sub(b, a);
    let l2 = dot(ab, ab);
    let t = if l2 > 0.0 { (dot(sub(p, a), ab) / l2).clamp(0.0, 1.0) } else { 0.0 };
    let q = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
    (t, norm(sub(p, q)))
}

/// A loop of the outline as a closed polyline: its vertices, and the arc length at each.
struct Outline {
    vertices: Vec<u32>,
    at: Vec<f64>,
    length: f64,
}

impl Outline {
    fn new(vertices: Vec<u32>, pos: &[[f64; 3]]) -> Self {
        let n = vertices.len();
        let mut at = Vec::with_capacity(n);
        let mut s = 0.0;
        for i in 0..n {
            at.push(s);
            s += norm(sub(pos[vertices[(i + 1) % n] as usize], pos[vertices[i] as usize]));
        }
        Outline { vertices, at, length: s }
    }

    /// Nearest point to `p`: (segment, parameter, distance, arc length).
    fn nearest(&self, p: [f64; 3], pos: &[[f64; 3]]) -> (usize, f64, f64, f64) {
        let n = self.vertices.len();
        let mut best = (0, 0.0, f64::INFINITY, 0.0);
        for i in 0..n {
            let (a, b) = (pos[self.vertices[i] as usize], pos[self.vertices[(i + 1) % n] as usize]);
            let (t, d) = closest_on_segment(p, a, b);
            if d < best.2 {
                best = (i, t, d, self.at[i] + t * norm(sub(b, a)));
            }
        }
        best
    }
}

/// The mesh being stitched: triangles with a half-edge index.
struct Stitch {
    positions: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
    half: HashMap<(u32, u32), usize>,
}

impl Stitch {
    fn index(&mut self) {
        self.half.clear();
        for (t, tri) in self.triangles.iter().enumerate() {
            for e in 0..3 {
                self.half.insert((tri[e], tri[(e + 1) % 3]), t);
            }
        }
    }

    /// Splits the triangle that has half-edge `u`→`v` at `points` (in order from `u` to `v`): a fan
    /// from its third corner.
    fn split(&mut self, u: u32, v: u32, points: &[u32]) -> Result<(), String> {
        if points.is_empty() {
            return Ok(());
        }
        let t = *self.half.get(&(u, v)).ok_or("an outline edge has no triangle")?;
        let tri = self.triangles[t];
        let x = (0..3).map(|e| tri[e]).find(|&c| c != u && c != v).ok_or("a degenerate triangle on the outline")?;
        for e in 0..3 {
            self.half.remove(&(tri[e], tri[(e + 1) % 3]));
        }
        let mut chain = Vec::with_capacity(points.len() + 2);
        chain.push(u);
        chain.extend_from_slice(points);
        chain.push(v);
        for k in 0..chain.len() - 1 {
            let nt = [chain[k], chain[k + 1], x];
            let id = if k == 0 {
                self.triangles[t] = nt;
                t
            } else {
                self.triangles.push(nt);
                self.triangles.len() - 1
            };
            for e in 0..3 {
                self.half.insert((nt[e], nt[(e + 1) % 3]), id);
            }
        }
        Ok(())
    }

    /// Fills the gap between chains `a` and `b` (as for [`Stitch::ladder`]) when the polygon they
    /// make is flat (within `flat`): ear clipping in its plane, the best-shaped ear first, which
    /// never lays triangles over each other. `None` when it is not flat or not simple.
    fn ears(&mut self, a: &[u32], b: &[u32], flat: f64) -> Option<usize> {
        if a.first() != b.first() || a.last() != b.last() || a.len() + b.len() < 5 {
            return None;
        }
        // Fill order: along b, then back along a.
        let mut poly: Vec<u32> = b.to_vec();
        poly.extend(a[1..a.len() - 1].iter().rev());
        let pos: Vec<[f64; 3]> = poly.iter().map(|&v| self.positions[v as usize]).collect();
        let n = poly.len();
        let mut nn = [0.0; 3];
        for i in 0..n {
            let (p, q) = (pos[i], pos[(i + 1) % n]);
            nn[0] += (p[1] - q[1]) * (p[2] + q[2]);
            nn[1] += (p[2] - q[2]) * (p[0] + q[0]);
            nn[2] += (p[0] - q[0]) * (p[1] + q[1]);
        }
        // The plane: least spread of the points (a long thin gap whose points wander a little off
        // it tilts the Newell normal); the Newell normal only says which side is out.
        let Some(mut nn2) = least_spread(&pos) else { return None };
        if dot(nn2, nn) < 0.0 {
            nn2 = nn2.map(|c| -c);
        }
        if norm(nn) <= 0.0 {
            return None;
        }
        let nn = nn2;
        let d0 = pos.iter().map(|&p| dot(p, nn)).sum::<f64>() / n as f64;
        if pos.iter().any(|&p| (dot(p, nn) - d0).abs() > flat) {
            return None;
        }
        let u = {
            let t = if nn[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
            let c = cross(nn, t);
            let lc = norm(c);
            c.map(|x| x / lc)
        };
        let w = cross(nn, u);
        let q: Vec<[f64; 2]> = pos.iter().map(|&p| [dot(p, u), dot(p, w)]).collect();
        let cross2 = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        let scale = q.iter().map(|p| p[0].abs().max(p[1].abs())).fold(0.0, f64::max).max(1e-9);
        let eps = 1e-12 * scale * scale;
        let mut idx: Vec<usize> = (0..n).collect();
        let mut tris = Vec::with_capacity(n - 2);
        while idx.len() > 3 {
            let m = idx.len();
            let mut best: Option<(usize, f64)> = None;
            for k in 0..m {
                let (i0, i1, i2) = (idx[(k + m - 1) % m], idx[k], idx[(k + 1) % m]);
                let (p0, p1, p2) = (q[i0], q[i1], q[i2]);
                let area = cross2(p0, p1, p2);
                if area <= eps {
                    continue;
                }
                let blocked = idx.iter().any(|&j| {
                    if j == i0 || j == i1 || j == i2 || q[j] == p0 || q[j] == p1 || q[j] == p2 {
                        return false;
                    }
                    let pj = q[j];
                    cross2(p0, p1, pj) >= -eps && cross2(p1, p2, pj) >= -eps && cross2(p2, p0, pj) >= -eps
                });
                if blocked {
                    continue;
                }
                // Shape: area over the squared longest edge.
                let e = |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2);
                let quality = area / e(p0, p1).max(e(p1, p2)).max(e(p2, p0));
                if best.is_none_or(|b| quality > b.1) {
                    best = Some((k, quality));
                }
            }
            let (k, _) = best?;
            let m = idx.len();
            tris.push([poly[idx[(k + m - 1) % m]], poly[idx[k]], poly[idx[(k + 1) % m]]]);
            idx.remove(k);
        }
        if cross2(q[idx[0]], q[idx[1]], q[idx[2]]) <= eps {
            return None;
        }
        tris.push([poly[idx[0]], poly[idx[1]], poly[idx[2]]]);
        let added = tris.len();
        self.triangles.extend(tris);
        Some(added)
    }

    /// Fills the gap between chain `a` (the result's side, along its boundary) and chain `b` (the
    /// outline's, against its boundary), which start together and end together, as a ladder: each
    /// vertex of `a` gets its foot on `b` (its nearest point, never behind the previous foot; an
    /// outline vertex within `snap` of it is used as it is), the outline's triangles are split
    /// there, and the rungs are joined by pairs of triangles. Returns the triangles added.
    fn ladder(&mut self, a: &[u32], b: &[u32], snap: f64) -> Result<usize, String> {
        let pos = |s: &Self, v: u32| s.positions[v as usize];
        let mut sb = vec![0.0];
        for k in 1..b.len() {
            sb.push(sb[k - 1] + norm(sub(pos(self, b[k]), pos(self, b[k - 1]))));
        }
        let total = *sb.last().unwrap();
        let mut ta = vec![0.0; a.len()];
        *ta.last_mut().unwrap() = total;
        let mut floor = 0.0f64;
        for i in 1..a.len() - 1 {
            let p = pos(self, a[i]);
            let (mut at, mut best) = (floor, f64::INFINITY);
            for k in 0..b.len() - 1 {
                let (t, d) = closest_on_segment(p, pos(self, b[k]), pos(self, b[k + 1]));
                if d < best {
                    best = d;
                    at = sb[k] + t * (sb[k + 1] - sb[k]);
                }
            }
            floor = at.max(floor).min(total);
            ta[i] = floor;
        }
        // The feet: the outline's own vertices, and new ones on its segments.
        let mut rungs: Vec<(f64, u32)> = b.iter().enumerate().map(|(k, &v)| (sb[k], v)).collect();
        let mut on_segment: HashMap<usize, Vec<(f64, u32)>> = HashMap::new();
        for &t in &ta[1..a.len() - 1] {
            if rungs.iter().any(|r| (r.0 - t).abs() <= snap) {
                continue;
            }
            let Some(k) = (0..b.len() - 1).find(|&k| sb[k] < t && t < sb[k + 1]) else { continue };
            let f = (t - sb[k]) / (sb[k + 1] - sb[k]);
            let (p, q) = (pos(self, b[k]), pos(self, b[k + 1]));
            let id = self.positions.len() as u32;
            self.positions.push([p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1]), p[2] + f * (q[2] - p[2])]);
            rungs.push((t, id));
            on_segment.entry(k).or_default().push((t, id));
        }
        // The chain runs against the outline: segment k is the outline's half-edge b[k+1] → b[k].
        let mut segs: Vec<usize> = on_segment.keys().copied().collect();
        segs.sort_unstable();
        for k in segs {
            let mut list = on_segment.remove(&k).unwrap();
            list.sort_by(|x, y| y.0.total_cmp(&x.0));
            let ids: Vec<u32> = list.iter().map(|x| x.1).collect();
            self.split(b[k + 1], b[k], &ids)?;
        }
        rungs.sort_by(|x, y| x.0.total_cmp(&y.0));
        let (tb, bb): (Vec<f64>, Vec<u32>) = rungs.into_iter().unzip();
        let (mut i, mut j) = (0, 0);
        let mut new = Vec::new();
        while i + 1 < a.len() || j + 1 < bb.len() {
            let advance_a = if i + 1 == a.len() {
                false
            } else if j + 1 == bb.len() {
                true
            } else {
                ta[i + 1] <= tb[j + 1] + 0.5 * snap
            };
            let t = if advance_a {
                i += 1;
                [a[i], a[i - 1], bb[j]]
            } else {
                j += 1;
                [bb[j - 1], bb[j], a[i]]
            };
            if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
                new.push(t);
            }
        }
        let added = new.len();
        self.triangles.extend(new);
        Ok(added)
    }
}

/// One graft of the kept triangles `keep` (all of one or more patches).
fn graft(mesh: &Mesh, original: &Original, keep: &[bool], h: f64, check: bool, ears: bool) -> Result<(Mesh, Stats), String> {
    let mut stats = Stats::default();
    let opos = &original.positions;
    // The outline of the kept triangles, for "is this point on the outline".
    let kept_tris: Vec<[u32; 3]> = (0..original.triangles.len()).filter(|&t| keep[t]).map(|t| original.triangles[t]).collect();
    let outline_loops = boundary_loops(&kept_tris)?;
    if outline_loops.is_empty() {
        return Err("the kept faces have no outline".into());
    }
    let mut segments = Vec::new();
    for lp in &outline_loops {
        for i in 0..lp.len() {
            segments.push([lp[i], lp[(i + 1) % lp.len()]]);
        }
    }
    let outline_shape = Polyline::new(opos.iter().map(|&p| v3(p)).collect(), Some(segments));

    // 1. Cut.
    let near = 0.25 * h;
    let off_outline = 1e-3 * h;
    let covered = |p: [f64; 3], n: [f64; 3]| -> bool {
        let (proj, (tri, _)) = original.shape.project_local_point_and_get_location(v3(p), false);
        let q = [proj.point.x, proj.point.y, proj.point.z];
        let t = tri as usize;
        let d = norm(sub(p, q));
        // On the face: covered whatever its normal (next to a rounded edge it points along the
        // rounding). Near it: only on its outer side (not the far side of a thin wall).
        keep[t]
            && d <= near
            && (d <= 1e-6 * h || dot(n, original.normals[t]) > 0.3)
            && outline_shape.distance_to_local_point(v3(q), false) > off_outline
    };
    let normals = vertex_normals(mesh);
    let cover: Vec<bool> = (0..mesh.positions.len()).map(|v| covered(mesh.positions[v], normals[v])).collect();
    let mut rest = Vec::with_capacity(mesh.triangles.len());
    for (t, tri) in mesh.triangles.iter().enumerate() {
        let gone = tri.iter().any(|&v| cover[v as usize]) || {
            let [a, b, c] = tri.map(|v| mesh.positions[v as usize]);
            let centroid = [(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0];
            let n = mesh.face_normal(t);
            let l = norm(n);
            l > 0.0 && covered(centroid, n.map(|c| c / l))
        };
        if gone {
            stats.removed += 1;
        } else {
            rest.push(*tri);
        }
    }
    if stats.removed == 0 {
        return Err("the result does not lie on the kept faces".into());
    }
    // Spikes: a triangle with two edges on the cut, whose corner between them is off the outline,
    // juts into the hole and would lie across the bridge to the outline. Cut it too (a few rounds:
    // a large triangle left by decimation can uncover another).
    for _ in 0..8 {
        let mut half: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::with_capacity(rest.len() * 3);
        for t in &rest {
            for e in 0..3 {
                half.insert((t[e], t[(e + 1) % 3]));
            }
        }
        let before = rest.len();
        rest.retain(|t| {
            let open: Vec<bool> = (0..3).map(|e| !half.contains(&(t[(e + 1) % 3], t[e]))).collect();
            // Corner c lies between edges c-1 → c and c → c+1.
            !(0..3).any(|c| {
                open[c] && open[(c + 2) % 3] && {
                    let p = mesh.positions[t[c] as usize];
                    outline_shape.distance_to_local_point(v3(p), false) > 0.02 * h
                }
            })
        });
        stats.removed += before - rest.len();
        if rest.len() == before {
            break;
        }
    }
    let holes = boundary_loops(&rest)?;

    // 2. Insert: the original's vertices go after the result's.
    let base = mesh.positions.len() as u32;
    let mut s = Stitch { positions: mesh.positions.clone(), triangles: rest, half: HashMap::new() };
    s.positions.extend_from_slice(opos);
    for t in &kept_tris {
        s.triangles.push(t.map(|v| v + base));
    }
    stats.inserted = kept_tris.len();
    let mut outlines: Vec<Outline> = outline_loops.iter().map(|lp| Outline::new(lp.iter().map(|&v| v + base).collect(), &s.positions)).collect();

    // 3. Match each hole to the nearest outline loop, one to one.
    if holes.len() != outlines.len() {
        return Err(format!("{} holes cut for {} outline loops", holes.len(), outlines.len()));
    }
    let mut owner = vec![usize::MAX; outlines.len()];
    for (k, hole) in holes.iter().enumerate() {
        let (best, mean) = (0..outlines.len())
            .map(|o| (o, hole.iter().map(|&v| outlines[o].nearest(s.positions[v as usize], &s.positions).2).sum::<f64>() / hole.len() as f64))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap();
        if mean > 1.5 * h {
            return Err(format!("a cut is {:.2} samples from the outline", mean / h));
        }
        if owner[best] != usize::MAX {
            return Err("two cuts meet one outline loop".into());
        }
        owner[best] = k;
    }

    let weld = 0.02 * h;
    let mut remap: Vec<u32> = (0..s.positions.len() as u32).collect();
    for (o, &k) in owner.iter().enumerate() {
        let hole = &holes[k];
        let m = hole.len();
        // Weld candidates: (hole index, segment, parameter, distance, arc length).
        let mut cand: Vec<Option<(usize, f64, f64, f64)>> = Vec::with_capacity(m);
        {
            let out = &outlines[o];
            let n = out.vertices.len();
            for &v in hole {
                let (seg, mut t, d, mut at) = out.nearest(s.positions[v as usize], &s.positions);
                if d > weld {
                    cand.push(None);
                    continue;
                }
                let len = out.at.get(seg + 1).copied().unwrap_or(out.length) - out.at[seg];
                if t * len < weld {
                    (t, at) = (0.0, out.at[seg]);
                } else if (1.0 - t) * len < weld {
                    (t, at) = (1.0, out.at[(seg + 1) % n]);
                }
                cand.push(Some((seg, t, d, at)));
            }
        }
        // Keep a subsequence that runs back along the outline once (the hole runs against it).
        let length = outlines[o].length;
        let start = (0..m).filter(|&i| cand[i].is_some()).min_by(|&a, &b| cand[a].unwrap().2.total_cmp(&cand[b].unwrap().2));
        let mut anchors: Vec<usize> = Vec::new();
        if let Some(i0) = start {
            anchors.push(i0);
            let mut prev = cand[i0].unwrap().3;
            let mut total = 0.0;
            for step in 1..m {
                let i = (i0 + step) % m;
                if let Some((_, _, _, at)) = cand[i] {
                    let back = (prev - at).rem_euclid(length);
                    if back > 1e-9 * length && back < 0.5 * length {
                        anchors.push(i);
                        total += back;
                        prev = at;
                    }
                }
            }
            if anchors.len() > 1 {
                total += (prev - cand[i0].unwrap().3).rem_euclid(length);
                if (total - length).abs() > 1e-6 * length {
                    // Not once around: keep the best one only.
                    anchors.truncate(1);
                }
            }
        }
        // New outline vertices where anchors land inside segments, then the anchors become them.
        let mut inserts: HashMap<usize, Vec<(f64, usize)>> = HashMap::new();
        let mut target = vec![u32::MAX; m];
        for &i in &anchors {
            let (seg, t, _, _) = cand[i].unwrap();
            let out = &outlines[o];
            let n = out.vertices.len();
            if t <= 0.0 {
                target[i] = out.vertices[seg];
            } else if t >= 1.0 {
                target[i] = out.vertices[(seg + 1) % n];
            } else {
                inserts.entry(seg).or_default().push((t, i));
            }
        }
        s.index();
        let mut new_vertices: Vec<u32> = Vec::new();
        let mut segs: Vec<usize> = inserts.keys().copied().collect();
        segs.sort_unstable();
        let mut after: HashMap<usize, Vec<u32>> = HashMap::new();
        for seg in segs {
            let list = inserts.get_mut(&seg).unwrap();
            list.sort_by(|a, b| a.0.total_cmp(&b.0));
            let out = &outlines[o];
            let n = out.vertices.len();
            let (u, v) = (out.vertices[seg], out.vertices[(seg + 1) % n]);
            let (pu, pv) = (s.positions[u as usize], s.positions[v as usize]);
            let mut ids = Vec::new();
            for &(t, i) in list.iter() {
                let id = s.positions.len() as u32;
                s.positions.push([pu[0] + t * (pv[0] - pu[0]), pu[1] + t * (pv[1] - pu[1]), pu[2] + t * (pv[2] - pu[2])]);
                ids.push(id);
                target[i] = id;
                new_vertices.push(id);
            }
            s.split(u, v, &ids)?;
            after.insert(seg, ids);
        }
        if !after.is_empty() {
            let out = &outlines[o];
            let mut vs = Vec::with_capacity(out.vertices.len() + new_vertices.len());
            for (seg, &v) in out.vertices.iter().enumerate() {
                vs.push(v);
                if let Some(ids) = after.get(&seg) {
                    vs.extend_from_slice(ids);
                }
            }
            outlines[o] = Outline::new(vs, &s.positions);
        }
        // Vertices added since (new outline vertices here, ladder feet of earlier loops) map to
        // themselves.
        while remap.len() < s.positions.len() {
            remap.push(remap.len() as u32);
        }
        for &i in &anchors {
            remap[hole[i] as usize] = target[i];
        }
        stats.welded += anchors.len();
        for t in s.triangles.iter_mut() {
            *t = t.map(|v| remap[v as usize]);
        }
        let hole: Vec<u32> = hole.iter().map(|&v| remap[v as usize]).collect();
        s.index();

        // 4. Bridge between consecutive anchors (or the whole loop when there are none).
        let out = &outlines[o];
        let n = out.vertices.len();
        let place: HashMap<u32, usize> = out.vertices.iter().enumerate().map(|(i, &v)| (v, i)).collect();
        let back = |from: usize, to: usize, whole: bool| -> Vec<u32> {
            let mut c = vec![out.vertices[from]];
            let mut i = from;
            loop {
                i = (i + n - 1) % n;
                c.push(out.vertices[i]);
                if i == to && (!whole || c.len() > n) {
                    break;
                }
            }
            c
        };
        if anchors.is_empty() {
            // No contact: zip the whole loop, from the hole vertex nearest to an outline vertex.
            let (hi, oi) = (0..m)
                .flat_map(|i| (0..n).map(move |j| (i, j)))
                .min_by(|a, b| {
                    let da = norm(sub(s.positions[hole[a.0] as usize], s.positions[out.vertices[a.1] as usize]));
                    let db = norm(sub(s.positions[hole[b.0] as usize], s.positions[out.vertices[b.1] as usize]));
                    da.total_cmp(&db)
                })
                .unwrap();
            let a: Vec<u32> = (0..=m).map(|k| hole[(hi + k) % m]).collect();
            let b = back(oi, oi, true);
            stats.bridged += s.ladder(&a, &b, 0.02 * h)?;
            continue;
        }
        let mut sorted = anchors.clone();
        sorted.sort_unstable();
        let rotate = sorted.iter().position(|&i| i == anchors[0]).unwrap();
        sorted.rotate_left(rotate);
        let k = sorted.len();
        let mut jobs = Vec::new();
        for (idx, &i) in sorted.iter().enumerate() {
            let j = sorted[(idx + 1) % k];
            let steps = if k == 1 { m } else { (j + m - i) % m };
            let a: Vec<u32> = (0..=steps).map(|q| hole[(i + q) % m]).collect();
            let from = place[&a[0]];
            let to = place[a.last().unwrap()];
            let b = back(from, to, k == 1);
            jobs.push((a, b));
        }
        for (a, b) in jobs {
            if a.len() == 2 && b.len() == 2 {
                continue; // already twins
            }
            if a.len() == 2 {
                // The result's edge runs along the outline: split it at the outline's vertices if
                // they lie on it, else fill the sliver between.
                let (pa, pb) = (s.positions[a[0] as usize], s.positions[a[1] as usize]);
                let inner = &b[1..b.len() - 1];
                if inner.iter().all(|&v| closest_on_segment(s.positions[v as usize], pa, pb).1 <= 0.05 * h) {
                    s.split(a[0], a[1], inner)?;
                    continue;
                }
            }
            let filled = if ears { s.ears(&a, &b, h) } else { None };
            stats.bridged += match filled {
                Some(added) => added,
                None => s.ladder(&a, &b, 0.02 * h)?,
            };
            s.index();
        }
    }

    let mut out = Mesh { positions: s.positions, triangles: s.triangles };
    out.compact();
    let r = out.validate(check);
    if !r.ok() {
        return Err(format!(
            "the stitched mesh failed validation (open edges {}, non-manifold {}/{}, misoriented {}, degenerate {}, crossings {:?}, folds {})",
            r.boundary_edges, r.nonmanifold_edges, r.nonmanifold_vertices, r.misoriented_edges, r.degenerate_triangles, r.self_intersections, r.folds
        ));
    }
    Ok((out, stats))
}

/// The unit direction in which points spread least (the normal of their best-fit plane): the
/// eigenvector of the smallest eigenvalue of their covariance, by Jacobi rotations.
fn least_spread(points: &[[f64; 3]]) -> Option<[f64; 3]> {
    let n = points.len() as f64;
    if n < 3.0 {
        return None;
    }
    let mut c = [0.0; 3];
    for p in points {
        for d in 0..3 {
            c[d] += p[d] / n;
        }
    }
    let mut a = [[0.0f64; 3]; 3];
    for p in points {
        let q = sub(*p, c);
        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += q[i] * q[j];
            }
        }
    }
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let (mut p, mut q, mut big) = (0, 1, 0.0);
        for i in 0..3 {
            for j in i + 1..3 {
                if a[i][j].abs() > big {
                    (p, q, big) = (i, j, a[i][j].abs());
                }
            }
        }
        if big <= 1e-30 * (a[0][0].abs() + a[1][1].abs() + a[2][2].abs()).max(1e-300) {
            break;
        }
        let theta = 0.5 * (2.0 * a[p][q]).atan2(a[q][q] - a[p][p]);
        let (sn, cs) = theta.sin_cos();
        // a ← Jᵀ a J, v ← v J for the rotation J in the (p, q) plane.
        for k in 0..3 {
            let (akp, akq) = (a[k][p], a[k][q]);
            a[k][p] = cs * akp - sn * akq;
            a[k][q] = sn * akp + cs * akq;
        }
        for k in 0..3 {
            let (apk, aqk) = (a[p][k], a[q][k]);
            a[p][k] = cs * apk - sn * aqk;
            a[q][k] = sn * apk + cs * aqk;
        }
        for k in 0..3 {
            let (vkp, vkq) = (v[k][p], v[k][q]);
            v[k][p] = cs * vkp - sn * vkq;
            v[k][q] = sn * vkp + cs * vkq;
        }
    }
    let k = (0..3).min_by(|&i, &j| a[i][i].total_cmp(&a[j][j]))?;
    let e = [v[0][k], v[1][k], v[2][k]];
    let l = norm(e);
    (l > 0.0).then(|| e.map(|x| x / l))
}

fn v3(p: [f64; 3]) -> Vector {
    Vector::new(p[0], p[1], p[2])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn least_spread_finds_the_plane_normal() {
        let pts: Vec<[f64; 3]> = (0..40).map(|i| [i as f64 * 0.3, 15.0 + if i % 2 == 0 { 0.05 } else { 0.0 }, (i as f64 * 0.7).sin() * 5.0]).collect();
        let n = least_spread(&pts).unwrap();
        assert!(n[1].abs() > 0.999, "{n:?}");
    }
}
