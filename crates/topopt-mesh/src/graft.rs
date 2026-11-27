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

/// How the gaps between the result and an outline are filled (each falls back to a ladder).
#[derive(Clone, Copy, PartialEq)]
enum Fill {
    /// Ear clipping in the gap's plane when it is flat, else in its ribbon coordinates.
    Ears,
    /// Ear clipping in ribbon coordinates.
    Ribbon,
    /// Ladders only.
    Ladder,
}

/// Grafts with each way of filling the gaps in turn, then with wider cuts, until one passes
/// validation. The first failure is the one reported.
fn graft_either(mesh: &Mesh, original: &Original, keep: &[bool], h: f64, check: bool) -> Result<(Mesh, Stats), String> {
    let mut first = None;
    for band in [0.0, 0.5 * h, 1.0 * h] {
        for mode in [Fill::Ears, Fill::Ribbon, Fill::Ladder] {
            match graft(mesh, original, keep, h, check, mode, band) {
                Ok(done) => return Ok(done),
                Err(e) => {
                    first.get_or_insert(e);
                }
            }
        }
    }
    Err(first.unwrap_or_default())
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

/// Ear clipping of a simple polygon `poly` given by its points `q` in a plane, counter-clockwise:
/// the best-shaped ear first (area over the squared longest edge). `None` when it is not simple.
fn ear_clip(poly: &[u32], q: &[[f64; 2]]) -> Option<Vec<[u32; 3]>> {
    let cross2 = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    let scale = q.iter().map(|p| p[0].abs().max(p[1].abs())).fold(0.0, f64::max).max(1e-9);
    let eps = 1e-12 * scale * scale;
    let n = q.len();
    if n < 3 {
        return None;
    }
    // Only a simple polygon: no two edges that are not neighbours may cross or touch.
    let orient = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        let v = cross2(a, b, c);
        if v > eps { 1 } else if v < -eps { -1 } else { 0 }
    };
    for i in 0..n {
        for j in i + 2..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            let (a, b, c, d) = (q[i], q[(i + 1) % n], q[j], q[(j + 1) % n]);
            let (o1, o2, o3, o4) = (orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b));
            if o1 * o2 <= 0 && o3 * o4 <= 0 && !(o1 == 0 && o2 == 0) {
                return None;
            }
        }
    }
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
    Some(tris)
}

/// Removes islands: pieces of the result the cut has cut off from the rest. Small ones (64
/// triangles or fewer: a sliver of a wall next to the outline), and any lying all on or next to
/// the kept faces (a ring around a hole the cut went round). Their loops would match no outline.
/// Returns the triangles removed.
fn drop_islands(mesh: &Mesh, original: &Original, keep: &[bool], h: f64, rest: &mut Vec<[u32; 3]>) -> usize {
    let tmp = Mesh { positions: mesh.positions.clone(), triangles: rest.clone() };
    let (comp, count) = tmp.components();
    if count < 2 {
        return 0;
    }
    let mut size = vec![0usize; count];
    for &c in &comp {
        size[c as usize] += 1;
    }
    let largest = (0..count).max_by_key(|&c| size[c]).unwrap_or(0);
    let mut drop: Vec<bool> = (0..count).map(|c| c != largest).collect();
    for (t, tri) in rest.iter().enumerate() {
        let c = comp[t] as usize;
        if !drop[c] || size[c] <= 64 {
            continue;
        }
        for &v in tri {
            let p = mesh.positions[v as usize];
            let (proj, (k, _)) = original.shape.project_local_point_and_get_location(v3(p), false);
            if !keep[k as usize] || norm(sub(p, [proj.point.x, proj.point.y, proj.point.z])) > h {
                drop[c] = false;
                break;
            }
        }
    }
    let before = rest.len();
    let mut k = 0;
    rest.retain(|_| {
        let c = comp[k] as usize;
        k += 1;
        !drop[c]
    });
    before - rest.len()
}

/// A vertex with more than one boundary half-edge leaving it (two holes touching there), if any.
fn pinched_vertex(triangles: &[[u32; 3]]) -> Option<u32> {
    let mut half: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    for t in triangles {
        for e in 0..3 {
            half.insert((t[e], t[(e + 1) % 3]));
        }
    }
    let mut out: HashMap<u32, u32> = HashMap::new();
    let mut found = None;
    for &(a, b) in &half {
        if !half.contains(&(b, a)) {
            let n = out.entry(a).or_insert(0);
            *n += 1;
            if *n > 1 {
                found = Some(found.map_or(a, |f: u32| f.min(a)));
            }
        }
    }
    found
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
    /// Per triangle: one of the original surface's (or a piece of one), which mending must not
    /// change.
    patch: Vec<bool>,
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
                self.patch.push(self.patch[t]);
                self.triangles.len() - 1
            };
            for e in 0..3 {
                self.half.insert((nt[e], nt[(e + 1) % 3]), id);
            }
        }
        Ok(())
    }

    /// Mends folds (neighbours folded flat onto each other) and slivers along the seam by edge
    /// flips: an edge of a folded or degenerate triangle is flipped when both new triangles face
    /// the way the triangles around them do; an edge between two vertices on one point is
    /// collapsed. The original surface's triangles are never changed. Returns the repairs made.
    fn mend(&mut self, tiny: f64) -> usize {
        let mut flips = 0;
        // A triangle and its mirror image (the same corners, the other way round) are a flap of
        // no thickness: both go, and their edges are left with the triangles on either side. One
        // of the original surface's stays, the other goes.
        let mut by_corners: HashMap<[u32; 3], Vec<usize>> = HashMap::new();
        for (t, tri) in self.triangles.iter().enumerate() {
            let mut k = *tri;
            k.sort_unstable();
            by_corners.entry(k).or_default().push(t);
        }
        let mut gone = vec![false; self.triangles.len()];
        for list in by_corners.values().filter(|l| l.len() == 2) {
            let (a, b) = (list[0], list[1]);
            let (ta, tb) = (self.triangles[a], self.triangles[b]);
            let rot = |t: [u32; 3], k: usize| [t[k], t[(k + 1) % 3], t[(k + 2) % 3]];
            let mirrored = (0..3).any(|k| rot(tb, k) == [ta[0], ta[2], ta[1]]);
            if !mirrored {
                continue;
            }
            match (self.patch[a], self.patch[b]) {
                (false, false) => {
                    gone[a] = true;
                    gone[b] = true;
                }
                (true, false) => gone[b] = true,
                (false, true) => gone[a] = true,
                (true, true) => {}
            }
        }
        if gone.iter().any(|&g| g) {
            let mut k = 0;
            let patch = std::mem::take(&mut self.patch);
            self.triangles.retain(|_| {
                k += 1;
                !gone[k - 1]
            });
            self.patch = patch.into_iter().zip(&gone).filter(|(_, g)| !**g).map(|(p, _)| p).collect();
            flips += gone.iter().filter(|&&g| g).count();
        }
        for _ in 0..200 {
            let tmp = Mesh { positions: self.positions.clone(), triangles: self.triangles.clone() };
            let mut bad: Vec<u32> = crate::mesh::folded_triangles(&tmp);
            bad.extend(crate::mesh::degenerate_triangles(&tmp));
            if bad.is_empty() {
                break;
            }
            // Two vertices on one point (a vertex of the result not welded to the outline point
            // it sits on): collapse the edge between them, keeping the original surface's vertex.
            let on_patch: std::collections::HashSet<u32> =
                self.triangles.iter().zip(&self.patch).filter(|(_, p)| **p).flat_map(|(t, _)| *t).collect();
            let mut collapse = None;
            'pairs: for &t in &bad {
                let tri = self.triangles[t as usize];
                for e in 0..3 {
                    let (a, b) = (tri[e], tri[(e + 1) % 3]);
                    if norm(sub(self.positions[a as usize], self.positions[b as usize])) <= tiny {
                        if !on_patch.contains(&a) {
                            collapse = Some((a, b));
                            break 'pairs;
                        }
                        if !on_patch.contains(&b) {
                            collapse = Some((b, a));
                            break 'pairs;
                        }
                    }
                }
            }
            if let Some((gone, into)) = collapse {
                let mut kept = Vec::with_capacity(self.triangles.len());
                let mut flags = Vec::with_capacity(self.triangles.len());
                for (t, &p) in self.triangles.iter().zip(&self.patch) {
                    let t = t.map(|v| if v == gone { into } else { v });
                    if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
                        kept.push(t);
                        flags.push(p);
                    }
                }
                self.triangles = kept;
                self.patch = flags;
                flips += 1;
                continue;
            }
            self.index();
            let mut star: HashMap<u32, Vec<usize>> = HashMap::new();
            for (k, t) in self.triangles.iter().enumerate() {
                for &v in t {
                    star.entry(v).or_default().push(k);
                }
            }
            let normal = |s: &Self, t: [u32; 3]| {
                let p = |v: u32| s.positions[v as usize];
                cross(sub(p(t[1]), p(t[0])), sub(p(t[2]), p(t[0])))
            };
            let mut done = false;
            'search: for &t in &bad {
                let tri = self.triangles[t as usize];
                for e in 0..3 {
                    let (u, v, a) = (tri[e], tri[(e + 1) % 3], tri[(e + 2) % 3]);
                    let Some(&t2) = self.half.get(&(v, u)) else { continue };
                    if self.patch[t as usize] || self.patch[t2] {
                        continue;
                    }
                    let Some(b) = self.triangles[t2].iter().copied().find(|&c| c != u && c != v) else { continue };
                    if a == b || self.half.contains_key(&(a, b)) || self.half.contains_key(&(b, a)) {
                        continue;
                    }
                    // The way the surface faces here: the triangles around u and v but these two.
                    let mut r = [0.0; 3];
                    let mut around: Vec<usize> = star.get(&u).into_iter().chain(star.get(&v)).flatten().copied().collect();
                    around.sort_unstable();
                    around.dedup();
                    for k in around {
                        if k != t as usize && k != t2 {
                            let n = normal(self, self.triangles[k]);
                            r = [r[0] + n[0], r[1] + n[1], r[2] + n[2]];
                        }
                    }
                    let (n1, n2) = (normal(self, [u, b, a]), normal(self, [b, v, a]));
                    let l = (norm(n1), norm(n2), norm(r));
                    if l.0 <= 1e-12 * l.2.max(1e-300) || l.1 <= 1e-12 * l.2.max(1e-300) || dot(n1, r) <= 0.2 * l.0 * l.2 || dot(n2, r) <= 0.2 * l.1 * l.2 {
                        continue;
                    }
                    self.triangles[t as usize] = [u, b, a];
                    self.triangles[t2] = [b, v, a];
                    flips += 1;
                    done = true;
                    break 'search;
                }
            }
            if !done {
                // No flip helps: collapse an edge of a bad triangle instead, moving a vertex of the
                // result (never one of the original surface's) onto its neighbour, when the mesh
                // stays manifold (the link condition) and no triangle turns over.
                if let Some((gone, into)) = self.collapsible(&bad, &on_patch, &star) {
                    let mut kept = Vec::with_capacity(self.triangles.len());
                    let mut flags = Vec::with_capacity(self.triangles.len());
                    for (t, &p) in self.triangles.iter().zip(&self.patch) {
                        let t = t.map(|v| if v == gone { into } else { v });
                        if t[0] != t[1] && t[1] != t[2] && t[0] != t[2] {
                            kept.push(t);
                            flags.push(p);
                        }
                    }
                    self.triangles = kept;
                    self.patch = flags;
                    flips += 1;
                    continue;
                }
                break;
            }
        }
        self.index();
        flips
    }

    /// Untangles the seam: the result's own vertices (`< base`, on no patch triangle) within
    /// `band` of the outline move to the mean of their neighbours, `rounds` times; one that lay on
    /// the original surface goes back onto it, but never onto a kept face (it would sit on top of
    /// the patch: such a vertex stays put). Returns the vertices moved.
    fn relax(&mut self, base: u32, original: &Original, keep: &[bool], near_outline: &dyn Fn([f64; 3]) -> bool, h: f64, rounds: usize) -> usize {
        let on_patch: std::collections::HashSet<u32> =
            self.triangles.iter().zip(&self.patch).filter(|(_, p)| **p).flat_map(|(t, _)| *t).collect();
        let mut neighbours: HashMap<u32, Vec<u32>> = HashMap::new();
        for t in &self.triangles {
            for k in 0..3 {
                neighbours.entry(t[k]).or_default().push(t[(k + 1) % 3]);
            }
        }
        let movable: Vec<u32> = neighbours
            .keys()
            .copied()
            .filter(|&v| v < base && !on_patch.contains(&v) && near_outline(self.positions[v as usize]))
            .collect();
        let on_surface: HashMap<u32, bool> = movable
            .iter()
            .map(|&v| (v, original.shape.distance_to_local_point(v3(self.positions[v as usize]), false) <= 1e-6 * h))
            .collect();
        let mut moved = 0;
        for _ in 0..rounds {
            let next: Vec<(u32, [f64; 3])> = movable
                .iter()
                .map(|&v| {
                    let ns = &neighbours[&v];
                    let n = ns.len() as f64;
                    let mean = ns.iter().fold([0.0; 3], |a, &u| {
                        let p = self.positions[u as usize];
                        [a[0] + p[0] / n, a[1] + p[1] / n, a[2] + p[2] / n]
                    });
                    (v, mean)
                })
                .collect();
            for (v, mut p) in next {
                if on_surface[&v] {
                    let (proj, (t, _)) = original.shape.project_local_point_and_get_location(v3(p), false);
                    if keep[t as usize] {
                        continue;
                    }
                    p = [proj.point.x, proj.point.y, proj.point.z];
                }
                self.positions[v as usize] = p;
                moved += 1;
            }
        }
        moved
    }

    /// An edge of one of the `bad` triangles to collapse, (gone, into): `gone` is the result's own
    /// vertex, the two share exactly the two neighbours across the edge, and every triangle
    /// around `gone` that stays keeps facing the way it did.
    fn collapsible(&self, bad: &[u32], on_patch: &std::collections::HashSet<u32>, star: &HashMap<u32, Vec<usize>>) -> Option<(u32, u32)> {
        let normal = |t: [u32; 3]| {
            let p = |v: u32| self.positions[v as usize];
            cross(sub(p(t[1]), p(t[0])), sub(p(t[2]), p(t[0])))
        };
        let ring = |v: u32| -> std::collections::HashSet<u32> {
            star.get(&v).into_iter().flatten().flat_map(|&t| self.triangles[t]).filter(|&u| u != v).collect()
        };
        for &t in bad {
            let tri = self.triangles[t as usize];
            for k in 0..3 {
                for (gone, into) in [(tri[k], tri[(k + 1) % 3]), (tri[(k + 1) % 3], tri[k])] {
                    if on_patch.contains(&gone) {
                        continue;
                    }
                    let common: Vec<u32> = ring(gone).intersection(&ring(into)).copied().collect();
                    let across: Vec<u32> = star[&gone]
                        .iter()
                        .map(|&s| self.triangles[s])
                        .filter(|s| s.contains(&into))
                        .filter_map(|s| s.iter().copied().find(|&c| c != gone && c != into))
                        .collect();
                    if common.len() != across.len() || !common.iter().all(|c| across.contains(c)) {
                        continue;
                    }
                    let ok = star[&gone].iter().map(|&s| self.triangles[s]).filter(|s| !s.contains(&into)).all(|s| {
                        let after = s.map(|v| if v == gone { into } else { v });
                        let (n0, n1) = (normal(s), normal(after));
                        let l = norm(n1);
                        l > 1e-9 * norm(n0).max(1e-300) && dot(n0, n1) > 0.0
                    });
                    if ok {
                        return Some((gone, into));
                    }
                }
            }
        }
        None
    }

    /// Whether splitting the triangle on half-edge `u`→`v` at `points` (a fan from its third
    /// corner) gives triangles that all face the way it does, none of them degenerate.
    fn fan_fits(&self, u: u32, v: u32, points: &[u32]) -> bool {
        let Some(&t) = self.half.get(&(u, v)) else { return false };
        let tri = self.triangles[t];
        let Some(x) = tri.iter().copied().find(|&c| c != u && c != v) else { return false };
        let p = |i: u32| self.positions[i as usize];
        let n0 = cross(sub(p(v), p(u)), sub(p(x), p(u)));
        let l0 = norm(n0);
        if l0 <= 0.0 {
            return false;
        }
        let mut chain = vec![u];
        chain.extend_from_slice(points);
        chain.push(v);
        chain.windows(2).all(|w| {
            let n = cross(sub(p(w[1]), p(w[0])), sub(p(x), p(w[0])));
            dot(n, n0) > 1e-3 * norm(n) * l0 && norm(n) > 1e-9 * l0
        })
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
        let tris = ear_clip(&poly, &q)?;
        // In 3D each triangle must face the plane's way: a gap that bends (a strip around a
        // cylinder) can pass the flatness test and still fold once projected.
        let p3 = |v: u32| self.positions[v as usize];
        if tris.iter().any(|t| {
            let n = cross(sub(p3(t[1]), p3(t[0])), sub(p3(t[2]), p3(t[0])));
            dot(n, nn) <= 0.5 * norm(n)
        }) {
            return None;
        }
        let added = tris.len();
        self.triangles.extend(tris);
        self.patch.resize(self.triangles.len(), false);
        Some(added)
    }

    /// Fills the gap between chains `a` and `b` (as for [`Stitch::ladder`]) by ear clipping in
    /// the gap's own coordinates: arc length along `b` and distance from it. That unrolls a strip
    /// along the outline whatever the face next to it is (a plane, a hole's wall), so it works
    /// where the gap is not flat. `None` when the gap is not simple in those coordinates.
    fn ribbon(&mut self, a: &[u32], b: &[u32]) -> Option<usize> {
        if a.first() != b.first() || a.last() != b.last() || a.len() + b.len() < 5 || b.len() < 2 {
            return None;
        }
        let p3 = |s: &Self, v: u32| s.positions[v as usize];
        let mut sb = vec![0.0];
        for k in 1..b.len() {
            sb.push(sb[k - 1] + norm(sub(p3(self, b[k]), p3(self, b[k - 1]))));
        }
        let mut poly: Vec<u32> = b.to_vec();
        poly.extend(a[1..a.len() - 1].iter().rev());
        let mut q: Vec<[f64; 2]> = sb.iter().map(|&s| [s, 0.0]).collect();
        for &v in a[1..a.len() - 1].iter().rev() {
            let p = p3(self, v);
            let (mut at, mut best) = (0.0, f64::INFINITY);
            for k in 0..b.len() - 1 {
                let (t, d) = closest_on_segment(p, p3(self, b[k]), p3(self, b[k + 1]));
                if d < best {
                    best = d;
                    at = sb[k] + t * (sb[k + 1] - sb[k]);
                }
            }
            q.push([at, best]);
        }
        let area: f64 = (0..q.len()).map(|i| {
            let (p, r) = (q[i], q[(i + 1) % q.len()]);
            p[0] * r[1] - r[0] * p[1]
        }).sum();
        if area < 0.0 {
            q.iter_mut().for_each(|p| p[1] = -p[1]);
        }
        let tris = ear_clip(&poly, &q)?;
        if tris.iter().any(|t| norm(cross(sub(p3(self, t[1]), p3(self, t[0])), sub(p3(self, t[2]), p3(self, t[0])))) <= 0.0) {
            return None;
        }
        let added = tris.len();
        self.triangles.extend(tris);
        self.patch.resize(self.triangles.len(), false);
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
        self.patch.resize(self.triangles.len(), false);
        Ok(added)
    }
}

/// One graft of the kept triangles `keep` (all of one or more patches).
fn graft(mesh: &Mesh, original: &Original, keep: &[bool], h: f64, check: bool, mode: Fill, band: f64) -> Result<(Mesh, Stats), String> {
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
            && (d <= 0.05 * h || dot(n, original.normals[t]) > 0.3)
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
    // A vertex of the result on the outline (within the weld distance) must end up on the cut, to
    // be welded: one inside the remaining surface lies between the outline's chords and the true
    // edge, and its triangles would overlap what fills the gap. Cut its triangles too.
    let weld = 0.1 * h;
    for _ in 0..8 {
        let mut on_cut = vec![false; mesh.positions.len()];
        {
            let mut half: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::with_capacity(rest.len() * 3);
            for t in &rest {
                for e in 0..3 {
                    half.insert((t[e], t[(e + 1) % 3]));
                }
            }
            for &(a, b) in &half {
                if !half.contains(&(b, a)) {
                    on_cut[a as usize] = true;
                    on_cut[b as usize] = true;
                }
            }
        }
        let mut seen = vec![false; mesh.positions.len()];
        let mut inner: Vec<u32> = Vec::new();
        for t in &rest {
            for &v in t {
                let v = v as usize;
                if !seen[v] && !on_cut[v] {
                    seen[v] = true;
                    if outline_shape.distance_to_local_point(v3(mesh.positions[v]), false) <= weld {
                        inner.push(v as u32);
                    }
                }
            }
        }
        if inner.is_empty() {
            break;
        }
        let inner: std::collections::HashSet<u32> = inner.into_iter().collect();
        let before = rest.len();
        rest.retain(|t| !t.iter().any(|v| inner.contains(v)));
        stats.removed += before - rest.len();
    }
    // A wider cut, when asked for: every triangle with a corner within `band` of the outline, so
    // that the gap is filled afresh rather than stitched to a ragged edge.
    if band > 0.0 {
        let near: Vec<bool> = mesh.positions.iter().map(|&p| outline_shape.distance_to_local_point(v3(p), false) <= band).collect();
        let before = rest.len();
        rest.retain(|t| !t.iter().any(|&v| near[v as usize]));
        stats.removed += before - rest.len();
    }
    // Where two cuts touch at a vertex, cut its triangles too, until the cut's loops are simple.
    let holes = {
        let mut tries = 0;
        loop {
            stats.removed += drop_islands(mesh, original, keep, h, &mut rest);
            match pinched_vertex(&rest) {
                None => break boundary_loops(&rest)?,
                Some(v) if tries < 64 => {
                    let before = rest.len();
                    rest.retain(|t| !t.contains(&v));
                    stats.removed += before - rest.len();
                    tries += 1;
                }
                Some(_) => return Err("the cut touches itself at a vertex".into()),
            }
        }
    };

    // 2. Insert: the original's vertices go after the result's.
    let base = mesh.positions.len() as u32;
    let mut s = Stitch { positions: mesh.positions.clone(), patch: vec![false; rest.len()], triangles: rest, half: HashMap::new() };
    s.positions.extend_from_slice(opos);
    for t in &kept_tris {
        s.triangles.push(t.map(|v| v + base));
        s.patch.push(true);
    }
    stats.inserted = kept_tris.len();
    let mut outlines: Vec<Outline> = outline_loops.iter().map(|lp| Outline::new(lp.iter().map(|&v| v + base).collect(), &s.positions)).collect();

    // 3. Match each hole to the nearest outline loop, one to one: the nearest pairs first. A
    // small hole left over (a triangle or two cut apart from the rest) is closed where it is.
    let mut pairs: Vec<(f64, usize, usize)> = Vec::new();
    for (k, hole) in holes.iter().enumerate() {
        for (o, out) in outlines.iter().enumerate() {
            let mean = hole.iter().map(|&v| out.nearest(s.positions[v as usize], &s.positions).2).sum::<f64>() / hole.len() as f64;
            pairs.push((mean, k, o));
        }
    }
    // Small holes last: a stray one must not take a loop from the cut that runs along it.
    pairs.sort_by(|a, b| (holes[a.1].len() <= 12).cmp(&(holes[b.1].len() <= 12)).then(a.0.total_cmp(&b.0)));
    let mut owner = vec![usize::MAX; outlines.len()];
    let mut matched = vec![false; holes.len()];
    for &(mean, k, o) in &pairs {
        if owner[o] == usize::MAX && !matched[k] {
            if mean > 10.0 * h {
                return Err(format!("a cut is {:.2} samples from the outline", mean / h));
            }
            owner[o] = k;
            matched[k] = true;
        }
    }
    if owner.contains(&usize::MAX) {
        return Err(format!("{} holes cut for {} outline loops", holes.len(), outlines.len()));
    }
    for (k, hole) in holes.iter().enumerate() {
        if matched[k] {
            continue;
        }
        if hole.len() > 12 {
            return Err(format!("{} holes cut for {} outline loops", holes.len(), outlines.len()));
        }
        // Fill it: against the hole's direction, ear-clipped in its best-fit plane.
        let poly: Vec<u32> = hole.iter().rev().copied().collect();
        let pos: Vec<[f64; 3]> = poly.iter().map(|&v| s.positions[v as usize]).collect();
        let n = pos.len();
        let mut nn = [0.0; 3];
        for i in 0..n {
            let (p, q) = (pos[i], pos[(i + 1) % n]);
            nn[0] += (p[1] - q[1]) * (p[2] + q[2]);
            nn[1] += (p[2] - q[2]) * (p[0] + q[0]);
            nn[2] += (p[0] - q[0]) * (p[1] + q[1]);
        }
        let mut axis = least_spread(&pos).ok_or("a small cut does not close")?;
        if dot(axis, nn) < 0.0 {
            axis = axis.map(|c| -c);
        }
        let t = if axis[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let u = cross(axis, t);
        let u = u.map(|c| c / norm(u));
        let w = cross(axis, u);
        let q: Vec<[f64; 2]> = pos.iter().map(|&p| [dot(p, u), dot(p, w)]).collect();
        let tris = ear_clip(&poly, &q).ok_or("a small cut does not close")?;
        stats.bridged += tris.len();
        s.triangles.extend(tris);
        s.patch.resize(s.triangles.len(), false);
    }

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
            if a.len() == 2 && s.fan_fits(a[0], a[1], &b[1..b.len() - 1]) {
                // The result's edge runs along the outline, whose vertices between lie on it or
                // inside its triangle: split the triangle at them.
                s.split(a[0], a[1], &b[1..b.len() - 1])?;
                continue;
            }
            let filled = match mode {
                Fill::Ears => s.ears(&a, &b, h).or_else(|| s.ribbon(&a, &b)),
                Fill::Ribbon => s.ribbon(&a, &b),
                Fill::Ladder => None,
            };
            stats.bridged += match filled {
                Some(added) => added,
                None => s.ladder(&a, &b, 0.02 * h)?,
            };
            s.index();
        }
    }

    s.mend(1e-6 * h);
    // Still folded or crossing next to the seam: untangle there and mend again.
    {
        let tmp = Mesh { positions: s.positions.clone(), triangles: s.triangles.clone() };
        let tangled = !crate::mesh::folded_triangles(&tmp).is_empty()
            || !crate::mesh::degenerate_triangles(&tmp).is_empty()
            || (check && !crate::intersect::find_self_intersections(&tmp).is_empty());
        if tangled {
            let near = |p: [f64; 3]| outline_shape.distance_to_local_point(v3(p), false) <= 2.5 * h;
            s.relax(base, original, keep, &near, h, 12);
            s.mend(1e-6 * h);
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
