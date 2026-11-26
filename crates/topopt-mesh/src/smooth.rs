//! Smoothing of the contour and the projection onto the original part surface.
//!
//! Smoothers (all take a per-vertex `constrain` hook, applied after every pass, that may replace a
//! vertex's position, e.g. to keep it on the original surface):
//! - [`taubin`]: λ|μ umbrella smoothing (Taubin 1995). Removes the voxel staircase; shrinks a
//!   little.
//! - [`hc_laplacian`]: the HC algorithm (Vollmer, Mencl & Müller 1999): Laplacian smoothing whose
//!   every step is pushed back towards the previous and the original positions, so it neither
//!   shrinks nor drifts.
//! - [`bilateral`]: bilateral mesh filtering on the triangles' normals (Zheng et al. 2011), then a
//!   vertex update that fits the filtered normals (Sun et al. 2007). The weights fall off with
//!   distance and with the difference of normals, so creases are kept.
//!
//! [`restore_volume`] moves the free vertices along their normals until the enclosed volume is back
//! to a target: the "rescaling" that makes any smoother volume-preserving.
//!
//! [`Projector`] finds the closest point on the original surface, and also on its sharp feature
//! edges and corners (dihedral angle above a threshold), so vertices along a crease of the part can
//! be pinned onto the crease itself instead of onto one of its faces.

use crate::mesh::{Mesh, sub, vertex_normals};
use parry3d_f64::math::Vector;
use parry3d_f64::query::{PointQuery, PointQueryWithLocation};
use parry3d_f64::shape::{Polyline, TriMesh};
use std::collections::HashMap;

/// A per-vertex constraint: `Some(p)` replaces vertex `v`'s new position.
pub type Constrain<'a> = &'a dyn Fn(usize, [f64; 3]) -> Option<[f64; 3]>;

/// Which smoother [`crate::extract`] uses.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Smoother {
    /// Taubin λ|μ.
    #[default]
    Taubin,
    /// HC Laplacian (Vollmer et al.).
    Hc,
    /// A few Taubin passes for the staircase, then bilateral normal filtering, which smooths the
    /// remaining ripples while keeping creases.
    Bilateral,
}

impl Smoother {
    pub fn parse(name: &str) -> Option<Smoother> {
        match name {
            "taubin" => Some(Smoother::Taubin),
            "hc" => Some(Smoother::Hc),
            "bilateral" => Some(Smoother::Bilateral),
            _ => None,
        }
    }
}

/// Unique neighbours of every vertex.
pub fn neighbours(mesh: &Mesh) -> Vec<Vec<u32>> {
    let mut nb: Vec<Vec<u32>> = vec![Vec::new(); mesh.positions.len()];
    for t in &mesh.triangles {
        for e in 0..3 {
            let (a, b) = (t[e], t[(e + 1) % 3]);
            nb[a as usize].push(b);
            nb[b as usize].push(a);
        }
    }
    for list in &mut nb {
        list.sort_unstable();
        list.dedup();
    }
    nb
}

fn mean_of(positions: &[[f64; 3]], list: &[u32]) -> [f64; 3] {
    let mut c = [0.0; 3];
    for &u in list {
        let q = positions[u as usize];
        for d in 0..3 {
            c[d] += q[d];
        }
    }
    let inv = 1.0 / list.len() as f64;
    [c[0] * inv, c[1] * inv, c[2] * inv]
}

/// `iterations` λ/μ pairs. Taubin's alternating shrink/inflate removes the voxel staircase while
/// keeping the volume nearly constant (λ = 0.5, μ = −0.53).
pub fn taubin(mesh: &mut Mesh, iterations: usize, lambda: f64, mu: f64, constrain: Constrain) {
    let nb = neighbours(mesh);
    let mut next = mesh.positions.clone();
    for _ in 0..iterations {
        for &factor in &[lambda, mu] {
            for v in 0..mesh.positions.len() {
                if nb[v].is_empty() {
                    continue;
                }
                let p = mesh.positions[v];
                let c = mean_of(&mesh.positions, &nb[v]);
                let moved = [p[0] + factor * (c[0] - p[0]), p[1] + factor * (c[1] - p[1]), p[2] + factor * (c[2] - p[2])];
                next[v] = constrain(v, moved).unwrap_or(moved);
            }
            std::mem::swap(&mut mesh.positions, &mut next);
        }
    }
}

/// The HC algorithm: one Laplacian step `p`, then the difference `b = p − (α·o + (1−α)·q)` to the
/// original `o` and previous `q` positions is pushed back, mixed with the neighbours' differences:
/// `p − (β·b + (1−β)·mean(b of neighbours))`. α = 0 lets the shape move freely but without
/// shrinking; α → 1 ties it to the original. Vollmer et al. suggest α ≈ 0…0.2, β ≈ 0.5…0.9.
pub fn hc_laplacian(mesh: &mut Mesh, iterations: usize, alpha: f64, beta: f64, constrain: Constrain) {
    let nb = neighbours(mesh);
    let original = mesh.positions.clone();
    let n = mesh.positions.len();
    let mut p = vec![[0.0; 3]; n];
    let mut b = vec![[0.0; 3]; n];
    for _ in 0..iterations {
        let q = mesh.positions.clone();
        for v in 0..n {
            p[v] = if nb[v].is_empty() { q[v] } else { mean_of(&q, &nb[v]) };
            for d in 0..3 {
                b[v][d] = p[v][d] - (alpha * original[v][d] + (1.0 - alpha) * q[v][d]);
            }
        }
        for v in 0..n {
            let mb = if nb[v].is_empty() { [0.0; 3] } else { mean_of(&b, &nb[v]) };
            let moved: [f64; 3] = std::array::from_fn(|d| p[v][d] - (beta * b[v][d] + (1.0 - beta) * mb[d]));
            mesh.positions[v] = constrain(v, moved).unwrap_or(moved);
        }
    }
}

/// Tangential relaxation: each vertex moves towards the mean of its neighbours by `factor`, but only
/// within its tangent plane, so the shape stays and the vertices spread evenly (the bilateral
/// vertex update moves them only along normals).
fn relax_tangentially(mesh: &mut Mesh, nb: &[Vec<u32>], factor: f64, constrain: Constrain) {
    let normals = vertex_normals(mesh);
    let pos = mesh.positions.clone();
    for v in 0..pos.len() {
        if nb[v].is_empty() {
            continue;
        }
        let c = mean_of(&pos, &nb[v]);
        let d = sub(c, pos[v]);
        let n = normals[v];
        let along = dot(d, n);
        let t = [d[0] - along * n[0], d[1] - along * n[1], d[2] - along * n[2]];
        let moved = [pos[v][0] + factor * t[0], pos[v][1] + factor * t[1], pos[v][2] + factor * t[2]];
        mesh.positions[v] = constrain(v, moved).unwrap_or(moved);
    }
}

/// Bilateral mesh filtering, in its robust normal-filtering form (Zheng, Fu, Au & Tai 2011): each
/// triangle's normal is replaced by a weighted mean of its neighbours' normals, the weights
/// falling off with the distance between centroids (`sigma_c`) and with the difference between
/// normals (`sigma_r`, e.g. 0.35), so normals on the two sides of a crease do not mix. Then the
/// vertices are moved to fit the filtered normals (Sun, Rosin, Martin & Langbein 2007: each vertex
/// moves to the mean of its projections onto its triangles' planes). Each vertex only moves
/// towards the planes of its own triangles, so the filter is stable on irregular meshes.
pub fn bilateral(mesh: &mut Mesh, iterations: usize, sigma_c: f64, sigma_r: f64, constrain: Constrain) {
    let nv = mesh.positions.len();
    let nf = mesh.triangles.len();
    let mut vfaces: Vec<Vec<u32>> = vec![Vec::new(); nv];
    for (f, t) in mesh.triangles.iter().enumerate() {
        for &v in t {
            vfaces[v as usize].push(f as u32);
        }
    }
    // Neighbouring triangles: those sharing a vertex.
    let fnb: Vec<Vec<u32>> = (0..nf)
        .map(|f| {
            let mut r: Vec<u32> = mesh.triangles[f].iter().flat_map(|&v| vfaces[v as usize].iter().copied()).collect();
            r.sort_unstable();
            r.dedup();
            r
        })
        .collect();
    let (two_c2, two_r2) = (2.0 * sigma_c * sigma_c, 2.0 * sigma_r * sigma_r);
    let nb = neighbours(mesh);
    for _ in 0..iterations {
        let pos = &mesh.positions;
        let mut normal = vec![[0.0; 3]; nf];
        let mut centroid = vec![[0.0; 3]; nf];
        let mut area = vec![0.0; nf];
        for (f, t) in mesh.triangles.iter().enumerate() {
            let [a, b, c] = t.map(|i| pos[i as usize]);
            let n = cross(sub(b, a), sub(c, a));
            let l = dot(n, n).sqrt();
            area[f] = 0.5 * l;
            normal[f] = if l > 0.0 { [n[0] / l, n[1] / l, n[2] / l] } else { [0.0; 3] };
            centroid[f] = std::array::from_fn(|d| (a[d] + b[d] + c[d]) / 3.0);
        }
        let filtered: Vec<[f64; 3]> = (0..nf)
            .map(|i| {
                let mut s = [0.0; 3];
                for &j in &fnb[i] {
                    let j = j as usize;
                    let dc = sub(centroid[i], centroid[j]);
                    let dn = sub(normal[i], normal[j]);
                    let w = area[j] * (-dot(dc, dc) / two_c2).exp() * (-dot(dn, dn) / two_r2).exp();
                    for d in 0..3 {
                        s[d] += w * normal[j][d];
                    }
                }
                let u = unit(s);
                if u == [0.0; 3] { normal[i] } else { u }
            })
            .collect();
        // Fit the vertices to the filtered normals (a few small steps, centroids refreshed).
        for _ in 0..3 {
            let pos = mesh.positions.clone();
            let cent: Vec<[f64; 3]> = mesh
                .triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|i| pos[i as usize]);
                    std::array::from_fn(|d| (a[d] + b[d] + c[d]) / 3.0)
                })
                .collect();
            for v in 0..nv {
                if vfaces[v].is_empty() {
                    continue;
                }
                let p = pos[v];
                let mut step = [0.0; 3];
                for &f in &vfaces[v] {
                    let (n, c) = (filtered[f as usize], cent[f as usize]);
                    let h = dot(n, sub(c, p));
                    for d in 0..3 {
                        step[d] += h * n[d];
                    }
                }
                let k = 1.0 / vfaces[v].len() as f64;
                let moved = [p[0] + k * step[0], p[1] + k * step[1], p[2] + k * step[2]];
                mesh.positions[v] = constrain(v, moved).unwrap_or(moved);
            }
        }
        relax_tangentially(mesh, &nb, 0.5, constrain);
    }
}

/// Moves the vertices with `free[v]` along their normals by a common offset until the mesh
/// encloses `target` volume again (Newton steps on the offset, with the exact derivative
/// dV/dδ = Σ ∂V/∂p·n). Returns the final relative volume error. Vertices that are not free (pinned
/// to the original surface) keep their positions, so they stay exact. The total offset is at most
/// `max_offset`: when few vertices are free, they cannot carry the whole correction without
/// growing bumps.
pub fn restore_volume(mesh: &mut Mesh, target: f64, free: &[bool], max_offset: f64) -> f64 {
    let mut total = 0.0f64;
    for _ in 0..8 {
        let volume = mesh.volume();
        let err = target - volume;
        if err.abs() <= 1e-9 * target.abs() {
            break;
        }
        let normals = vertex_normals(mesh);
        // ∂V/∂p for each vertex: V = Σ a·(b×c)/6 over the triangles (a, b, c).
        let mut grad = vec![[0.0; 3]; mesh.positions.len()];
        for t in &mesh.triangles {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            for (v, g) in [(t[0], cross(b, c)), (t[1], cross(c, a)), (t[2], cross(a, b))] {
                for d in 0..3 {
                    grad[v as usize][d] += g[d] / 6.0;
                }
            }
        }
        let slope: f64 = (0..grad.len()).filter(|&v| free[v]).map(|v| dot(grad[v], normals[v])).sum();
        if slope <= 0.0 {
            break;
        }
        let offset = (total + err / slope).clamp(-max_offset, max_offset) - total;
        if offset == 0.0 {
            break;
        }
        total += offset;
        for v in 0..mesh.positions.len() {
            if free[v] {
                for d in 0..3 {
                    mesh.positions[v][d] += offset * normals[v][d];
                }
            }
        }
    }
    (mesh.volume() - target) / target.abs().max(1e-300)
}

/// Inserts vertices on the sharp edges of `surface`. Where the contour rounds or chamfers a crease
/// of the original part, its mesh edges run from one face of the crease to the other without a
/// vertex on the crease itself. Every mesh edge that crosses the bisector plane of a sharp edge
/// within `reach` of it is split at the crease (the crossing point, projected onto the sharp
/// edge), so the surface can bend exactly there. Splitting keeps the mesh closed and manifold.
/// Returns the new vertices (they lie on the sharp edges).
pub fn insert_crease_vertices(mesh: &mut Mesh, surface: &Projector, reach: f64) -> Vec<u32> {
    let mut edge_tris: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
    let key = |a: u32, b: u32| (a.min(b), a.max(b));
    for (t, tri) in mesh.triangles.iter().enumerate() {
        for e in 0..3 {
            edge_tris.entry(key(tri[e], tri[(e + 1) % 3])).or_default().push(t);
        }
    }
    // Decide the splits on the original edges first.
    let mut splits: Vec<((u32, u32), [f64; 3])> = Vec::new();
    let mut keys: Vec<_> = edge_tris.keys().copied().collect();
    keys.sort_unstable();
    for (u, w) in keys {
        let (pu, pw) = (mesh.positions[u as usize], mesh.positions[w as usize]);
        let mid = [(pu[0] + pw[0]) / 2.0, (pu[1] + pw[1]) / 2.0, (pu[2] + pw[2]) / 2.0];
        let Some(e) = surface.nearest_edge(mid) else { continue };
        if e.distance > reach {
            continue;
        }
        // Only the contour's own rounding of the crease: both ends lie on the part's surface (a
        // strut that merely passes near the crease is left alone). Work with the ends as they
        // will be once pinned onto the surface: the contour may run a little inside it.
        let ((qu, du), (qw, dw)) = (surface.project(pu), surface.project(pw));
        if du > 0.75 * reach || dw > 0.75 * reach {
            continue;
        }
        let (pu, pw) = (qu, qw);
        let along = unit(sub(pw, pu));
        if dot(along, e.direction).abs() > 0.9 {
            continue; // runs along the crease, not across it
        }
        let [n0, n1] = e.normals;
        let bis = unit(sub(n0, n1));
        if bis == [0.0; 3] {
            continue; // a boundary edge: no two sides
        }
        let (su, sw) = (dot(bis, sub(pu, e.point)), dot(bis, sub(pw, e.point)));
        if su * sw >= 0.0 {
            continue;
        }
        let t = su / (su - sw);
        let x = std::array::from_fn(|d| pu[d] + t * (pw[d] - pu[d]));
        let Some(on) = surface.nearest_edge(x) else { continue };
        // The new vertex may not jump far off the edge it splits (that would fold its triangles).
        let len = dist(pu, pw);
        if on.distance > reach.min(len) || !(0.05..=0.95).contains(&t) {
            continue;
        }
        splits.push(((u, w), on.point));
    }
    let mut added = Vec::with_capacity(splits.len());
    for ((u, w), point) in splits {
        let Some(tris) = edge_tris.get(&key(u, w)).cloned() else { continue };
        if tris.len() != 2 {
            continue;
        }
        let c = mesh.positions.len() as u32;
        mesh.positions.push(point);
        added.push(c);
        edge_tris.remove(&key(u, w));
        for t in tris {
            let tri = mesh.triangles[t];
            // Rotate so the split edge is (tri[0], tri[1]) in this triangle's winding.
            let r = (0..3).find(|&e| key(tri[e], tri[(e + 1) % 3]) == key(u, w)).unwrap();
            let [a, b, apex] = [tri[r], tri[(r + 1) % 3], tri[(r + 2) % 3]];
            let t2 = mesh.triangles.len();
            mesh.triangles[t] = [a, c, apex];
            mesh.triangles.push([c, b, apex]);
            // Update the edge → triangles map: (b, apex) moves from t to t2.
            if let Some(list) = edge_tris.get_mut(&key(b, apex)) {
                for x in list.iter_mut() {
                    if *x == t {
                        *x = t2;
                    }
                }
            }
            edge_tris.entry(key(a, c)).or_default().push(t);
            edge_tris.entry(key(c, b)).or_default().push(t2);
            edge_tris.entry(key(c, apex)).or_default().extend([t, t2]);
        }
    }
    added
}

/// Flips mesh edges so that the creases of the original surface become mesh edges. A contour
/// triangle pair often straddles a crease diagonally: its two apexes lie on the crease while the
/// shared edge cuts across it, so the surface sags between the crease vertices. An edge (a, b)
/// whose opposite vertices c and d are both `on_crease`, and whose midpoint of c–d lies on a
/// sharp edge of `surface` (within `tol`), is replaced by c–d when the flip keeps the mesh
/// manifold and does not turn a triangle over. Returns the number of flips.
pub fn flip_to_creases(mesh: &mut Mesh, on_crease: &[bool], surface: &Projector, tol: f64) -> usize {
    let mut flips = 0;
    for _ in 0..8 {
        let mut edge_tris: HashMap<(u32, u32), Vec<usize>> = HashMap::new();
        let mut neighbours: HashMap<u32, Vec<u32>> = HashMap::new();
        for (t, tri) in mesh.triangles.iter().enumerate() {
            for e in 0..3 {
                let (a, b) = (tri[e], tri[(e + 1) % 3]);
                edge_tris.entry((a.min(b), a.max(b))).or_default().push(t);
                neighbours.entry(a).or_default().push(b);
            }
        }
        let has_edge = |a: u32, b: u32, nb: &HashMap<u32, Vec<u32>>| nb.get(&a).is_some_and(|l| l.contains(&b)) || nb.get(&b).is_some_and(|l| l.contains(&a));
        let mut touched = vec![false; mesh.triangles.len()];
        let mut round = 0;
        let mut keys: Vec<_> = edge_tris.keys().copied().collect();
        keys.sort_unstable();
        for (a, b) in keys {
            let ts = &edge_tris[&(a, b)];
            let [t1, t2] = ts[..] else { continue };
            if touched[t1] || touched[t2] {
                continue;
            }
            // Orient: t1 = (a, b, c) and t2 = (b, a, d) in their own winding.
            let apex = |t: usize, from: u32, to: u32| -> Option<u32> {
                let tri = mesh.triangles[t];
                (0..3).find(|&e| tri[e] == from && tri[(e + 1) % 3] == to).map(|e| tri[(e + 2) % 3])
            };
            let (a, b) = if apex(t1, a, b).is_some() { (a, b) } else { (b, a) };
            let (Some(c), Some(d)) = (apex(t1, a, b), apex(t2, b, a)) else { continue };
            if !on_crease[c as usize] || !on_crease[d as usize] || (on_crease[a as usize] && on_crease[b as usize]) || c == d {
                continue;
            }
            if has_edge(c, d, &neighbours) {
                continue;
            }
            let [pa, pb, pc, pd] = [a, b, c, d].map(|v| mesh.positions[v as usize]);
            let mid = [(pc[0] + pd[0]) / 2.0, (pc[1] + pd[1]) / 2.0, (pc[2] + pd[2]) / 2.0];
            if surface.nearest_edge(mid).is_none_or(|e| e.distance > tol) {
                continue;
            }
            // New triangles (a, d, c) and (d, b, c); each must face the way the old pair did.
            let old = [cross(sub(pb, pa), sub(pc, pa)), cross(sub(pa, pb), sub(pd, pb))];
            let n1 = cross(sub(pd, pa), sub(pc, pa));
            let n2 = cross(sub(pb, pd), sub(pc, pd));
            let len = |v: [f64; 3]| dot(v, v).sqrt();
            if len(n1) < 1e-12 * dot(sub(pc, pd), sub(pc, pd)) || len(n2) < 1e-12 * dot(sub(pc, pd), sub(pc, pd)) {
                continue;
            }
            let agree = |n: [f64; 3]| old.iter().any(|o| dot(unit(n), unit(*o)) > 0.2);
            if !agree(n1) || !agree(n2) {
                continue;
            }
            mesh.triangles[t1] = [a, d, c];
            mesh.triangles[t2] = [d, b, c];
            touched[t1] = true;
            touched[t2] = true;
            neighbours.entry(c).or_default().push(d);
            round += 1;
        }
        flips += round;
        if round == 0 {
            break;
        }
    }
    flips
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    if l > 0.0 { [a[0] / l, a[1] / l, a[2] / l] } else { [0.0; 3] }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(a, b);
    dot(d, d).sqrt()
}

fn v3(p: [f64; 3]) -> Vector {
    Vector::new(p[0], p[1], p[2])
}

/// The nearest point on a feature, with the distance to it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FeatureHit {
    pub point: [f64; 3],
    pub distance: f64,
    /// Unit direction of the feature edge (zero for a corner).
    pub direction: [f64; 3],
    /// Unit normals of the two faces that meet at the edge (equal for a boundary edge).
    pub normals: [[f64; 3]; 2],
}

/// Closest-point queries on the original surface (a triangle mesh, from `model.mesh` or a spec's
/// domain), on its sharp edges and on its corners.
pub struct Projector {
    mesh: TriMesh,
    face_normals: Vec<[f64; 3]>,
    edges: Option<Polyline>,
    /// Per feature segment: its direction and its two faces' normals.
    edge_info: Vec<([f64; 3], [[f64; 3]; 2])>,
    corners: Vec<[f64; 3]>,
}

impl Projector {
    /// A projector whose sharp edges are those with a dihedral angle above 30°.
    pub fn new(positions: &[[f64; 3]], triangles: &[[u32; 3]]) -> Result<Self, String> {
        Self::with_feature_angle(positions, triangles, 30.0)
    }

    /// `feature_angle` in degrees: edges whose two faces' normals differ by more are sharp
    /// (0 or less: no features). Corners are the vertices where one, or three or more, sharp edges
    /// end, or where two meet at an angle sharper than `feature_angle`.
    pub fn with_feature_angle(positions: &[[f64; 3]], triangles: &[[u32; 3]], feature_angle: f64) -> Result<Self, String> {
        let verts = positions.iter().map(|&p| v3(p)).collect();
        let mesh = TriMesh::new(verts, triangles.to_vec()).map_err(|e| format!("surface mesh: {e:?}"))?;
        // Weld by position: CAD meshes may duplicate vertices along face boundaries.
        let mut weld: HashMap<[u64; 3], u32> = HashMap::new();
        let ids: Vec<u32> = positions
            .iter()
            .map(|p| {
                let n = weld.len() as u32;
                *weld.entry(p.map(|c| (c + 0.0).to_bits())).or_insert(n)
            })
            .collect();
        let mut rep = vec![[0.0; 3]; weld.len()];
        for (i, &id) in ids.iter().enumerate() {
            rep[id as usize] = positions[i];
        }
        let mut face_normals: Vec<[f64; 3]> =
            triangles.iter().map(|t| unit(cross(sub(positions[t[1] as usize], positions[t[0] as usize]), sub(positions[t[2] as usize], positions[t[0] as usize])))).collect();
        // Outward normals whatever the file's winding: a negative enclosed volume means the
        // triangles run inside out.
        let volume: f64 = triangles
            .iter()
            .map(|t| dot(positions[t[0] as usize], cross(positions[t[1] as usize], positions[t[2] as usize])) / 6.0)
            .sum();
        if volume < 0.0 {
            face_normals.iter_mut().for_each(|n| *n = n.map(|c| -c));
        }
        let mut edge_faces: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
        for (f, t) in triangles.iter().enumerate() {
            for e in 0..3 {
                let (a, b) = (ids[t[e] as usize], ids[t[(e + 1) % 3] as usize]);
                if a != b {
                    edge_faces.entry((a.min(b), a.max(b))).or_default().push(f as u32);
                }
            }
        }
        let cos_limit = feature_angle.to_radians().cos();
        let mut seg_vertices: Vec<[u32; 2]> = Vec::new();
        let mut edge_info = Vec::new();
        if feature_angle > 0.0 {
            let mut list: Vec<_> = edge_faces.iter().collect();
            list.sort_by_key(|(k, _)| **k);
            for (&(a, b), faces) in list {
                let normals = match faces.as_slice() {
                    [f] => [face_normals[*f as usize]; 2],
                    [f, g] => [face_normals[*f as usize], face_normals[*g as usize]],
                    _ => continue, // non-manifold: not a crease we can follow
                };
                if faces.len() == 2 && dot(normals[0], normals[1]) >= cos_limit {
                    continue;
                }
                seg_vertices.push([a, b]);
                edge_info.push((unit(sub(rep[b as usize], rep[a as usize])), normals));
            }
        }
        // Corners: where the sharp edges end, branch or turn sharply.
        let mut at: HashMap<u32, Vec<usize>> = HashMap::new();
        for (s, &[a, b]) in seg_vertices.iter().enumerate() {
            at.entry(a).or_default().push(s);
            at.entry(b).or_default().push(s);
        }
        let mut corners = Vec::new();
        let mut keys: Vec<_> = at.keys().copied().collect();
        keys.sort_unstable();
        for v in keys {
            let segs = &at[&v];
            let corner = match segs.as_slice() {
                [s, t] => {
                    // Directions pointing away from v.
                    let away = |s: usize| {
                        let [a, b] = seg_vertices[s];
                        let other = if a == v { b } else { a };
                        unit(sub(rep[other as usize], rep[v as usize]))
                    };
                    // A straight crease has the two directions opposite (dot −1).
                    dot(away(*s), away(*t)) > -cos_limit
                }
                _ => true,
            };
            if corner {
                corners.push(rep[v as usize]);
            }
        }
        let edges = if seg_vertices.is_empty() { None } else { Some(Polyline::new(rep.iter().map(|&p| v3(p)).collect(), Some(seg_vertices))) };
        Ok(Projector { mesh, face_normals, edges, edge_info, corners })
    }

    /// The closest point on the surface and the distance to it.
    pub fn project(&self, p: [f64; 3]) -> ([f64; 3], f64) {
        let q = self.mesh.project_local_point(v3(p), false).point;
        let c = [q.x, q.y, q.z];
        (c, dist(c, p))
    }

    /// The closest point on the surface, the distance, and the outward normal of the triangle it
    /// lies on.
    pub fn project_with_normal(&self, p: [f64; 3]) -> ([f64; 3], f64, [f64; 3]) {
        let (proj, (tri, _)) = self.mesh.project_local_point_and_get_location(v3(p), false);
        let c = [proj.point.x, proj.point.y, proj.point.z];
        (c, dist(c, p), self.face_normals.get(tri as usize).copied().unwrap_or([0.0; 3]))
    }

    /// The closest point on a sharp edge, if the surface has any.
    pub fn nearest_edge(&self, p: [f64; 3]) -> Option<FeatureHit> {
        let edges = self.edges.as_ref()?;
        let (proj, (seg, _)) = edges.project_local_point_and_get_location(v3(p), false);
        let c = [proj.point.x, proj.point.y, proj.point.z];
        let (direction, normals) = self.edge_info[seg as usize];
        Some(FeatureHit { point: c, distance: dist(c, p), direction, normals })
    }

    /// The closest corner, if the surface has any.
    pub fn nearest_corner(&self, p: [f64; 3]) -> Option<FeatureHit> {
        self.corners
            .iter()
            .map(|&c| (c, dist(c, p)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(c, d)| FeatureHit { point: c, distance: d, direction: [0.0; 3], normals: [[0.0; 3]; 2] })
    }

    pub fn corners(&self) -> &[[f64; 3]] {
        &self.corners
    }

    /// Number of sharp edge segments.
    pub fn feature_segments(&self) -> usize {
        self.edge_info.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube() -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
        let c = |i: usize| [(i & 1) as f64 * 10.0, ((i >> 1) & 1) as f64 * 10.0, ((i >> 2) & 1) as f64 * 10.0];
        let p = vec![c(0), c(1), c(3), c(2), c(4), c(5), c(7), c(6)];
        let t = vec![[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]];
        (p, t)
    }

    #[test]
    fn cube_has_twelve_edges_and_eight_corners() {
        let (p, t) = cube();
        let pr = Projector::new(&p, &t).unwrap();
        // The face diagonals are flat (0°): only the 12 cube edges are sharp.
        assert_eq!(pr.feature_segments(), 12);
        assert_eq!(pr.corners().len(), 8);
        let e = pr.nearest_edge([5.0, -1.0, -1.0]).unwrap();
        assert!(dist(e.point, [5.0, 0.0, 0.0]) < 1e-12);
        assert!(dot(e.direction, [1.0, 0.0, 0.0]).abs() > 1.0 - 1e-12);
        let c = pr.nearest_corner([11.0, 11.0, 9.0]).unwrap();
        assert!(dist(c.point, [10.0, 10.0, 10.0]) < 1e-12);
        let (_, d, n) = pr.project_with_normal([5.0, 5.0, 12.0]);
        assert!((d - 2.0).abs() < 1e-12 && (n[2] - 1.0).abs() < 1e-12);
        assert_eq!(Projector::with_feature_angle(&p, &t, 0.0).unwrap().feature_segments(), 0);
    }

    #[test]
    fn restore_volume_hits_the_target() {
        let (p, t) = cube();
        let mut m = Mesh { positions: p, triangles: t };
        let target = m.volume() * 1.05;
        let err = restore_volume(&mut m, target, &vec![true; 8], f64::INFINITY);
        assert!(err.abs() < 1e-6, "{err}");
    }
}
