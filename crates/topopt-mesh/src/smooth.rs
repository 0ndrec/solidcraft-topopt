//! Taubin λ|μ smoothing (umbrella operator) with per-vertex constraints, and the projection onto
//! the original part surface used as that constraint.

use crate::mesh::Mesh;
use parry3d_f64::math::Vector;
use parry3d_f64::query::PointQuery;
use parry3d_f64::shape::TriMesh;

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

/// `iterations` λ/μ pairs. After every half-step `constrain(v, p)` may replace a vertex position
/// (`Some`), e.g. to keep it on the original surface. Taubin's alternating shrink/inflate removes
/// the voxel staircase while keeping the volume nearly constant (λ = 0.5, μ = −0.53).
pub fn taubin(mesh: &mut Mesh, iterations: usize, lambda: f64, mu: f64, constrain: &dyn Fn(usize, [f64; 3]) -> Option<[f64; 3]>) {
    let nb = neighbours(mesh);
    let mut next = mesh.positions.clone();
    for _ in 0..iterations {
        for &factor in &[lambda, mu] {
            for v in 0..mesh.positions.len() {
                let list = &nb[v];
                if list.is_empty() {
                    continue;
                }
                let p = mesh.positions[v];
                let mut c = [0.0; 3];
                for &u in list {
                    let q = mesh.positions[u as usize];
                    for d in 0..3 {
                        c[d] += q[d];
                    }
                }
                let inv = 1.0 / list.len() as f64;
                let moved = [
                    p[0] + factor * (c[0] * inv - p[0]),
                    p[1] + factor * (c[1] * inv - p[1]),
                    p[2] + factor * (c[2] * inv - p[2]),
                ];
                next[v] = constrain(v, moved).unwrap_or(moved);
            }
            std::mem::swap(&mut mesh.positions, &mut next);
        }
    }
}

/// Closest-point queries on a triangle mesh (the original part, from `model.mesh`).
pub struct Projector {
    mesh: TriMesh,
}

impl Projector {
    pub fn new(positions: &[[f64; 3]], triangles: &[[u32; 3]]) -> Result<Self, String> {
        let verts = positions.iter().map(|p| Vector::new(p[0], p[1], p[2])).collect();
        let mesh = TriMesh::new(verts, triangles.to_vec()).map_err(|e| format!("surface mesh: {e:?}"))?;
        Ok(Projector { mesh })
    }

    /// The closest point on the surface and the distance to it.
    pub fn project(&self, p: [f64; 3]) -> ([f64; 3], f64) {
        let q = self.mesh.project_local_point(Vector::new(p[0], p[1], p[2]), false).point;
        let c = [q.x, q.y, q.z];
        let d = ((c[0] - p[0]).powi(2) + (c[1] - p[1]).powi(2) + (c[2] - p[2]).powi(2)).sqrt();
        (c, d)
    }
}
