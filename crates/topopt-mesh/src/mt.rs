//! Marching tetrahedra on the Freudenthal (Kuhn) decomposition.
//!
//! Every grid cube is split into six tetrahedra around its main diagonal 0–7. The split is the
//! same in every cube, so neighbouring cubes' face diagonals match and the tetrahedra form a
//! conforming mesh. Inside a tetrahedron the linear interpolant's zero set is a plane, so the
//! surface piece there is one planar triangle or one planar convex quad. Hence:
//!
//! - one vertex per crossed grid edge, shared through an edge-keyed map: the mesh is indexed and
//!   **watertight**;
//! - with no sample exactly at the iso value ([`crate::field::Field::nudge`]), the surface is a
//!   closed **2-manifold** and **does not self-intersect** (pieces in different tetrahedra meet
//!   only on shared faces);
//! - orientation comes from a table built once from exact corner positions: every triangle's
//!   normal points from `v > iso` to `v < iso` (outward from the solid).

use crate::field::Field;
use crate::mesh::Mesh;
use std::collections::HashMap;

/// Cube corner `c` is at `(c & 1, (c >> 1) & 1, (c >> 2) & 1)`.
const TETS: [[usize; 4]; 6] = [[0, 1, 3, 7], [0, 1, 5, 7], [0, 2, 3, 7], [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 6, 7]];

fn corner(c: usize) -> [f64; 3] {
    [(c & 1) as f64, ((c >> 1) & 1) as f64, ((c >> 2) & 1) as f64]
}

/// For each tetrahedron and each inside-mask (bit `v` set when tet vertex `v` is inside), the
/// polygon as a list of tet-local edges `(a, b)`, oriented outward.
type Table = Vec<Vec<Vec<(usize, usize)>>>;

fn build_table() -> Table {
    let mut table = vec![vec![Vec::new(); 16]; 6];
    for (t, tet) in TETS.iter().enumerate() {
        for mask in 1..15usize {
            let inside: Vec<usize> = (0..4).filter(|v| mask >> v & 1 == 1).collect();
            let outside: Vec<usize> = (0..4).filter(|v| mask >> v & 1 == 0).collect();
            let mut poly: Vec<(usize, usize)> = match inside.len() {
                1 => outside.iter().map(|&o| (inside[0], o)).collect(),
                3 => inside.iter().map(|&i| (i, outside[0])).collect(),
                2 => {
                    let (a, b, c, d) = (inside[0], inside[1], outside[0], outside[1]);
                    vec![(a, c), (a, d), (b, d), (b, c)]
                }
                _ => unreachable!(),
            };
            // Orientation from midpoints: the polygon is planar for any crossing positions, so its
            // orientation does not depend on where along the edges the vertices lie.
            let p: Vec<[f64; 3]> = poly
                .iter()
                .map(|&(a, b)| {
                    let (pa, pb) = (corner(tet[a]), corner(tet[b]));
                    [(pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0, (pa[2] + pb[2]) / 2.0]
                })
                .collect();
            let mut normal = [0.0; 3];
            for i in 0..p.len() {
                let (u, w) = (p[i], p[(i + 1) % p.len()]);
                normal[0] += (u[1] - w[1]) * (u[2] + w[2]);
                normal[1] += (u[2] - w[2]) * (u[0] + w[0]);
                normal[2] += (u[0] - w[0]) * (u[1] + w[1]);
            }
            let mean = |vs: &[usize]| {
                let mut m = [0.0; 3];
                for &v in vs {
                    let c = corner(tet[v]);
                    for d in 0..3 {
                        m[d] += c[d] / vs.len() as f64;
                    }
                }
                m
            };
            let (mi, mo) = (mean(&inside), mean(&outside));
            let out = [mo[0] - mi[0], mo[1] - mi[1], mo[2] - mi[2]];
            if normal[0] * out[0] + normal[1] * out[1] + normal[2] * out[2] < 0.0 {
                poly.reverse();
            }
            table[t][mask] = poly;
        }
    }
    table
}

/// The iso-surface `{v = iso}` of the field, oriented with normals towards smaller values.
pub fn marching_tetrahedra(field: &Field, iso: f64) -> Mesh {
    let table = build_table();
    let n = field.n;
    let v = &field.values;
    let stride = [1usize, n[0], n[0] * n[1]];
    let corner_offset: Vec<usize> = (0..8).map(|c| (c & 1) * stride[0] + ((c >> 1) & 1) * stride[1] + ((c >> 2) & 1) * stride[2]).collect();
    let mut mesh = Mesh::default();
    let mut edge_vertex: HashMap<(u32, u32), u32> = HashMap::new();
    let mut poly_ids = Vec::with_capacity(4);
    for k in 0..n[2] - 1 {
        for j in 0..n[1] - 1 {
            for i in 0..n[0] - 1 {
                let base = field.idx(i, j, k);
                let cube: [usize; 8] = std::array::from_fn(|c| base + corner_offset[c]);
                let inside: [bool; 8] = std::array::from_fn(|c| v[cube[c]] > iso);
                if inside.iter().all(|&x| x) || inside.iter().all(|&x| !x) {
                    continue;
                }
                for (t, tet) in TETS.iter().enumerate() {
                    let mask = (0..4).fold(0usize, |m, q| m | ((inside[tet[q]] as usize) << q));
                    if mask == 0 || mask == 15 {
                        continue;
                    }
                    poly_ids.clear();
                    for &(a, b) in &table[t][mask] {
                        let (ga, gb) = (cube[tet[a]], cube[tet[b]]);
                        let key = if ga < gb { (ga as u32, gb as u32) } else { (gb as u32, ga as u32) };
                        let id = *edge_vertex.entry(key).or_insert_with(|| {
                            let (va, vb) = (v[ga], v[gb]);
                            let s = (iso - va) / (vb - va);
                            let (pa, pb) = (field.pos_of(ga), field.pos_of(gb));
                            mesh.positions.push([pa[0] + s * (pb[0] - pa[0]), pa[1] + s * (pb[1] - pa[1]), pa[2] + s * (pb[2] - pa[2])]);
                            (mesh.positions.len() - 1) as u32
                        });
                        poly_ids.push(id);
                    }
                    mesh.triangles.push([poly_ids[0], poly_ids[1], poly_ids[2]]);
                    if poly_ids.len() == 4 {
                        mesh.triangles.push([poly_ids[0], poly_ids[2], poly_ids[3]]);
                    }
                }
            }
        }
    }
    mesh
}

impl Field {
    /// Position of sample index `s`.
    #[inline]
    pub fn pos_of(&self, s: usize) -> [f64; 3] {
        let [i, j, k] = self.ijk(s);
        self.pos(i, j, k)
    }
}
