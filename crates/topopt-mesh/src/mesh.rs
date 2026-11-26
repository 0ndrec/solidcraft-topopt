//! Indexed triangle mesh, its checks and binary STL output.

use std::collections::HashMap;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    pub positions: Vec<[f64; 3]>,
    pub triangles: Vec<[u32; 3]>,
}

/// What [`Mesh::validate`] found. `ok()` means: closed, 2-manifold, consistently oriented with
/// outward normals (positive volume), no degenerate triangles and, when checked, no
/// self-intersections.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MeshReport {
    pub vertices: usize,
    pub triangles: usize,
    pub edges: usize,
    /// Edges used by one triangle.
    pub boundary_edges: usize,
    /// Edges used by more than two triangles.
    pub nonmanifold_edges: usize,
    /// Edges whose two triangles run them in the same direction.
    pub misoriented_edges: usize,
    /// Vertices whose triangles form more than one fan (pinch points).
    pub nonmanifold_vertices: usize,
    pub unused_vertices: usize,
    pub degenerate_triangles: usize,
    /// Connected components (by shared edges).
    pub components: usize,
    /// V − E + F.
    pub euler: i64,
    /// Total genus Σ(2 − χᵢ)/2 over components (holes through the part).
    pub genus: i64,
    pub volume: f64,
    pub area: f64,
    /// `None` when not checked.
    pub self_intersections: Option<usize>,
    /// Distinct vertices that become the same point once rounded to f32 (STL), which a reader
    /// that welds equal corners would merge.
    pub f32_collisions: usize,
}

impl MeshReport {
    pub fn ok(&self) -> bool {
        self.triangles > 0
            && self.boundary_edges == 0
            && self.nonmanifold_edges == 0
            && self.misoriented_edges == 0
            && self.nonmanifold_vertices == 0
            && self.degenerate_triangles == 0
            && self.volume > 0.0
            && self.self_intersections.unwrap_or(0) == 0
            && self.f32_collisions == 0
    }
}

pub(crate) fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub(crate) fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Area-weighted vertex normals (unit; zero for unused vertices).
pub fn vertex_normals(m: &Mesh) -> Vec<[f64; 3]> {
    let mut n = vec![[0.0; 3]; m.positions.len()];
    for t in 0..m.triangles.len() {
        let f = m.face_normal(t);
        for &v in &m.triangles[t] {
            for d in 0..3 {
                n[v as usize][d] += f[d];
            }
        }
    }
    for v in &mut n {
        let l = norm(*v);
        if l > 0.0 {
            *v = [v[0] / l, v[1] / l, v[2] / l];
        }
    }
    n
}

/// Triangles with (near) zero area, relative to the mean edge length.
pub fn degenerate_triangles(m: &Mesh) -> Vec<u32> {
    let scale = m.mean_edge().max(1e-300);
    (0..m.triangles.len())
        .filter(|&t| {
            let tri = m.triangles[t];
            norm(m.face_normal(t)) / 2.0 <= 1e-12 * scale * scale || tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2]
        })
        .map(|t| t as u32)
        .collect()
}

impl Mesh {
    /// Twice-area normal of triangle `t`.
    pub fn face_normal(&self, t: usize) -> [f64; 3] {
        let [a, b, c] = self.triangles[t].map(|i| self.positions[i as usize]);
        cross(sub(b, a), sub(c, a))
    }

    pub fn volume(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| self.positions[i as usize]);
                dot(a, cross(b, c)) / 6.0
            })
            .sum()
    }

    pub fn area(&self) -> f64 {
        (0..self.triangles.len()).map(|t| norm(self.face_normal(t)) / 2.0).sum()
    }

    /// Mean edge length.
    pub fn mean_edge(&self) -> f64 {
        let mut sum = 0.0;
        for t in &self.triangles {
            for e in 0..3 {
                sum += norm(sub(self.positions[t[e] as usize], self.positions[t[(e + 1) % 3] as usize]));
            }
        }
        sum / (3 * self.triangles.len()).max(1) as f64
    }

    /// Drops vertices no triangle uses and renumbers.
    pub fn compact(&mut self) {
        let mut map = vec![u32::MAX; self.positions.len()];
        let mut positions = Vec::new();
        for t in &mut self.triangles {
            for v in t.iter_mut() {
                if map[*v as usize] == u32::MAX {
                    map[*v as usize] = positions.len() as u32;
                    positions.push(self.positions[*v as usize]);
                }
                *v = map[*v as usize];
            }
        }
        self.positions = positions;
    }

    /// Component id of every triangle (triangles joined by shared edges) and the count.
    pub fn components(&self) -> (Vec<u32>, usize) {
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
            for e in 0..3 {
                let (a, b) = (tri[e], tri[(e + 1) % 3]);
                let key = (a.min(b), a.max(b));
                match first.get(&key) {
                    Some(&o) => {
                        let (ra, rb) = (find(&mut parent, o), find(&mut parent, t as u32));
                        parent[ra as usize] = rb;
                    }
                    None => {
                        first.insert(key, t as u32);
                    }
                }
            }
        }
        let mut ids: HashMap<u32, u32> = HashMap::new();
        let mut comp = vec![0u32; self.triangles.len()];
        for t in 0..self.triangles.len() {
            let r = find(&mut parent, t as u32);
            let n = ids.len() as u32;
            comp[t] = *ids.entry(r).or_insert(n);
        }
        (comp, ids.len())
    }

    /// Full topological and geometric check. `self_intersections`: also test triangle pairs
    /// (spatial hashing + segment/triangle tests); costs more than the rest together.
    pub fn validate(&self, self_intersections: bool) -> MeshReport {
        let mut r = MeshReport { vertices: self.positions.len(), triangles: self.triangles.len(), ..Default::default() };
        // Edges: (uses, orientation sum).
        let mut edges: HashMap<(u32, u32), (u32, i32)> = HashMap::with_capacity(self.triangles.len() * 2);
        let mut used = vec![false; self.positions.len()];
        let scale = self.mean_edge().max(1e-300);
        for (t, tri) in self.triangles.iter().enumerate() {
            for e in 0..3 {
                let (a, b) = (tri[e], tri[(e + 1) % 3]);
                used[a as usize] = true;
                let entry = edges.entry((a.min(b), a.max(b))).or_insert((0, 0));
                entry.0 += 1;
                entry.1 += if a < b { 1 } else { -1 };
            }
            if norm(self.face_normal(t)) / 2.0 <= 1e-12 * scale * scale || tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                r.degenerate_triangles += 1;
            }
        }
        r.edges = edges.len();
        for &(uses, orient) in edges.values() {
            match uses {
                1 => r.boundary_edges += 1,
                2 => {
                    if orient != 0 {
                        r.misoriented_edges += 1
                    }
                }
                _ => r.nonmanifold_edges += 1,
            }
        }
        r.unused_vertices = used.iter().filter(|&&u| !u).count();
        r.nonmanifold_vertices = self.count_pinched_vertices();
        let (comp, ncomp) = self.components();
        r.components = ncomp;
        r.euler = (r.vertices - r.unused_vertices) as i64 - r.edges as i64 + r.triangles as i64;
        // Genus per component from its own V, E, F.
        let mut cv: Vec<std::collections::HashSet<u32>> = vec![Default::default(); ncomp];
        let mut ce = vec![0i64; ncomp];
        let mut cf = vec![0i64; ncomp];
        for (t, tri) in self.triangles.iter().enumerate() {
            let c = comp[t] as usize;
            cf[c] += 1;
            for &v in tri {
                cv[c].insert(v);
            }
        }
        // Edges per component: in a closed manifold each edge is shared by two triangles of one
        // component, so E = 3F/2 per component.
        for (t, _) in self.triangles.iter().enumerate() {
            ce[comp[t] as usize] += 3;
        }
        r.genus = (0..ncomp).map(|c| (2 - (cv[c].len() as i64 - ce[c] / 2 + cf[c])) / 2).sum();
        r.volume = self.volume();
        r.area = self.area();
        let mut seen = std::collections::HashSet::with_capacity(self.positions.len());
        for (v, p) in self.positions.iter().enumerate() {
            if used[v] && !seen.insert(p.map(|c| (c as f32).to_bits())) {
                r.f32_collisions += 1;
            }
        }
        if self_intersections {
            r.self_intersections = Some(crate::intersect::count_self_intersections(self));
        }
        r
    }

    /// Vertices whose incident triangles do not form a single fan.
    fn count_pinched_vertices(&self) -> usize {
        // For each vertex: map "next" around the vertex from each incident triangle (v, a, b): a → b.
        let mut next: Vec<Vec<(u32, u32)>> = vec![Vec::new(); self.positions.len()];
        for tri in &self.triangles {
            for e in 0..3 {
                let (v, a, b) = (tri[e], tri[(e + 1) % 3], tri[(e + 2) % 3]);
                next[v as usize].push((a, b));
            }
        }
        let mut pinched = 0;
        for fan in &next {
            if fan.is_empty() {
                continue;
            }
            let map: HashMap<u32, u32> = fan.iter().copied().collect();
            if map.len() != fan.len() {
                pinched += 1; // the same spoke starts two triangles: not a disc
                continue;
            }
            // Follow the cycle from the first spoke; a single fan visits every triangle.
            let start = fan[0].0;
            let mut cur = start;
            let mut steps = 0;
            loop {
                match map.get(&cur) {
                    Some(&n) => {
                        cur = n;
                        steps += 1;
                    }
                    None => break,
                }
                if cur == start || steps > fan.len() {
                    break;
                }
            }
            if steps != fan.len() || cur != start {
                pinched += 1;
            }
        }
        pinched
    }

    /// Indexed binary PLY (mm): 12 bytes per vertex and 13 per triangle, about 19 bytes per
    /// triangle for a closed mesh against STL's 50, and corners are shared by index, so nothing
    /// depends on a reader welding equal coordinates.
    pub fn to_ply(&self, comment: &str) -> Vec<u8> {
        let header = format!(
            "ply\nformat binary_little_endian 1.0\ncomment {}\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nelement face {}\nproperty list uchar int vertex_indices\nend_header\n",
            comment.replace('\n', " "),
            self.positions.len(),
            self.triangles.len()
        );
        let mut out = Vec::with_capacity(header.len() + 12 * self.positions.len() + 13 * self.triangles.len());
        out.extend_from_slice(header.as_bytes());
        for p in &self.positions {
            for c in p {
                out.extend_from_slice(&(*c as f32).to_le_bytes());
            }
        }
        for t in &self.triangles {
            out.push(3);
            for &v in t {
                out.extend_from_slice(&(v as i32).to_le_bytes());
            }
        }
        out
    }

    /// Binary STL (mm), facet normals from the geometry.
    pub fn to_stl(&self, header: &str) -> Vec<u8> {
        let mut out = Vec::with_capacity(84 + 50 * self.triangles.len());
        let mut h = [0u8; 80];
        for (d, s) in h.iter_mut().zip(header.bytes()) {
            *d = s;
        }
        out.extend_from_slice(&h);
        out.extend_from_slice(&(self.triangles.len() as u32).to_le_bytes());
        for t in 0..self.triangles.len() {
            let n = self.face_normal(t);
            let l = norm(n).max(1e-300);
            for c in n {
                out.extend_from_slice(&((c / l) as f32).to_le_bytes());
            }
            for &v in &self.triangles[t] {
                for c in self.positions[v as usize] {
                    out.extend_from_slice(&(c as f32).to_le_bytes());
                }
            }
            out.extend_from_slice(&[0, 0]);
        }
        out
    }
}
