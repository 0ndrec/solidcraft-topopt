//! Reading triangle meshes: binary and ASCII STL, and the binary PLY this crate writes.
//! Corners at exactly the same coordinates are merged, so a closed STL gives a closed indexed mesh.

use crate::mesh::Mesh;
use std::collections::HashMap;

/// Merges equal corners (bit-exact f32) into an indexed mesh.
struct Welder {
    index: HashMap<[u32; 3], u32>,
    mesh: Mesh,
}

impl Welder {
    fn new() -> Self {
        Welder { index: HashMap::new(), mesh: Mesh::default() }
    }
    fn vertex(&mut self, p: [f32; 3]) -> u32 {
        let key = p.map(f32::to_bits);
        *self.index.entry(key).or_insert_with(|| {
            self.mesh.positions.push(p.map(|c| c as f64));
            (self.mesh.positions.len() - 1) as u32
        })
    }
    fn triangle(&mut self, t: [[f32; 3]; 3]) {
        let ids = t.map(|p| self.vertex(p));
        // Triangles that collapse when welded carry no area.
        if ids[0] != ids[1] && ids[1] != ids[2] && ids[0] != ids[2] {
            self.mesh.triangles.push(ids);
        }
    }
}

/// Reads an STL file's bytes, binary or ASCII. `scale` converts its units to millimetres.
pub fn read_stl(bytes: &[u8], scale: f64) -> Result<Mesh, String> {
    let mut w = Welder::new();
    let binary_len = |n: usize| 84 + 50 * n;
    let is_binary = bytes.len() >= 84 && {
        let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        bytes.len() == binary_len(n)
    };
    if is_binary {
        let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        for t in 0..n {
            let base = 84 + 50 * t + 12;
            let f = |k: usize| f32::from_le_bytes(bytes[base + 4 * k..base + 4 * k + 4].try_into().unwrap());
            w.triangle([[f(0), f(1), f(2)], [f(3), f(4), f(5)], [f(6), f(7), f(8)]]);
        }
    } else {
        let text = std::str::from_utf8(bytes).map_err(|_| "not a binary STL (wrong size) and not ASCII text".to_string())?;
        if !text.trim_start().starts_with("solid") {
            return Err("not an STL file".into());
        }
        let mut corners: Vec<[f32; 3]> = Vec::with_capacity(3);
        for (line_no, line) in text.lines().enumerate() {
            let mut it = line.split_whitespace();
            if it.next() == Some("vertex") {
                let mut p = [0f32; 3];
                for c in &mut p {
                    *c = it.next().and_then(|s| s.parse().ok()).ok_or_else(|| format!("line {}: bad vertex", line_no + 1))?;
                }
                corners.push(p);
                if corners.len() == 3 {
                    w.triangle([corners[0], corners[1], corners[2]]);
                    corners.clear();
                }
            }
        }
    }
    if w.mesh.triangles.is_empty() {
        return Err("the STL has no triangles".into());
    }
    if scale != 1.0 {
        for p in &mut w.mesh.positions {
            *p = p.map(|c| c * scale);
        }
    }
    Ok(w.mesh)
}

/// Reads the binary little-endian PLY that [`Mesh::to_ply`] writes (vertex x y z float, face list
/// uchar int). Other PLY layouts are refused.
pub fn read_ply(bytes: &[u8]) -> Result<Mesh, String> {
    let end = bytes.windows(11).position(|w| w == b"end_header\n").ok_or("not a PLY file")? + 11;
    let header = std::str::from_utf8(&bytes[..end]).map_err(|_| "bad PLY header")?;
    if !header.contains("format binary_little_endian 1.0") {
        return Err("only binary little-endian PLY is read".into());
    }
    let count = |name: &str| -> Result<usize, String> {
        header
            .lines()
            .find_map(|l| l.strip_prefix(&format!("element {name} ")))
            .and_then(|n| n.trim().parse().ok())
            .ok_or(format!("PLY without {name} count"))
    };
    let (nv, nf) = (count("vertex")?, count("face")?);
    let mut at = end;
    let mut mesh = Mesh::default();
    let f32_at = |i: usize| -> Result<f32, String> { Ok(f32::from_le_bytes(bytes.get(i..i + 4).ok_or("truncated PLY")?.try_into().unwrap())) };
    for _ in 0..nv {
        mesh.positions.push([f32_at(at)? as f64, f32_at(at + 4)? as f64, f32_at(at + 8)? as f64]);
        at += 12;
    }
    for _ in 0..nf {
        if *bytes.get(at).ok_or("truncated PLY")? != 3 {
            return Err("only triangle faces are read".into());
        }
        let i = |k: usize| -> Result<u32, String> {
            let v = i32::from_le_bytes(bytes.get(at + 1 + 4 * k..at + 5 + 4 * k).ok_or("truncated PLY")?.try_into().unwrap());
            if v < 0 || v as usize >= nv { Err("PLY face index out of range".into()) } else { Ok(v as u32) }
        };
        mesh.triangles.push([i(0)?, i(1)?, i(2)?]);
        at += 13;
    }
    Ok(mesh)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tetra() -> Mesh {
        Mesh {
            positions: vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            triangles: vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
        }
    }

    #[test]
    fn stl_and_ply_round_trips() {
        let m = tetra();
        assert!(m.validate(true).ok());
        let from_stl = read_stl(&m.to_stl("t"), 1.0).unwrap();
        assert_eq!(from_stl.positions.len(), 4, "corners are welded");
        assert!(from_stl.validate(true).ok());
        assert!((from_stl.volume() - m.volume()).abs() < 1e-9);
        let from_ply = read_ply(&m.to_ply("t")).unwrap();
        assert_eq!(from_ply, m);
        // ASCII, in centimetres.
        let mut ascii = String::from("solid t\n");
        for t in &m.triangles {
            ascii.push_str("facet normal 0 0 0\nouter loop\n");
            for &v in t {
                let p = m.positions[v as usize];
                ascii.push_str(&format!("vertex {} {} {}\n", p[0] / 10.0, p[1] / 10.0, p[2] / 10.0));
            }
            ascii.push_str("endloop\nendfacet\n");
        }
        ascii.push_str("endsolid t\n");
        let cm = read_stl(ascii.as_bytes(), 10.0).unwrap();
        assert!((cm.volume() - m.volume()).abs() < 1e-6);
        assert!(read_stl(b"garbage", 1.0).is_err());
    }
}
