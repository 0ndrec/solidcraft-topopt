//! Meshes a field given as JSON, for experiments with fields prepared elsewhere (filtered,
//! upsampled, distance fields…).
//!
//!   cargo run --release -p topopt-mesh --example field_mesh -- FIELD.json OUT.ply
//!
//! FIELD.json: {"dims": [nx, ny, nz], "h": h, "origin": [x, y, z], "values": [cell values,
//! index i + nx·(j + ny·k)], "anchors": [bool per cell] (optional), "box": [[lo], [hi]] (optional
//! original surface), "smoother": "taubin" | "hc" | "bilateral",
//! "iterations": n, "iso": 0.5, "triangles": n, "sharp_edges": true}

use serde_json::Value;
use topopt_mesh::{Field, Options, Projector, Smoother, extract};

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let v: Value = serde_json::from_slice(&std::fs::read(&a[0]).unwrap()).unwrap();
    let dims: [usize; 3] = std::array::from_fn(|d| v["dims"][d].as_u64().unwrap() as usize);
    let origin: [f64; 3] = std::array::from_fn(|d| v["origin"][d].as_f64().unwrap());
    let h = v["h"].as_f64().unwrap();
    let values: Vec<f64> = v["values"].as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect();
    let field = Field::from_cells(dims, h, origin, &values, 1);
    let anchors = v["anchors"].as_array().map(|l| field.cell_mask(dims, 1, &l.iter().map(|x| x.as_bool().unwrap()).collect::<Vec<_>>()));
    let surface = v["box"].as_array().map(|b| {
        let lo: [f64; 3] = std::array::from_fn(|d| b[0][d].as_f64().unwrap());
        let hi: [f64; 3] = std::array::from_fn(|d| b[1][d].as_f64().unwrap());
        let c = |i: usize| [if i & 1 == 1 { hi[0] } else { lo[0] }, if i & 2 == 2 { hi[1] } else { lo[1] }, if i & 4 == 4 { hi[2] } else { lo[2] }];
        let p = vec![c(0), c(1), c(3), c(2), c(4), c(5), c(7), c(6)];
        let t = vec![[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]];
        Projector::new(&p, &t).unwrap()
    });
    let opts = Options {
        smoother: Smoother::parse(v["smoother"].as_str().unwrap_or("taubin")).unwrap(),
        smooth_iterations: v["iterations"].as_u64().unwrap_or(10) as usize,
        iso: v["iso"].as_f64().unwrap_or(0.5),
        target_triangles: Some(v["triangles"].as_u64().unwrap_or(30_000) as usize),
        sharp_edges: v["sharp_edges"].as_bool().unwrap_or(true),
        ..Default::default()
    };
    let t = std::time::Instant::now();
    let (mesh, rep) = extract(field, anchors.as_deref(), surface.as_ref(), &opts).unwrap();
    std::fs::write(&a[1], mesh.to_ply("field_mesh")).unwrap();
    println!(
        "{{\"triangles\": {}, \"volume\": {:.1}, \"genus\": {}, \"smooth_iterations\": {}, \"sharp_edges\": {}, \"seconds\": {:.2}}}",
        rep.mesh.triangles,
        rep.mesh.volume,
        rep.mesh.genus,
        rep.smooth_iterations,
        rep.sharp_edges,
        t.elapsed().as_secs_f64()
    );
}
