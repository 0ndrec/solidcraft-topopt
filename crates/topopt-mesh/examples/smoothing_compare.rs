//! Compares the smoothers and pinning modes on shapes with known answers.
//!
//!   cargo run --release -p topopt-mesh --example smoothing_compare
//!
//! - A solid box with its own surface for pinning: how far the true edges and corners are from
//!   the output surface (sharp-edge pinning on and off).
//! - A free sphere and torus (no pinning): volume error and RMS distance to the true surface for
//!   each smoother, with and without volume preservation.

use topopt_mesh::{Field, Mesh, Options, Projector, Smoother, extract};

fn cells(dims: [usize; 3], h: f64, inside: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
    let mut c = Vec::with_capacity(dims[0] * dims[1] * dims[2]);
    for k in 0..dims[2] {
        for j in 0..dims[1] {
            for i in 0..dims[0] {
                // Fraction of 4³ sub-samples inside: a smooth density like an optimizer's.
                let mut n = 0;
                for s in 0..64 {
                    let o = [(s % 4) as f64, ((s / 4) % 4) as f64, (s / 16) as f64].map(|x| (x + 0.5) / 4.0);
                    if inside([(i as f64 + o[0]) * h, (j as f64 + o[1]) * h, (k as f64 + o[2]) * h]) <= 0.0 {
                        n += 1;
                    }
                }
                c.push(n as f64 / 64.0);
            }
        }
    }
    c
}

fn box_mesh(lo: [f64; 3], hi: [f64; 3]) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let c = |i: usize| [if i & 1 == 1 { hi[0] } else { lo[0] }, if i & 2 == 2 { hi[1] } else { lo[1] }, if i & 4 == 4 { hi[2] } else { lo[2] }];
    (
        vec![c(0), c(1), c(3), c(2), c(4), c(5), c(7), c(6)],
        vec![[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]],
    )
}

fn surface_distance(m: &Mesh, pts: &[[f64; 3]]) -> (f64, f64) {
    let pr = Projector::with_feature_angle(&m.positions, &m.triangles, 0.0).unwrap();
    let d: Vec<f64> = pts.iter().map(|&p| pr.project(p).1).collect();
    (d.iter().sum::<f64>() / d.len() as f64, d.iter().copied().fold(0.0, f64::max))
}

fn main() {
    let h = 1.0;
    // 1. A box, its edges and corners.
    let (lo, hi) = ([0.0; 3], [20.0, 12.0, 8.0]);
    let dims = [20, 12, 8];
    let c = vec![1.0; dims[0] * dims[1] * dims[2]];
    let (bp, bt) = box_mesh(lo, hi);
    let proj = Projector::new(&bp, &bt).unwrap();
    let mut edge_pts = Vec::new();
    for t in 0..=40 {
        let f = t as f64 / 40.0;
        for &(y, z) in &[(lo[1], lo[2]), (hi[1], lo[2]), (lo[1], hi[2]), (hi[1], hi[2])] {
            edge_pts.push([lo[0] + f * (hi[0] - lo[0]), y, z]);
        }
    }
    let exact = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
    println!("box 20×12×8 (h = 1): distance from the true edges and corners to the result");
    println!("{:<34} {:>10} {:>10} {:>10} {:>9} {:>6} {:>7}", "", "edge mean", "edge max", "corner max", "volume", "edge v", "corner");
    for smoother in [Smoother::Taubin, Smoother::Hc, Smoother::Bilateral] {
        for sharp in [false, true] {
            let f = Field::from_cells(dims, h, lo, &c, 1);
            let opts = Options { smoother, sharp_edges: sharp, target_triangles: None, ..Default::default() };
            let (m, rep) = extract(f, None, Some(&proj), &opts).unwrap();
            let (em, ex) = surface_distance(&m, &edge_pts);
            let (_, cx) = surface_distance(&m, &bp);
            println!(
                "{:<34} {:>10.4} {:>10.4} {:>10.4} {:>+8.2}% {:>6} {:>7}",
                format!("{smoother:?}, sharp edges {}", if sharp { "on" } else { "off" }),
                em,
                ex,
                cx,
                100.0 * (rep.mesh.volume - exact) / exact,
                rep.edge_vertices,
                rep.corner_vertices
            );
        }
    }

    // 2. Free shapes: volume and accuracy.
    let r = 6.3;
    let sphere = |p: [f64; 3]| ((p[0] - 10.0).powi(2) + (p[1] - 10.0).powi(2) + (p[2] - 10.0).powi(2)).sqrt() - r;
    let (big, small) = (7.0, 2.6);
    let torus = |p: [f64; 3]| {
        let q = ((p[0] - 12.0).powi(2) + (p[1] - 12.0).powi(2)).sqrt() - big;
        (q * q + (p[2] - 5.0).powi(2)).sqrt() - small
    };
    let shapes: [(&str, [usize; 3], &dyn Fn([f64; 3]) -> f64, f64); 2] = [
        ("sphere r 6.3", [20, 20, 20], &sphere, 4.0 / 3.0 * std::f64::consts::PI * r.powi(3)),
        ("torus 7 / 2.6", [24, 24, 10], &torus, 2.0 * std::f64::consts::PI.powi(2) * big * small * small),
    ];
    for (name, dims, sdf, exact) in shapes {
        println!("\n{name}: volume error and RMS distance to the true surface (in h)");
        let c = cells(dims, h, sdf);
        for smoother in [Smoother::Taubin, Smoother::Hc, Smoother::Bilateral] {
            for pv in [false, true] {
                let f = Field::from_cells(dims, h, [0.0; 3], &c, 1);
                let opts = Options { smoother, preserve_volume: pv, target_triangles: None, ..Default::default() };
                let (m, rep) = extract(f, None, None, &opts).unwrap();
                let rms = (m.positions.iter().map(|&p| sdf(p).powi(2)).sum::<f64>() / m.positions.len() as f64).sqrt();
                println!(
                    "  {:<34} volume {:>+6.2}%   rms {:.4}   (contour {:+.2}%)",
                    format!("{smoother:?}, preserve volume {}", if pv { "on" } else { "off" }),
                    100.0 * (rep.mesh.volume - exact) / exact,
                    rms / h,
                    100.0 * (rep.contour_volume - exact) / exact
                );
            }
        }
    }
}
