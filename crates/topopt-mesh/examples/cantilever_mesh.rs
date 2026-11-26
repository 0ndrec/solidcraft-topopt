//! End to end: optimize the top3d cantilever with topopt-core, turn the densities into a mesh,
//! write a binary STL.
//!
//!   cargo run --release -p topopt-mesh --example cantilever_mesh -- OUT.stl [NX NY NZ] [H_mm] [TRIANGLES]

use std::time::Instant;
use topopt_core::{Grid, Optimizer, Problem, Settings};
use topopt_mesh::{Field, Options, Projector, extract};

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let out = a.first().cloned().unwrap_or_else(|| "cantilever.stl".into());
    let num = |i: usize, d: f64| a.get(i).and_then(|s| s.parse::<f64>().ok()).unwrap_or(d);
    let (nx, ny, nz) = (num(1, 48.0) as usize, num(2, 24.0) as usize, num(3, 24.0) as usize);
    let h = num(4, 2.0);
    let target = num(5, 30000.0) as usize;
    let (l, w, t) = (nx as f64 * h, ny as f64 * h, nz as f64 * h);

    let mut pr = Problem::new(Grid::new(nx, ny, nz, h));
    pr.material.e = 70_000.0;
    pr.material.nu = 0.33;
    pr.fix_box([0.0, 0.0, 0.0], [0.0, w, t], [true; 3]);
    pr.load_box([l, 0.0, 0.0], [l, w, 0.0], [0.0, 0.0, -1000.0]);
    // Solid pads where the part is held and loaded: anchors for clean-up, and real mounting faces.
    pr.set_kind_box([0.0, 0.0, 0.0], [h, w, t], topopt_core::ElementKind::Solid);
    pr.set_kind_box([l - h, 0.0, 0.0], [l, w, 2.0 * h], topopt_core::ElementKind::Solid);
    let settings = Settings { volfrac: 0.2, ..Default::default() };

    let t0 = Instant::now();
    let mut opt = Optimizer::new(pr.clone(), settings).unwrap();
    let outcome = opt.run(|it| {
        if it.iter % 20 == 0 {
            eprintln!("iter {:3}  c {:.4}  vol {:.3}  beta {}  gray {:.3}", it.iter, it.compliance, it.volume, it.beta, it.grayness);
        }
        true
    });
    let last = *opt.history.last().unwrap();
    eprintln!("{outcome:?}: {} iterations, {:.1} s, compliance {:.4} N·mm, grayness {:.4}", last.iter, t0.elapsed().as_secs_f64(), last.compliance, last.grayness);

    let t1 = Instant::now();
    let dims = [nx, ny, nz];
    let pad = 1;
    let field = Field::from_cells(dims, h, [0.0; 3], &opt.x_phys, pad);
    let solid: Vec<bool> = (0..pr.grid.nel()).map(|e| pr.kind(e) == topopt_core::ElementKind::Solid).collect();
    let anchors = field.cell_mask(dims, pad, &solid);
    // The design domain's own surface (a box), for pinning mounting faces.
    let bx = [[0.0, 0.0, 0.0], [l, 0.0, 0.0], [l, w, 0.0], [0.0, w, 0.0], [0.0, 0.0, t], [l, 0.0, t], [l, w, t], [0.0, w, t]];
    let bt: Vec<[u32; 3]> = vec![[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]];
    let box_surface = Projector::new(&bx, &bt).unwrap();
    let opts = Options { target_triangles: Some(target), ..Default::default() };
    let (mesh, rep) = extract(field, Some(&anchors), Some(&box_surface), &opts).expect("mesh");
    eprintln!("mesh in {:.2} s: {:?}", t1.elapsed().as_secs_f64(), rep);
    let design_volume = last.volume * (l * w * t);
    eprintln!(
        "material volume from densities ≈ {:.0} mm³ (design part only), mesh volume {:.0} mm³",
        design_volume, rep.mesh.volume
    );
    std::fs::write(&out, mesh.to_stl("topopt-mesh cantilever")).unwrap();
    eprintln!("wrote {out}: {} triangles", mesh.triangles.len());
}
