//! Mesh guarantees on analytic shapes and on random fields.

use topopt_mesh::{Field, Options, Projector, extract, marching_tetrahedra};

/// Cell densities of a solid given by a signed distance (negative inside), smoothed over one cell.
fn cells_from_sdf(dims: [usize; 3], h: f64, sdf: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
    let mut cells = Vec::with_capacity(dims[0] * dims[1] * dims[2]);
    for k in 0..dims[2] {
        for j in 0..dims[1] {
            for i in 0..dims[0] {
                let p = [(i as f64 + 0.5) * h, (j as f64 + 0.5) * h, (k as f64 + 0.5) * h];
                cells.push((0.5 - sdf(p) / h).clamp(0.0, 1.0));
            }
        }
    }
    cells
}

fn sphere_sdf(c: [f64; 3], r: f64) -> impl Fn([f64; 3]) -> f64 {
    move |p| ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt() - r
}

/// A UV sphere with exactly shared poles, outward.
fn uv_sphere(c: [f64; 3], r: f64, seg: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let rings = seg / 2;
    let mut pos = vec![[c[0], c[1], c[2] + r]];
    for ring in 1..rings {
        let th = std::f64::consts::PI * ring as f64 / rings as f64;
        for s in 0..seg {
            let ph = 2.0 * std::f64::consts::PI * s as f64 / seg as f64;
            pos.push([c[0] + r * th.sin() * ph.cos(), c[1] + r * th.sin() * ph.sin(), c[2] + r * th.cos()]);
        }
    }
    pos.push([c[0], c[1], c[2] - r]);
    let south = (pos.len() - 1) as u32;
    let at = |ring: usize, s: usize| (1 + (ring - 1) * seg + s % seg) as u32;
    let mut tri = Vec::new();
    for s in 0..seg {
        tri.push([0, at(1, s), at(1, s + 1)]);
        tri.push([south, at(rings - 1, s + 1), at(rings - 1, s)]);
    }
    for ring in 1..rings - 1 {
        for s in 0..seg {
            tri.push([at(ring, s), at(ring + 1, s), at(ring + 1, s + 1)]);
            tri.push([at(ring, s), at(ring + 1, s + 1), at(ring, s + 1)]);
        }
    }
    (pos, tri)
}

#[test]
fn sphere_contour_is_closed_manifold_genus_zero() {
    let (h, r) = (1.0, 7.3);
    let dims = [20, 20, 20];
    let cells = cells_from_sdf(dims, h, sphere_sdf([10.0, 10.0, 10.0], r));
    let mut f = Field::from_cells(dims, h, [0.0; 3], &cells, 1);
    f.nudge(0.5, 0.02);
    let m = marching_tetrahedra(&f, 0.5);
    let rep = m.validate(true);
    println!("{rep:?}");
    assert!(rep.ok(), "{rep:?}");
    assert_eq!(rep.components, 1);
    assert_eq!(rep.genus, 0);
    let exact = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
    assert!((rep.volume - exact).abs() / exact < 0.03, "volume {} vs {exact}", rep.volume);
}

#[test]
fn torus_has_genus_one() {
    let dims = [30, 30, 12];
    let (big, small) = (9.0, 3.2);
    let cells = cells_from_sdf(dims, 1.0, |p| {
        let (x, y, z) = (p[0] - 15.0, p[1] - 15.0, p[2] - 6.0);
        ((x * x + y * y).sqrt() - big).hypot(z) - small
    });
    let f = Field::from_cells(dims, 1.0, [0.0; 3], &cells, 1);
    let (m, rep) = extract(f, None, None, &Options { target_triangles: Some(1500), blur_passes: 0, ..Default::default() }).unwrap();
    let exact = 2.0 * std::f64::consts::PI.powi(2) * big * small * small;
    println!("torus volumes: exact {exact:.1}, contour {:.1}, smoothed {:.1}, final {:.1}", rep.contour_volume, rep.smoothed_volume, rep.mesh.volume);
    assert!(rep.mesh.ok());
    assert_eq!(rep.mesh.genus, 1);
    assert!((rep.mesh.volume - exact).abs() / exact < 0.03, "torus volume {} vs {exact}", rep.mesh.volume);
    assert_eq!(rep.mesh.components, 1);
    assert!(m.triangles.len() <= 1500 * 3 / 2);
}

/// Deterministic pseudo-random numbers in [0, 1).
fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    (*seed >> 11) as f64 / (1u64 << 53) as f64
}

#[test]
fn random_fields_always_give_valid_meshes() {
    for seed0 in 0..60u64 {
        let mut seed = seed0 * 7919 + 1;
        let n = 4 + (seed0 % 7) as usize;
        let dims = [n, n + 1, n + 2];
        let cells: Vec<f64> = (0..dims[0] * dims[1] * dims[2]).map(|_| lcg(&mut seed)).collect();
        // Raw contour, no clean-up: may have many shells, must still be valid.
        let mut f = Field::from_cells(dims, 1.0, [0.0; 3], &cells, 1);
        f.nudge(0.5, 0.02);
        let m = marching_tetrahedra(&f, 0.5);
        if !m.triangles.is_empty() {
            let rep = m.validate(true);
            assert!(rep.ok(), "seed {seed0}, raw contour: {rep:?}");
        }
        // Full pipeline: one shell, smoothed, decimated, still valid.
        let f = Field::from_cells(dims, 1.0, [0.0; 3], &cells, 1);
        let opts = Options { blur_passes: (seed0 % 2) as usize, target_triangles: Some(200), ..Default::default() };
        match extract(f, None, None, &opts) {
            Ok((_, rep)) => {
                assert!(rep.mesh.ok(), "seed {seed0}: {rep:?}");
                assert_eq!(rep.mesh.components, 1, "seed {seed0}: islands or cavities left: {rep:?}");
            }
            Err(e) => assert!(e.contains("no material"), "seed {seed0}: {e}"),
        }
    }
}

#[test]
fn pinned_vertices_land_on_the_original_surface() {
    // A ball domain, fully solid: the result should be the ball itself, smoothed onto the true
    // sphere instead of the voxel staircase.
    let (h, r, c) = (1.0, 6.4, [8.0, 8.0, 8.0]);
    let dims = [16, 16, 16];
    let sdf = sphere_sdf(c, r);
    let cells: Vec<f64> = cells_from_sdf(dims, h, &sdf).iter().map(|&v| if v >= 0.5 { 1.0 } else { 0.0 }).collect();
    let (sp, st) = uv_sphere(c, r, 96);
    let proj = Projector::new(&sp, &st).unwrap();
    let f = Field::from_cells(dims, h, [0.0; 3], &cells, 1);
    let opts = Options { blur_passes: 0, target_triangles: Some(2000), ..Default::default() };
    let (m, rep) = extract(f, None, Some(&proj), &opts).unwrap();
    println!("{rep:?}");
    assert!(rep.mesh.ok());
    assert!(rep.pinned_vertices > 100, "pinned {}", rep.pinned_vertices);
    assert!(rep.collapses > 0, "pinned surfaces must still decimate");
    // Every output vertex is on the tessellated sphere (radius error well under a cell).
    let worst = m.positions.iter().map(|&p| (sdf(p)).abs()).fold(0.0, f64::max);
    assert!(worst < 0.05 * h, "worst radial error {worst}");
    let exact = 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
    assert!((rep.mesh.volume - exact).abs() / exact < 0.02, "volume {} vs {exact}", rep.mesh.volume);
    // STL: 84-byte header + 50 bytes per triangle.
    assert_eq!(m.to_stl("test").len(), 84 + 50 * m.triangles.len());
}
