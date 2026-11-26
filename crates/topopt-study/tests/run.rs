//! The study pipeline end to end on a box "part", with parry3d standing in for `view.pick`.

use parry3d_f64::math::Vector;
use parry3d_f64::query::{PointQuery, Ray, RayCast};
use parry3d_f64::shape::TriMesh;
use topopt_study::voxel::{plan_grid, voxelize};
use topopt_study::{BodyInfo, Caster, Condition, FaceRef, Hit, Hooks, Study, run};

/// A box as a triangle mesh whose faces are numbered like a CAD body's:
/// 0 −z, 1 +z, 2 −y, 3 +x, 4 +y, 5 −x (two triangles each).
struct BoxPart {
    mesh: TriMesh,
    positions: Vec<[f64; 3]>,
    triangles: Vec<[u32; 3]>,
}

impl BoxPart {
    fn new(lo: [f64; 3], hi: [f64; 3]) -> Self {
        let c = |i: usize| [if i & 1 == 1 { hi[0] } else { lo[0] }, if i & 2 == 2 { hi[1] } else { lo[1] }, if i & 4 == 4 { hi[2] } else { lo[2] }];
        // Corner numbering as the unit cube 0..7 with bits x, y, z, mapped to a CCW box order.
        let corners = [c(0), c(1), c(3), c(2), c(4), c(5), c(7), c(6)];
        let triangles: Vec<[u32; 3]> =
            vec![[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]];
        let verts = corners.iter().map(|p| Vector::new(p[0], p[1], p[2])).collect();
        BoxPart { mesh: TriMesh::new(verts, triangles.clone()).unwrap(), positions: corners.to_vec(), triangles }
    }
}

impl Caster for BoxPart {
    fn cast(&mut self, rays: &[([f64; 3], [f64; 3])]) -> Result<Vec<Option<Hit>>, String> {
        let (lo, hi) = (self.positions[0], self.positions[6]);
        Ok(rays
            .iter()
            .map(|(o, d)| {
                let ray = Ray::new(Vector::new(o[0], o[1], o[2]), Vector::new(d[0], d[1], d[2]));
                let t = self.mesh.cast_local_ray(&ray, 1e6, false)?;
                let p = [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]];
                // Which box face the point is on (parry reports no triangle index for TriMesh).
                let near = |a: f64, b: f64| (a - b).abs() < 1e-7;
                let (face, normal) = if near(p[2], lo[2]) {
                    (0, [0.0, 0.0, -1.0])
                } else if near(p[2], hi[2]) {
                    (1, [0.0, 0.0, 1.0])
                } else if near(p[1], lo[1]) {
                    (2, [0.0, -1.0, 0.0])
                } else if near(p[0], hi[0]) {
                    (3, [1.0, 0.0, 0.0])
                } else if near(p[1], hi[1]) {
                    (4, [0.0, 1.0, 0.0])
                } else {
                    (5, [-1.0, 0.0, 0.0])
                };
                Some(Hit { point: p, normal, face })
            })
            .collect())
    }
}

struct Quiet;
impl Hooks for Quiet {
    fn progress(&mut self, _: f32, _: &str) {}
    fn cancelled(&mut self) -> bool {
        false
    }
    fn elapsed(&mut self) -> f64 {
        0.0
    }
}

#[test]
fn voxelized_box_has_the_right_volume_and_face_samples() {
    let mut part = BoxPart::new([0.0, 0.0, 0.0], [80.0, 40.0, 20.0]);
    let (grid, warnings) = plan_grid([0.0, 0.0, 0.0], [80.0, 40.0, 20.0], 0.0, 8000);
    println!("grid {}×{}×{} h {:.3} {warnings:?}", grid.nx, grid.ny, grid.nz, grid.h);
    assert!(grid.nel() <= 8000);
    assert!([grid.nx, grid.ny, grid.nz].iter().all(|n| n % 2 == 0));
    let v = voxelize(&mut part, &grid, |_| {}).unwrap();
    let vol = v.inside.iter().filter(|&&i| i).count() as f64 * grid.h.powi(3);
    println!("voxel volume {vol:.0} (exact 64000), {} ray calls, {} hits", v.calls, v.hits.len());
    assert!((vol - 64_000.0).abs() / 64_000.0 < 0.15);
    assert_eq!(v.unbalanced_rays, 0);
    for face in 0..6 {
        assert!(v.hits.iter().any(|h| h.face == face), "face {face} never hit");
    }
}

#[test]
fn bracket_study_runs_end_to_end() {
    let mut part = BoxPart::new([0.0, 0.0, 0.0], [80.0, 40.0, 20.0]);
    let face = |i| FaceRef { body: "Block".into(), index: i, point: None, label: None };
    let mut study = Study { body: Some("Block".into()), ..Default::default() };
    study.settings.max_elements = 6000;
    study.settings.max_iterations = 60;
    study.mesh.triangles = 4000;
    study.add(Condition::Support { id: 0, faces: vec![face(5)], dirs: [true; 3] }, None);
    study.add(Condition::Load { id: 0, faces: vec![face(3)], mode: "vector".into(), force: [0.0, 0.0, -1000.0], magnitude: 0.0, case: 1 }, None);
    let surface = (part.positions.clone(), part.triangles.clone());
    let info = BodyInfo { lo: [0.0; 3], hi: [80.0, 40.0, 20.0], volume: 64_000.0 };
    let out = run(&study, info, &mut part, &mut [], Some((&surface.0, &surface.1)), &mut Quiet).unwrap();
    let r = &out.result;
    println!("{r:?}");
    assert!(out.mesh.validate(true).ok());
    // Coarse grid (h ≈ 2.5 mm on a 20 mm plate): keeping the CAD faces trims the half-cells the
    // voxels added, so the part comes out lighter than the target.
    assert!(r.volume_fraction > 0.18 && r.volume_fraction < 0.33, "volume fraction {}", r.volume_fraction);
    assert!(r.mass_g > 0.0 && r.compliance > 0.0);
    // The clamped face stays a face of the result: some vertices lie on x = 0.
    let minx = out.mesh.positions.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
    let on = out.mesh.positions.iter().filter(|p| p[0].abs() < 1e-6).count();
    let near = out.mesh.positions.iter().filter(|p| p[0] < 1.0).count();
    println!("min x {minx}, on x=0: {on}, x<1: {near}");
    assert!(on > 10);

    // A load "into the face" of magnitude 1000 on the top face equals a −z vector of 1000.
    let mut by_normal = study.clone();
    by_normal.conditions.retain(|c| matches!(c, Condition::Support { .. }));
    by_normal.add(Condition::Load { id: 0, faces: vec![face(1)], mode: "normal".into(), force: [0.0; 3], magnitude: 1000.0, case: 1 }, None);
    by_normal.settings.max_iterations = 1;
    let (grid, _) = plan_grid([0.0; 3], [80.0, 40.0, 20.0], 0.0, 6000);
    let vox = voxelize(&mut part, &grid, |_| {}).unwrap();
    let built = topopt_study::build::build(&by_normal, grid, &vox, &vec![false; grid.nel()]).unwrap();
    let total: [f64; 3] = built.problem.load_cases[0].forces.iter().fold([0.0; 3], |a, f| [a[0] + f.f[0], a[1] + f.f[1], a[2] + f.f[2]]);
    assert!((total[2] + 1000.0).abs() < 1e-6 && total[0].abs() < 1e-9, "{total:?}");
}



/// Area of the mesh's triangles lying in the plane `axis = value`, within `tol`.
fn area_on_plane(mesh: &topopt_study::Mesh, axis: usize, value: f64, tol: f64) -> f64 {
    mesh.triangles
        .iter()
        .filter(|t| t.iter().all(|&i| (mesh.positions[i as usize][axis] - value).abs() < tol))
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
            let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
            0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
        })
        .sum()
}

#[test]
fn a_kept_surface_stays_whole_and_exact() {
    let face = |i| FaceRef { body: "Block".into(), index: i, point: None, label: None };
    let mut study = Study { body: Some("Block".into()), ..Default::default() };
    study.settings.max_elements = 6000;
    study.settings.max_iterations = 40;
    study.mesh.triangles = 6000;
    study.add(Condition::Support { id: 0, faces: vec![face(5)], dirs: [true; 3] }, None);
    study.add(Condition::Load { id: 0, faces: vec![face(3)], mode: "vector".into(), force: [0.0, 0.0, -1000.0], magnitude: 0.0, case: 1 }, None);
    let info = BodyInfo { lo: [0.0; 3], hi: [80.0, 40.0, 20.0], volume: 64_000.0 };
    let solve = |study: &Study| {
        let mut part = BoxPart::new([0.0, 0.0, 0.0], [80.0, 40.0, 20.0]);
        let surface = (part.positions.clone(), part.triangles.clone());
        run(study, info, &mut part, &mut [], Some((&surface.0, &surface.1)), &mut Quiet).unwrap()
    };
    // The +y side (face 4, 80 × 20 mm) of the plain optimization has openings.
    let plain = solve(&study);
    let open = area_on_plane(&plain.mesh, 1, 40.0, 0.01);
    // Kept, even with "keep the body's faces exact" off: the whole face, exactly in its plane.
    let mut kept = study.clone();
    kept.mesh.keep_faces = false;
    let id = kept.add(Condition::Keep { id: 0, faces: vec![face(4)], thickness: 3.0 }, None);
    assert_eq!(kept.conditions.iter().find(|c| c.id() == id).unwrap().summary().0, "Keep surface");
    assert!(kept.pins_faces());
    let out = solve(&kept);
    assert!(out.mesh.validate(true).ok());
    // Exactly in its plane, but for strips along its edges where it meets the other faces.
    let whole = area_on_plane(&out.mesh, 1, 40.0, 0.01);
    let exact = area_on_plane(&out.mesh, 1, 40.0, 1e-9);
    let beyond = out.mesh.positions.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max);
    println!("+y face: {open:.0} mm² plain, {whole:.0} mm² kept, {exact:.0} mm² exactly (of 1600); max y {beyond}; {} triangles (plain {}); {:?}", out.result.triangles, plain.result.triangles, out.result.warnings);
    // How far the face moved: the distance from points of the original face to the mesh.
    let tm = TriMesh::new(
        out.mesh.positions.iter().map(|p| Vector::new(p[0], p[1], p[2])).collect(),
        out.mesh.triangles.clone(),
    )
    .unwrap();
    let mut worst = (0.0f64, [0.0; 2]);
    let mut over = 0;
    let (mut nx, mut total) = (0, 0);
    for i in 0..=160 {
        for k in 0..=40 {
            let (x, z) = (i as f64 * 0.5, k as f64 * 0.5);
            let d = tm.distance_to_local_point(Vector::new(x, 40.0, z), false);
            total += 1;
            if d > 0.01 {
                over += 1;
                if x > 1.0 && x < 79.0 && z > 1.0 && z < 19.0 {
                    nx += 1;
                }
            }
            if d > worst.0 {
                worst = (d, [x, z]);
            }
        }
    }
    println!("face to mesh: worst {:.4} mm at x {} z {}; {over} of {total} points over 0.01 mm, {nx} of them more than 1 mm from an edge", worst.0, worst.1[0], worst.1[1]);
    assert_eq!(nx, 0, "the kept face moved away from its edges");
    assert!(worst.0 < 0.15, "the kept face moved {} mm", worst.0);
    assert!(exact > 0.95 * 1600.0, "exactly in the plane: {exact}");
    assert!(open < 0.9 * 1600.0, "the plain optimization left the face whole ({open}): the test proves nothing");
    assert!(beyond <= 40.0 + 1e-9, "the kept face moved outward to y = {beyond}");
}

