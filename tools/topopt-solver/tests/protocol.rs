//! The solver protocol end to end: a job written in parts the way the plugin writes it, solved by
//! `serve(once)`, read back the way the plugin reads it.

use parry3d_f64::math::Vector;
use parry3d_f64::query::{Ray, RayCast};
use parry3d_f64::shape::TriMesh;
use topopt_study::solver::{self, JobHeader, PROTOCOL, ResultFile};
use topopt_study::{BodyInfo, Caster, Condition, Density, FaceRef, Hit, Hooks, Study, prepare};
use topopt_solver::{Folder, serve};

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
fn job_round_trip_through_the_folder() {
    let root = std::env::temp_dir().join(format!("topopt-solver-test-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    let folder = Folder { root: root.clone() };

    let mut part = BoxPart::new([0.0, 0.0, 0.0], [80.0, 40.0, 20.0]);
    let face = |i| FaceRef { body: "Block".into(), index: i, point: None, label: None };
    let mut study = Study { body: Some("Block".into()), ..Default::default() };
    study.settings.max_elements = 6000;
    study.settings.max_iterations = 30;
    study.add(Condition::Support { id: 0, faces: vec![face(5)], dirs: [true; 3] }, None);
    study.add(Condition::Load { id: 0, faces: vec![face(3)], mode: "vector".into(), force: [0.0, 0.0, -1000.0], magnitude: 0.0, case: 1 }, None);
    let info = BodyInfo { lo: [0.0; 3], hi: [80.0, 40.0, 20.0], volume: 64_000.0 };
    let surface = (part.positions.clone(), part.triangles.clone());
    let mut job = prepare(&study, info, &mut part, &mut [], Some(topopt_study::SurfaceRef::new(&surface.0, &surface.1)), &mut Quiet).unwrap();
    job.id = "j1".into();
    job.mesh_path = "topopt/Block.ply".into();

    // As the plugin writes it: parts first, header last (its writes publish together anyway).
    let bytes = serde_json::to_vec(&job).unwrap();
    let parts = solver::split(&bytes);
    for (k, p) in parts.iter().enumerate() {
        let path = root.join(solver::job_part("j1", k));
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(path, p).unwrap();
    }
    let header = JobHeader { protocol: PROTOCOL, id: "j1".into(), parts: parts.len(), bytes: bytes.len() };
    std::fs::write(root.join(solver::job_header("j1")), serde_json::to_vec(&header).unwrap()).unwrap();
    assert_eq!(folder.ready_jobs(), vec!["j1".to_string()]);

    assert_eq!(serve(&folder, true, &mut |_, _| {}), 1);

    // As the plugin reads it.
    let file: ResultFile = serde_json::from_slice(&std::fs::read(root.join(solver::result_file("j1"))).unwrap()).unwrap();
    assert!(file.error.is_none(), "{:?}", file.error);
    let r = file.result.as_ref().unwrap();
    assert_eq!(r.stl, "topopt/Block.ply");
    assert!(r.triangles > 100 && r.volume_fraction > 0.1, "{r:?}");
    let ply = std::fs::read(root.join("topopt/Block.ply")).unwrap();
    assert!(ply.starts_with(b"ply
"));
    let mut density = Vec::new();
    for k in 0..file.density_parts {
        let part = std::fs::read(root.join(solver::density_part("j1", k))).unwrap();
        assert!(part.len() <= solver::PART_BYTES);
        density.extend(part);
    }
    assert_eq!(density.len(), file.density_bytes);
    let d: Density = serde_json::from_slice(&density).unwrap();
    assert_eq!(d.values.len(), job.problem.grid.nel());
    // The job files are gone, the heartbeat is there.
    assert!(!root.join(solver::job_header("j1")).exists());
    assert!(root.join(solver::ALIVE).exists());
    println!("{r:?}");
    let _ = std::fs::remove_dir_all(&root);
}
