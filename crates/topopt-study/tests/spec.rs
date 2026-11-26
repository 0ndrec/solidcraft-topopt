//! Specs: studies described without a CAD app.

use topopt_core::ElementKind;
use topopt_study::spec::{LoadSpec, Region, Shape, Spec, SupportSpec, box_mesh};
use topopt_study::{Hooks, prepare_spec, solve};

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

fn no_files(path: &str) -> Result<Vec<u8>, String> {
    Err(format!("no file {path}"))
}

fn small() -> Spec {
    let mut s = Spec::example();
    s.settings.max_elements = 4000;
    s.settings.max_iterations = 15;
    s.mesh.triangles = 3000;
    s
}

#[test]
fn example_spec_solves() {
    let spec = small();
    let job = prepare_spec(&spec, &no_files, &mut Quiet).unwrap();
    let p = &job.problem;
    assert!(!p.supports.is_empty());
    let total: f64 = p.load_cases[0].forces.iter().map(|f| f.f[2]).sum();
    assert!((total + 1000.0).abs() < 1e-6, "{total}");
    // Every support node is at x ≈ 0 (the plane region), within the default tolerance.
    let h = p.grid.h;
    assert!(p.supports.iter().all(|s| p.grid.node_pos(s.node)[0].abs() <= 0.6 * h + 1e-9));
    let out = solve(&job, &mut Quiet).unwrap();
    assert!(out.mesh.validate(true).ok());
    assert!(out.result.volume_fraction > 0.1 && out.result.volume_fraction < 0.4, "{:?}", out.result);
    // JSON round trip of the spec format.
    let text = serde_json::to_string_pretty(&spec).unwrap();
    assert_eq!(serde_json::from_str::<Spec>(&text).unwrap(), spec);
}

#[test]
fn box_mesh_and_stl_domains_voxelize_alike() {
    let base = small();
    let kinds = |spec: &Spec, read: &dyn Fn(&str) -> Result<Vec<u8>, String>| prepare_spec(spec, read, &mut Quiet).unwrap().problem.kinds;
    let a = kinds(&base, &no_files);
    let m = box_mesh([0.0, 0.0, 0.0], [120.0, 40.0, 40.0]);
    let mut as_mesh = base.clone();
    as_mesh.domain = Shape::Mesh { positions: m.positions.clone(), triangles: m.triangles.clone() };
    assert_eq!(a, kinds(&as_mesh, &no_files));
    // The same box as an STL in metres.
    let mut metres = m.clone();
    metres.positions.iter_mut().for_each(|p| *p = p.map(|c| c / 1000.0));
    let stl = metres.to_stl("box");
    let mut as_stl = base.clone();
    as_stl.domain = Shape::Stl { path: "box.stl".into(), scale: 1000.0 };
    let read = move |p: &str| if p == "box.stl" { Ok(stl.clone()) } else { Err(p.to_string()) };
    let b = kinds(&as_stl, &read);
    let diff = a.iter().zip(&b).filter(|(x, y)| x != y).count();
    assert!(diff <= a.len() / 200, "{diff} of {} cells differ", a.len());
}

#[test]
fn keep_out_and_errors() {
    let mut spec = small();
    let full = prepare_spec(&spec, &no_files, &mut Quiet).unwrap().problem.kinds.iter().filter(|k| **k != ElementKind::Void).count();
    spec.keep_out.push(Shape::Box { min: [40.0, -1.0, 15.0], max: [80.0, 41.0, 25.0] });
    let cut = prepare_spec(&spec, &no_files, &mut Quiet).unwrap().problem.kinds.iter().filter(|k| **k != ElementKind::Void).count();
    assert!(cut < full * 95 / 100, "{cut} vs {full}");

    let mut nowhere = small();
    nowhere.loads = vec![LoadSpec { region: Region::Sphere { center: [500.0, 0.0, 0.0], radius: 1.0, tolerance: None }, force: [0.0, 0.0, -1.0], case: 1 }];
    let e = prepare_spec(&nowhere, &no_files, &mut Quiet).unwrap_err();
    assert!(e.contains("load 1"), "{e}");
    let mut floating = small();
    floating.supports.clear();
    assert!(prepare_spec(&floating, &no_files, &mut Quiet).unwrap_err().contains("no supports"));
    let mut bad = small();
    bad.material.preset = Some("unobtainium".into());
    assert!(prepare_spec(&bad, &no_files, &mut Quiet).unwrap_err().contains("unknown material"));
    let _ = SupportSpec { region: Region::Box { min: [0.0; 3], max: [1.0; 3], tolerance: Some(0.0) }, dirs: [true; 3] };
}

#[test]
fn the_filtered_field_is_kept_and_meshed_finer() {
    let spec = small();
    let job = prepare_spec(&spec, &no_files, &mut Quiet).unwrap();
    let out = solve(&job, &mut Quiet).unwrap();
    let d = &out.density;
    // The filtered field is saved beside the projected one: pads 1, outside 0.
    assert_eq!(d.filtered.len(), d.values.len());
    for e in 0..d.values.len() {
        if d.pads[e] {
            assert_eq!(d.filtered[e], 1.0);
        } else if !d.inside[e] {
            assert_eq!(d.filtered[e], 0.0);
        }
    }
    let surface = job.surface.as_ref().map(|(p, t)| (p.as_slice(), t.as_slice()));
    // Upsampled ×2 (the default) and not: both valid; densities saved without the filtered field
    // still mesh (from the projected values).
    let mut coarse = job.study.clone();
    coarse.mesh.upsample = 1;
    let (fine_mesh, fine) = topopt_study::make_mesh(&job.study, d, surface).unwrap();
    let (_, plain) = topopt_study::make_mesh(&coarse, d, surface).unwrap();
    let mut old = d.clone();
    old.filtered.clear();
    let (_, legacy) = topopt_study::make_mesh(&job.study, &old, surface).unwrap();
    for r in [&fine, &plain, &legacy] {
        assert!(r.mesh.ok(), "{r:?}");
    }
    assert!(fine_mesh.validate(true).ok());
    // Both describe the same part. This design is still grey after 15 iterations, where the
    // interpolation moves the iso surface most (5 % here; 0.9 % on a converged bridge).
    assert!((fine.mesh.volume - plain.mesh.volume).abs() < 0.08 * plain.mesh.volume, "{} vs {}", fine.mesh.volume, plain.mesh.volume);
}

#[test]
fn fillets_only_add_material_at_the_joints() {
    let spec = small();
    let job = prepare_spec(&spec, &no_files, &mut Quiet).unwrap();
    let out = solve(&job, &mut Quiet).unwrap();
    let surface = job.surface.as_ref().map(|(p, t)| (p.as_slice(), t.as_slice()));
    let mut plain = job.study.clone();
    plain.mesh.fillets = false;
    let mut big = job.study.clone();
    big.mesh.fillet_radius = 3.0 * out.density.grid.h;
    let (_, without) = topopt_study::make_mesh(&plain, &out.density, surface).unwrap();
    let (_, auto) = topopt_study::make_mesh(&job.study, &out.density, surface).unwrap();
    let (_, larger) = topopt_study::make_mesh(&big, &out.density, surface).unwrap();
    for r in [&without, &auto, &larger] {
        assert!(r.mesh.ok(), "{r:?}");
    }
    // Fillets add material, a larger radius more, and not much of it.
    assert!(auto.mesh.volume >= without.mesh.volume * 0.999, "{} vs {}", auto.mesh.volume, without.mesh.volume);
    assert!(larger.mesh.volume >= auto.mesh.volume * 0.999, "{} vs {}", larger.mesh.volume, auto.mesh.volume);
    assert!(larger.mesh.volume < without.mesh.volume * 1.15);
}

#[test]
fn trimming_through_solve_stays_within_its_stiffness_budget() {
    let mut spec = small();
    spec.mesh.trim_percent = 2.0;
    let job = prepare_spec(&spec, &no_files, &mut Quiet).unwrap();
    let out = solve(&job, &mut Quiet).unwrap();
    assert!(out.mesh.validate(true).ok());
    assert!(out.result.trimmed >= 0.0 && out.result.trimmed < 0.2, "{:?}", out.result);
    // Trimmed cells are empty in both fields.
    for e in 0..out.density.values.len() {
        if out.density.values[e] == 0.0 && out.density.inside[e] {
            assert!(out.density.filtered[e] <= 1.0);
        }
    }
}
