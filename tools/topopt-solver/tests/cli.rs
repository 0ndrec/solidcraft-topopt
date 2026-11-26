//! The command line as a pipeline uses it: a spec with an STL domain next to it, `check`, `solve`
//! with JSON progress, the files it writes, and the exit codes.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use topopt_mesh::{read_ply, read_stl};
use topopt_solver::{Event, Input, RunOptions, run};
use topopt_study::RunResult;
use topopt_study::spec::{Shape, Spec, box_mesh};

fn exe() -> Command {
    Command::new(env!("CARGO_BIN_EXE_topopt-solver"))
}

fn temp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("topopt-cli-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The example, small and fast, with its domain as an STL file beside the spec.
fn write_spec(dir: &Path) -> PathBuf {
    let mut spec = Spec::example();
    spec.name = Some("beam".into());
    spec.settings.max_elements = 4000;
    spec.settings.max_iterations = 15;
    spec.mesh.triangles = 3000;
    std::fs::create_dir_all(dir.join("parts")).unwrap();
    std::fs::write(dir.join("parts/beam.stl"), box_mesh([0.0; 3], [120.0, 40.0, 40.0]).to_stl("beam")).unwrap();
    spec.domain = Shape::Stl { path: "parts/beam.stl".into(), scale: 1.0 };
    let path = dir.join("beam.json");
    std::fs::write(&path, serde_json::to_vec_pretty(&spec).unwrap()).unwrap();
    path
}

fn code(o: &Output) -> i32 {
    o.status.code().unwrap_or(-1)
}

#[test]
fn solve_writes_a_closed_mesh_and_reports_json() {
    let dir = temp("solve");
    let spec = write_spec(&dir);
    let out = dir.join("out");
    let o = exe().arg("solve").arg(&spec).arg("--out").arg(&out).args(["--progress", "json", "--threads", "2"]).output().unwrap();
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert_eq!(code(&o), 0, "{stderr}");

    // stdout: the summary; stderr: one JSON event per line, started first, finished last.
    let result: RunResult = serde_json::from_slice(&o.stdout).unwrap();
    assert!(result.volume_fraction > 0.1 && result.triangles > 100, "{result:?}");
    let events: Vec<serde_json::Value> = stderr.lines().map(|l| serde_json::from_str(l).unwrap()).collect();
    assert_eq!(events.first().unwrap()["event"], "started");
    assert_eq!(events.first().unwrap()["threads"], 2);
    assert_eq!(events.last().unwrap()["event"], "finished");
    assert!(events.iter().filter(|e| e["event"] == "progress").count() > 10);
    let fractions: Vec<f64> = events.iter().filter_map(|e| e["fraction"].as_f64()).collect();
    assert!(fractions.windows(2).all(|w| w[1] >= w[0] - 1e-6), "progress goes forward");

    // The files, named after the spec.
    let mesh = read_ply(&std::fs::read(out.join("beam.ply")).unwrap()).unwrap();
    assert!(mesh.validate(true).ok());
    assert!((mesh.volume() - result.mesh_volume_mm3).abs() < 1e-3 * result.mesh_volume_mm3);
    let saved: RunResult = serde_json::from_slice(&std::fs::read(out.join("beam.result.json")).unwrap()).unwrap();
    assert_eq!(saved, result);
    let density: topopt_study::Density = serde_json::from_slice(&std::fs::read(out.join("beam.density.json")).unwrap()).unwrap();
    assert_eq!(density.values.len(), density.grid.nel());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn check_example_stl_output_and_errors() {
    let dir = temp("misc");
    let spec = write_spec(&dir);

    // check: the problem without solving.
    let o = exe().arg("check").arg(&spec).output().unwrap();
    assert_eq!(code(&o), 0, "{}", String::from_utf8_lossy(&o.stderr));
    let report: serde_json::Value = serde_json::from_slice(&o.stdout).unwrap();
    assert!(report["elements"].as_u64().unwrap() > 1000);
    assert!(report["fixed_nodes"].as_u64().unwrap() > 0 && report["loaded_nodes"].as_u64().unwrap() > 0);

    // example: a spec that parses back.
    let o = exe().arg("example").output().unwrap();
    assert_eq!(code(&o), 0);
    assert_eq!(serde_json::from_slice::<Spec>(&o.stdout).unwrap(), Spec::example());

    // From standard input, mesh only, as STL, with an iteration cap.
    let mut text = std::fs::read_to_string(&spec).unwrap();
    text = text.replace("parts/beam.stl", &dir.join("parts/beam.stl").to_string_lossy().replace('\\', "/"));
    let stl = dir.join("beam.stl");
    let mut child = exe()
        .args(["solve", "-", "--progress", "none", "--max-iterations", "5", "--mesh"])
        .arg(&stl)
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    use std::io::Write;
    child.stdin.take().unwrap().write_all(text.as_bytes()).unwrap();
    let o = child.wait_with_output().unwrap();
    assert_eq!(code(&o), 0, "{}", String::from_utf8_lossy(&o.stderr));
    assert!(o.stderr.is_empty(), "--progress none is silent");
    let result: RunResult = serde_json::from_slice(&o.stdout).unwrap();
    assert!(result.iterations <= 5, "{result:?}");
    assert!(read_stl(&std::fs::read(&stl).unwrap(), 1.0).unwrap().validate(true).ok());
    assert!(!dir.join("beam.density.json").exists(), "only the files asked for");

    // Errors: bad input → 2, a cancel file already there → 3, unknown option → 2.
    std::fs::write(dir.join("bad.json"), b"{\"domain\": 5}").unwrap();
    assert_eq!(code(&exe().arg("solve").arg(dir.join("bad.json")).output().unwrap()), 2);
    assert_eq!(code(&exe().arg("solve").arg(dir.join("missing.json")).output().unwrap()), 2);
    assert_eq!(code(&exe().args(["solve", "--bogus"]).output().unwrap()), 2);
    std::fs::write(dir.join("stop"), b"").unwrap();
    let o = exe().arg("solve").arg(&spec).arg("--no-files").arg("--cancel-file").arg(dir.join("stop")).output().unwrap();
    assert_eq!(code(&o), 3, "{}", String::from_utf8_lossy(&o.stderr));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn library_run_in_memory() {
    // Embedded use: no files, progress to a closure, the mesh in memory.
    let mut spec = Spec::example();
    spec.settings.max_elements = 3000;
    spec.settings.max_iterations = 10;
    let input = Input::Spec { spec, base: PathBuf::new() };
    let mut events = Vec::new();
    let out = run(&input, &RunOptions::default(), &mut |e: &Event| events.push(e.clone())).unwrap();
    assert!(out.files.is_empty());
    assert!(out.mesh.validate(true).ok());
    assert!(matches!(events.first(), Some(Event::Started { .. })));
    assert!(matches!(events.last(), Some(Event::Finished { .. })));
    // A spec without loads fails cleanly, with a Failed event.
    let mut bad = Spec::example();
    bad.loads.clear();
    let mut last = None;
    let err = run(&Input::Spec { spec: bad, base: PathBuf::new() }, &RunOptions::default(), &mut |e: &Event| last = Some(e.clone()));
    assert!(err.is_err());
    assert!(matches!(last, Some(Event::Failed { .. })));
}
