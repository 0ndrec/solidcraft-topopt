//! The `topopt-cuda-solver` command line: the GPU run gives the CPU run's answer; GPU options and
//! the fallback. GPU-only checks print "skipped" on machines without one.

use serde_json::Value;
use std::process::Command;

const EXE: &str = env!("CARGO_BIN_EXE_topopt-cuda-solver");

fn example() -> String {
    format!("{}/../../examples/cantilever.json", env!("CARGO_MANIFEST_DIR"))
}

fn has_gpu() -> bool {
    let ok = topopt_cuda::devices().is_ok_and(|(_, d)| !d.is_empty());
    if !ok {
        eprintln!("skipped: no CUDA GPU");
    }
    ok
}

fn solve(extra: &[&str]) -> (i32, Value, String) {
    let mut args = vec!["solve", "--no-files", "--progress", "none", "--max-iterations", "6"];
    args.extend_from_slice(extra);
    let ex = example();
    args.push(&ex);
    let out = Command::new(EXE).args(&args).output().unwrap();
    let json = serde_json::from_slice(&out.stdout).unwrap_or(Value::Null);
    (out.status.code().unwrap_or(-1), json, String::from_utf8_lossy(&out.stderr).into_owned())
}

#[test]
fn the_gpu_run_gives_the_cpu_answer() {
    if !has_gpu() {
        return;
    }
    let (code_cpu, cpu, err) = solve(&["--cpu"]);
    assert_eq!(code_cpu, 0, "{err}");
    let (code_gpu, gpu, err) = solve(&["--require-gpu"]);
    assert_eq!(code_gpu, 0, "{err}");
    assert!(cpu["solver"].as_str().unwrap().starts_with("CPU"), "{}", cpu["solver"]);
    assert!(gpu["solver"].as_str().unwrap().contains("CUDA"), "{}", gpu["solver"]);
    let (a, b) = (cpu["compliance"].as_f64().unwrap(), gpu["compliance"].as_f64().unwrap());
    assert!((a - b).abs() / a < 1e-5, "compliance {a} against {b}");
    assert_eq!(cpu["iterations"], gpu["iterations"]);
    assert_eq!(cpu["triangles"], gpu["triangles"]);
}

#[test]
fn a_missing_gpu_falls_back_to_the_cpu_or_fails_when_required() {
    let (code, out, err) = solve(&["--device", "99"]);
    assert_eq!(code, 0, "{err}");
    assert!(out["solver"].as_str().unwrap().starts_with("CPU"));
    assert!(out["warnings"].as_array().unwrap().iter().any(|w| w.as_str().unwrap().contains("GPU could not be used")), "{out}");
    let (code, _, err) = solve(&["--device", "99", "--require-gpu"]);
    assert_eq!(code, 1, "{err}");
    assert!(err.is_empty() || err.contains("GPU"), "{err}");
}

#[test]
fn devices_and_help() {
    let out = Command::new(EXE).arg("devices").output().unwrap();
    let text = String::from_utf8_lossy(&out.stdout);
    if has_gpu() {
        assert!(out.status.success() && text.contains("0: "), "{text}");
    }
    let help = Command::new(EXE).arg("--help").output().unwrap();
    let help = String::from_utf8_lossy(&help.stdout);
    assert!(help.contains("topopt-cuda-solver solve INPUT.json") && help.contains("--require-gpu"), "{help}");
    let bad = Command::new(EXE).args(["--device", "x", "version"]).output().unwrap();
    assert_eq!(bad.status.code(), Some(2));
}
