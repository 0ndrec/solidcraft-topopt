//! TopOpt Solver — topology optimization on all CPU cores.
//!
//!   topopt-solver solve INPUT.json [options]   one study → mesh, densities, summary
//!   topopt-solver check INPUT.json             build the problem only, print a report
//!   topopt-solver remesh INPUT.json --density D.json --mesh OUT   mesh saved densities again
//!   topopt-solver example                      print an example spec
//!   topopt-solver watch [FOLDER] [--once]      serve the TopOpt plugin (its "results" folder)
//!   topopt-solver version
//!
//! `topopt-solver FOLDER` (and no arguments at all) means `watch`, as before. See
//! docs/topopt-solver.md for the spec format, outputs, progress lines and exit codes.

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::exit;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use topopt_solver::{Event, Folder, Input, RunOptions, VERSION, check, mesh_bytes, remesh, run, serve, write_atomic};

/// Exit codes.
const OK: i32 = 0;
const FAILED: i32 = 1;
const USAGE: i32 = 2;
const CANCELLED: i32 = 3;

const HELP: &str = "\
TopOpt Solver — topology optimization on all CPU cores

usage:
  topopt-solver solve INPUT.json [options]   solve one study (a spec, or a TopOpt plugin job)
  topopt-solver check INPUT.json             build the problem only; print a JSON report
  topopt-solver remesh INPUT.json --density FILE --mesh FILE [--trim PERCENT]
                                             make the mesh again from a run's densities with
                                             INPUT's mesh settings (no optimization); --trim
                                             removes thin dead weight costing at most PERCENT
                                             of the stiffness (overrides mesh.trim_percent)
  topopt-solver example                      print an example spec (JSON) to start from
  topopt-solver watch [FOLDER] [--once]      solve the TopOpt plugin's jobs in its results folder
  topopt-solver version

INPUT may be - to read standard input (relative file paths then resolve from the current folder).

solve options:
  --out DIR              write NAME.ply, NAME.density.json, NAME.result.json into DIR
                         (the default when no output option is given: the current folder)
  --mesh FILE            the mesh (.ply or .stl by extension)
  --density FILE         the densities (JSON)
  --result FILE          the summary (JSON); it is also printed to standard output
  --name NAME            the NAME in --out files (default: the spec's name, else the file name)
  --no-files             write no files; only print the summary
  --progress MODE        json | text | none: progress on standard error (default: text)
  --time-limit S         stop the optimization after S seconds (keeps the best design so far)
  --max-iterations N     at most N optimization iterations
  --cancel-file FILE     stop (exit code 3) as soon as FILE exists
  --threads N            worker threads (default: all cores)

watch options:
  --once                 solve the jobs already there, then exit
  --threads N

exit codes: 0 done, 1 the run failed, 2 bad arguments or input, 3 cancelled
";

fn usage(msg: &str) -> ! {
    eprintln!("topopt-solver: {msg}\n(topopt-solver --help for usage)");
    exit(USAGE)
}

fn set_threads(n: &str) {
    let n: usize = n.parse().ok().filter(|&n| n > 0).unwrap_or_else(|| usage(&format!("--threads {n}: not a positive number")));
    // SAFETY: still single-threaded here; rayon reads this when its pool starts.
    unsafe { std::env::set_var("RAYON_NUM_THREADS", n.to_string()) };
}

fn read_input(arg: &str) -> Input {
    let parsed = if arg == "-" {
        let mut bytes = Vec::new();
        std::io::stdin().read_to_end(&mut bytes).unwrap_or_else(|e| usage(&format!("reading standard input: {e}")));
        Input::from_json(&bytes, std::env::current_dir().unwrap_or_default()).map_err(|e| format!("standard input: {e}"))
    } else {
        Input::from_file(Path::new(arg))
    };
    parsed.unwrap_or_else(|e| usage(&e))
}

fn print_json<T: serde::Serialize>(value: &T) {
    let mut out = std::io::stdout().lock();
    let _ = serde_json::to_writer_pretty(&mut out, value);
    let _ = writeln!(out);
}

#[derive(Clone, Copy, PartialEq)]
enum Progress {
    Json,
    Text,
    None,
}

/// Progress on standard error: every event as one JSON line, or readable lines at most once per
/// percent or message change.
fn reporter(mode: Progress) -> impl FnMut(&Event) {
    let mut last: (i32, String) = (-1, String::new());
    move |e: &Event| {
        let mut err = std::io::stderr().lock();
        match mode {
            Progress::None => {}
            Progress::Json => {
                let _ = serde_json::to_writer(&mut err, e);
                let _ = writeln!(err);
            }
            Progress::Text => match e {
                Event::Started { name, threads, version } => {
                    let _ = writeln!(err, "TopOpt Solver {version}: {name} on {threads} threads");
                }
                Event::Progress { fraction, message, elapsed } => {
                    let pct = (fraction * 100.0).floor() as i32;
                    if pct != last.0 || *message != last.1 {
                        let _ = writeln!(err, "{pct:3}%  {elapsed:7.1} s  {message}");
                        last = (pct, message.clone());
                    }
                }
                Event::Finished { result, files } => {
                    let _ = writeln!(err, "done: {} in {:.1} s, {} iterations", result.outcome, result.seconds, result.iterations);
                    for f in files {
                        let _ = writeln!(err, "  wrote {f}");
                    }
                }
                Event::Failed { message } => {
                    let _ = writeln!(err, "failed: {message}");
                }
            },
        }
    }
}

fn solve_cmd(args: &[String]) -> i32 {
    let mut input = None;
    let (mut out, mut mesh, mut density, mut result, mut name) = (None, None, None, None, None);
    let mut no_files = false;
    let mut progress = Progress::Text;
    let mut opts = RunOptions::default();
    let mut cancel_file: Option<PathBuf> = None;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut value = |flag: &str| it.next().cloned().unwrap_or_else(|| usage(&format!("{flag} needs a value")));
        match a.as_str() {
            "--out" => out = Some(PathBuf::from(value(a))),
            "--mesh" => mesh = Some(PathBuf::from(value(a))),
            "--density" => density = Some(PathBuf::from(value(a))),
            "--result" => result = Some(PathBuf::from(value(a))),
            "--name" => name = Some(value(a)),
            "--no-files" => no_files = true,
            "--progress" => {
                progress = match value(a).as_str() {
                    "json" => Progress::Json,
                    "text" => Progress::Text,
                    "none" => Progress::None,
                    m => usage(&format!("--progress {m}: use json, text or none")),
                }
            }
            "--time-limit" => {
                let v = value(a);
                opts.time_limit = Some(v.parse().ok().filter(|&t: &f64| t > 0.0).unwrap_or_else(|| usage(&format!("--time-limit {v}: not a positive number"))));
            }
            "--max-iterations" => {
                let v = value(a);
                opts.max_iterations = Some(v.parse().ok().filter(|&n| n > 0).unwrap_or_else(|| usage(&format!("--max-iterations {v}: not a positive number"))));
            }
            "--cancel-file" => cancel_file = Some(PathBuf::from(value(a))),
            "--threads" => set_threads(&value(a)),
            f if f.starts_with("--") => usage(&format!("unknown option {f}")),
            path if input.is_none() => input = Some(path.to_string()),
            extra => usage(&format!("unexpected argument {extra}")),
        }
    }
    let input = read_input(&input.unwrap_or_else(|| usage("solve needs an INPUT.json (or - for standard input)")));
    if !no_files {
        if out.is_none() && mesh.is_none() && density.is_none() && result.is_none() {
            out = Some(PathBuf::from("."));
        }
        if let Some(dir) = &out {
            let n = name.unwrap_or_else(|| input.name());
            mesh.get_or_insert_with(|| dir.join(format!("{n}.ply")));
            density.get_or_insert_with(|| dir.join(format!("{n}.density.json")));
            result.get_or_insert_with(|| dir.join(format!("{n}.result.json")));
        }
        (opts.mesh_out, opts.density_out, opts.result_out) = (mesh, density, result);
    }
    if let Some(path) = cancel_file {
        let flag = Arc::new(AtomicBool::new(false));
        opts.cancel = Some(flag.clone());
        std::thread::spawn(move || loop {
            if path.exists() {
                flag.store(true, Ordering::Relaxed);
                return;
            }
            std::thread::sleep(std::time::Duration::from_millis(200));
        });
    }
    match run(&input, &opts, &mut reporter(progress)) {
        Ok(out) => {
            print_json(&out.result);
            OK
        }
        Err(e) if e == "cancelled" => CANCELLED,
        Err(e) => {
            // The progress stream already reported it, unless it is off.
            if progress == Progress::None {
                eprintln!("topopt-solver: {e}");
            }
            FAILED
        }
    }
}

fn remesh_cmd(args: &[String]) -> i32 {
    let (mut input, mut density, mut mesh) = (None, None, None);
    let mut trim: Option<f64> = None;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        let mut value = |flag: &str| it.next().cloned().unwrap_or_else(|| usage(&format!("{flag} needs a value")));
        match a.as_str() {
            "--density" => density = Some(PathBuf::from(value(a))),
            "--mesh" => mesh = Some(PathBuf::from(value(a))),
            "--threads" => set_threads(&value(a)),
            "--trim" => {
                let v = value(a);
                trim = Some(v.parse().ok().filter(|&t: &f64| (0.0..=10.0).contains(&t)).unwrap_or_else(|| usage(&format!("--trim {v}: the stiffness loss allowed, in percent (0…10)"))));
            }
            f if f.starts_with("--") => usage(&format!("unknown option {f}")),
            path if input.is_none() => input = Some(path.to_string()),
            extra => usage(&format!("unexpected argument {extra}")),
        }
    }
    let input = read_input(&input.unwrap_or_else(|| usage("remesh needs an INPUT.json")));
    let density_path = density.unwrap_or_else(|| usage("remesh needs --density FILE"));
    let mesh_path = mesh.unwrap_or_else(|| usage("remesh needs --mesh FILE"));
    let bytes = std::fs::read(&density_path).unwrap_or_else(|e| usage(&format!("{}: {e}", density_path.display())));
    let density: topopt_study::Density =
        serde_json::from_slice(&bytes).unwrap_or_else(|e| usage(&format!("{}: not a density file: {e}", density_path.display())));
    match remesh(&input, &density, trim) {
        Ok((mesh, report)) => {
            if let Err(e) = write_atomic(&mesh_path, &mesh_bytes(&mesh, &mesh_path, &format!("TopOpt {}", input.name()))) {
                eprintln!("topopt-solver: {}: {e}", mesh_path.display());
                return FAILED;
            }
            print_json(&report);
            OK
        }
        Err(e) => {
            eprintln!("topopt-solver: {e}");
            FAILED
        }
    }
}

fn check_cmd(args: &[String]) -> i32 {
    let [path] = args else { usage("check needs exactly one INPUT.json") };
    match check(&read_input(path)) {
        Ok(report) => {
            print_json(&report);
            OK
        }
        Err(e) => {
            eprintln!("topopt-solver: {e}");
            USAGE
        }
    }
}

fn watch_cmd(args: &[String]) -> i32 {
    let mut folder = None;
    let mut once = false;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--once" => once = true,
            "--threads" => set_threads(it.next().map(String::as_str).unwrap_or_else(|| usage("--threads needs a value"))),
            f if f.starts_with("--") => usage(&format!("unknown option {f}")),
            path if folder.is_none() => folder = Some(PathBuf::from(path)),
            extra => usage(&format!("unexpected argument {extra}")),
        }
    }
    let root = folder.unwrap_or_else(|| std::env::current_dir().expect("current folder"));
    if !root.is_dir() {
        usage(&format!("{} is not a folder", root.display()));
    }
    let folder = Folder { root: root.clone() };
    println!("TopOpt Solver {VERSION} — {} threads", topopt_solver::threads());
    println!("Watching {}", root.display());
    println!("This must be the folder the TopOpt plugin uses as \"results\". Close this window to stop.");
    let mut last = std::collections::HashMap::<String, i32>::new();
    let solved = serve(&folder, once, &mut |id, e| match e {
        Event::Started { .. } => println!("job {id}: started"),
        Event::Progress { fraction, message, .. } => {
            // One line per 10 %, so the window stays readable.
            let step = (fraction * 10.0).floor() as i32;
            if last.insert(id.to_string(), step) != Some(step) {
                println!("job {id}: {:3}%  {message}", step * 10);
            }
        }
        Event::Finished { result, .. } => println!("job {id}: done — {} in {:.1} s", result.outcome, result.seconds),
        Event::Failed { message } => println!("job {id}: failed — {message}"),
    });
    if once {
        println!("{solved} job(s) solved");
    }
    OK
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let code = match args.first().map(String::as_str) {
        Some("solve") => solve_cmd(&args[1..]),
        Some("check") => check_cmd(&args[1..]),
        Some("remesh") => remesh_cmd(&args[1..]),
        Some("watch") => watch_cmd(&args[1..]),
        Some("example") => {
            print_json(&topopt_solver::spec::Spec::example());
            OK
        }
        Some("version" | "--version" | "-V") => {
            println!("topopt-solver {VERSION}");
            OK
        }
        Some("help" | "-h" | "--help") => {
            print!("{HELP}");
            OK
        }
        // `topopt-solver FOLDER [--once]` and no arguments: watch, as in 0.5.
        _ => watch_cmd(&args),
    };
    exit(code)
}
