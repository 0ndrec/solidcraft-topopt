//! TopOpt CUDA Solver: TopOpt Solver with its linear systems solved on an NVIDIA GPU
//! (`topopt-cuda`). The same commands, files and folder protocol as `topopt-solver`, plus:
//!
//!   topopt-cuda-solver devices                 list the GPUs the driver sees
//!   --device N, --cpu, --require-gpu, --gpu-min-elements N   (anywhere on the line)
//!
//! Without a usable GPU (no NVIDIA driver, no GPU, too little memory) it solves on the CPU and
//! says so in the run's warnings. See docs/topopt-cuda.md.

use std::process::exit;
use topopt_cuda::{CudaOptions, devices, factory, version_string};
use topopt_solver::cli::{App, main_with};

const NAME: &str = "topopt-cuda-solver";

const GPU_HELP: &str = "\
GPU (anywhere on the line):
  topopt-cuda-solver devices    list the NVIDIA GPUs the driver sees
  --device N             the GPU to use, as numbered by `devices` (default 0)
  --cpu                  solve on the CPU only (as topopt-solver)
  --require-gpu          fail instead of falling back to the CPU when the GPU cannot be used
  --gpu-min-elements N   problems with fewer elements stay on the CPU (default 5000)
";

fn usage(msg: &str) -> ! {
    eprintln!("{NAME}: {msg}\n({NAME} --help for usage)");
    exit(2)
}

fn devices_cmd() -> i32 {
    match devices() {
        Ok((version, list)) => {
            println!("NVIDIA driver, CUDA {}", version_string(version));
            if list.is_empty() {
                println!("no CUDA GPU");
                return 1;
            }
            for d in list {
                println!(
                    "{}: {} — {} MB, compute capability {}.{}, {} multiprocessors",
                    d.index,
                    d.name,
                    d.memory >> 20,
                    d.compute_capability.0,
                    d.compute_capability.1,
                    d.multiprocessors
                );
            }
            0
        }
        Err(e) => {
            println!("{e}");
            1
        }
    }
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut opts = CudaOptions::default();
    let mut cpu = false;
    let mut i = 0;
    while i < args.len() {
        let number = |args: &[String], flag: &str| -> usize {
            let v = args.get(i + 1).unwrap_or_else(|| usage(&format!("{flag} needs a value")));
            v.parse().unwrap_or_else(|_| usage(&format!("{flag} {v}: not a number")))
        };
        match args[i].as_str() {
            "--device" => {
                opts.device = number(&args, "--device");
                args.drain(i..i + 2);
            }
            "--gpu-min-elements" => {
                opts.min_elements = number(&args, "--gpu-min-elements");
                args.drain(i..i + 2);
            }
            "--cpu" => {
                cpu = true;
                args.remove(i);
            }
            "--require-gpu" => {
                opts.fallback = false;
                args.remove(i);
            }
            _ => i += 1,
        }
    }
    if args.first().map(String::as_str) == Some("devices") {
        exit(devices_cmd());
    }
    let runs_on = if cpu {
        "all CPU cores (--cpu)".to_string()
    } else {
        match devices() {
            Ok((version, list)) => match list.get(opts.device) {
                Some(d) => format!("{} (CUDA {}, {} MB) and all CPU cores", d.name, version_string(version), d.memory >> 20),
                None => format!("all CPU cores (no CUDA GPU {})", opts.device),
            },
            Err(e) => format!("all CPU cores ({e})"),
        }
    };
    let app = App {
        name: NAME,
        title: "TopOpt CUDA Solver",
        runs_on,
        solver: if cpu { None } else { Some(factory(opts)) },
        extra_help: GPU_HELP.into(),
    };
    main_with(app, args)
}
