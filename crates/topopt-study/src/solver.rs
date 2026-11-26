//! The file protocol between the plugin and TopOpt Solver, a native multi-threaded program the
//! user runs beside SolidCraft. They share only the plugin's `results` folder:
//!
//! | File (inside `results`) | Written by | What |
//! |---|---|---|
//! | `topopt/solver/alive.json` | solver, every second | [`Alive`]: time, threads, current job, progress |
//! | `topopt/solver/jobs/<id>.partN` | plugin | the [`Job`](crate::Job) as JSON, in parts |
//! | `topopt/solver/jobs/<id>.json` | plugin | [`JobHeader`]: how many parts and bytes |
//! | `<job.mesh_path>` (`topopt/<body>.ply`) | solver | the mesh, as the plugin itself would write it |
//! | `topopt/solver/results/<id>.density.partN` | solver | the densities ([`Density`](crate::Density)) for Re-mesh |
//! | `topopt/solver/results/<id>.json` | solver, last | [`ResultFile`]: the run's summary or its error |
//!
//! Everything the plugin reads or writes is split into parts of at most [`PART_BYTES`], under the
//! host's 1 MiB limit on one plugin file read or write. A plugin run's writes are only published
//! when the run ends, so a job is complete when its header and all its parts exist with the
//! announced total size.

use crate::RunResult;
use serde::{Deserialize, Serialize};

/// One part of a file the plugin reads or writes (below the host's 1 MiB limit).
pub const PART_BYTES: usize = 1_000_000;
/// The solver is taken as running when its heartbeat is at most this old.
pub const ALIVE_SECONDS: u64 = 10;
pub const PROTOCOL: u32 = 1;

pub const ALIVE: &str = "topopt/solver/alive.json";
pub const JOBS: &str = "topopt/solver/jobs";
pub const RESULTS: &str = "topopt/solver/results";

pub fn job_header(id: &str) -> String {
    format!("{JOBS}/{id}.json")
}
pub fn job_part(id: &str, k: usize) -> String {
    format!("{JOBS}/{id}.part{k}")
}
pub fn job_cancel(id: &str) -> String {
    format!("{JOBS}/{id}.cancel")
}
pub fn result_file(id: &str) -> String {
    format!("{RESULTS}/{id}.json")
}
pub fn density_part(id: &str, k: usize) -> String {
    format!("{RESULTS}/{id}.density.part{k}")
}

/// The solver's heartbeat.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Alive {
    pub protocol: u32,
    pub version: String,
    /// Unix time, seconds.
    pub time: u64,
    pub threads: usize,
    /// The job being solved, if any.
    pub job: Option<String>,
    pub progress: f32,
    pub message: String,
}

impl Alive {
    pub fn online(&self, now: u64) -> bool {
        self.protocol == PROTOCOL && now <= self.time + ALIVE_SECONDS
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct JobHeader {
    pub protocol: u32,
    pub id: String,
    pub parts: usize,
    pub bytes: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct ResultFile {
    pub protocol: u32,
    pub id: String,
    pub result: Option<RunResult>,
    pub error: Option<String>,
    pub density_parts: usize,
    pub density_bytes: usize,
    pub threads: usize,
}

/// A job the plugin handed over, kept in the study until its result is taken.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct Pending {
    pub id: String,
    /// Unix time it was sent.
    pub sent: u64,
}

/// What the panel knows about the solver and the pending job (read-only checks).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SolverView {
    /// The heartbeat, when it is fresh.
    pub alive: Option<Alive>,
    /// The pending job's result file, when it has appeared.
    pub ready: Option<ResultFile>,
}

/// Splits bytes into parts the plugin may read or write.
pub fn split(bytes: &[u8]) -> Vec<&[u8]> {
    if bytes.is_empty() { vec![bytes] } else { bytes.chunks(PART_BYTES).collect() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parts_and_heartbeat() {
        let bytes: Vec<u8> = (0..2_500_000u32).map(|i| (i % 251) as u8).collect();
        let parts = split(&bytes);
        assert_eq!(parts.len(), 3);
        assert!(parts.iter().all(|p| p.len() <= PART_BYTES));
        assert_eq!(parts.concat(), bytes);
        let a = Alive { protocol: PROTOCOL, time: 100, ..Default::default() };
        assert!(a.online(105) && !a.online(111));
        assert!(!Alive { time: 100, ..Default::default() }.online(100), "other protocol");
    }
}
