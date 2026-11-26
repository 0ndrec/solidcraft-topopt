//! What the panel tells the user: the next step, the body, the resolution the settings give, how
//! long a run will take, and the last result. Plain sentences, no internal terms.

use crate::solver::SolverView;
use crate::study::{Condition, Study, fmt_n};
use crate::{BodyInfo, voxel};

/// Seconds of optimization per element and iteration inside the plugin sandbox (Cranelift),
/// measured: 65 536 elements at about 2.2 s per iteration.
const SECONDS_PER_ELEMENT_ITERATION: f64 = 34e-6;
/// Iterations a study typically needs (β-continuation to 16 at 40 per step, then settling).
const TYPICAL_ITERATIONS: usize = 170;

/// Rough wall time of Optimize, seconds, inside the plugin (`threads` = 0) or in TopOpt Solver.
pub fn estimate_seconds(study: &Study, elements: usize, threads: usize) -> f64 {
    let iterations = study.settings.max_iterations.min(TYPICAL_ITERATIONS) as f64;
    // Native code is about 1.85× the sandbox on one core; measured 2.6× more on 6 cores.
    let speedup = if threads == 0 { 1.0 } else { 1.85 * (1.0 + 0.32 * (threads.min(16) as f64 - 1.0)) };
    (5.0 + elements as f64 * iterations * SECONDS_PER_ELEMENT_ITERATION / speedup).min(study.settings.time_limit + 20.0)
}

pub fn fmt_duration(s: f64) -> String {
    if s < 60.0 {
        format!("{:.0} s", s.max(1.0))
    } else {
        format!("{:.0} min", (s / 60.0).round().max(1.0))
    }
}

/// The one thing to do next.
pub fn next_step(study: &Study) -> String {
    next_step_with(study, &SolverView::default())
}

/// The one thing to do next, knowing what TopOpt Solver is doing.
pub fn next_step_with(study: &Study, solver: &SolverView) -> String {
    if let Some(p) = &study.pending {
        if let Some(r) = &solver.ready {
            return match &r.error {
                Some(e) => format!("The solver failed: {e}"),
                None => "4. Click Create body (the solver has finished)".into(),
            };
        }
        return match &solver.alive {
            Some(a) if a.job.as_deref() == Some(p.id.as_str()) => {
                format!("Solver: {} ({:.0} %) — Refresh to update", a.message, a.progress * 100.0)
            }
            Some(_) => "Solver: the job is queued — Refresh to update".into(),
            None => "Waiting for TopOpt Solver: start it, or Optimize again to run here".into(),
        };
    }
    if let Some(e) = &study.last_error {
        return format!("The last run failed: {e}");
    }
    let supports = study.conditions.iter().any(|c| matches!(c, Condition::Support { .. }));
    let loads = study.conditions.iter().any(|c| matches!(c, Condition::Load { .. }));
    if study.body.is_none() {
        return "1. Pick the design body".into();
    }
    if !supports {
        return "2. Add a support".into();
    }
    if !loads {
        return "2. Add a load".into();
    }
    let foreign = study.problems();
    if !foreign.is_empty() {
        return format!("Fix: {}", foreign.join("; "));
    }
    match &study.result {
        None => "3. Click Optimize".into(),
        Some(r) if !r.inserted => "4. Click Create body".into(),
        Some(_) => "Done (Optimize again after changes)".into(),
    }
}

/// Key/value rows for the status table. `body` is the measured design body, when there is one.
pub fn rows(study: &Study, body: Option<&BodyInfo>) -> Vec<(String, String)> {
    rows_with(study, body, &SolverView::default())
}

/// The status rows, knowing what TopOpt Solver is doing.
pub fn rows_with(study: &Study, body: Option<&BodyInfo>, solver: &SolverView) -> Vec<(String, String)> {
    let mut r: Vec<(String, String)> = vec![("Next".into(), next_step_with(study, solver))];
    let threads = match (&solver.alive, study.settings.use_solver) {
        (Some(a), true) => a.threads.max(1),
        _ => 0,
    };
    r.push((
        "Runs on".into(),
        match (&solver.alive, study.settings.use_solver) {
            (Some(a), true) => format!("TopOpt Solver, {} threads", a.threads),
            (Some(_), false) => "SolidCraft, 1 core (solver off in Advanced)".into(),
            (None, _) => "SolidCraft, 1 core (TopOpt Solver not running)".into(),
        },
    ));
    match (&study.body, body) {
        (Some(name), Some(b)) => {
            let size: Vec<String> = (0..3).map(|d| fmt_n(((b.hi[d] - b.lo[d]) * 10.0).round() / 10.0)).collect();
            r.push(("Design body".into(), format!("{name}, {} mm", size.join("×"))));
            r.push(("Volume".into(), format!("{} cm³, keep {:.0} %", fmt_n((b.volume / 100.0).round() / 10.0), study.settings.volfrac * 100.0)));
            let s = &study.settings;
            let (grid, warnings) = voxel::plan_grid(b.lo, b.hi, s.element_size, s.max_elements);
            let inside_share = if b.volume > 0.0 {
                (b.volume / ((b.hi[0] - b.lo[0]) * (b.hi[1] - b.lo[1]) * (b.hi[2] - b.lo[2]))).clamp(0.05, 1.0)
            } else {
                1.0
            };
            let elements = (grid.nel() as f64 * inside_share) as usize;
            let thinnest = grid.nx.min(grid.ny).min(grid.nz);
            let _ = thinnest;
            r.push(("Resolution".into(), format!("{} cells of {:.2} mm", fmt_n(elements as f64), grid.h)));
            r.push(("Optimize takes".into(), format!("about {}", fmt_duration(estimate_seconds(study, elements, threads)))));
            for w in warnings {
                r.push(("Note".into(), w));
            }
        }
        (Some(name), None) => r.push(("Design body".into(), format!("{name} (not found in the design)"))),
        (None, _) => r.push(("Design body".into(), "not chosen".into())),
    }
    if !study.keep_out.is_empty() {
        r.push(("Keep empty".into(), study.keep_out.join(", ")));
    }
    let supports = study.conditions.iter().filter(|c| matches!(c, Condition::Support { .. })).count();
    let loads: Vec<&Condition> = study.conditions.iter().filter(|c| matches!(c, Condition::Load { .. })).collect();
    let total: f64 = loads.iter().map(|c| c.force_magnitude()).sum();
    let kept = study.conditions.iter().filter(|c| matches!(c, Condition::Keep { .. })).count();
    r.push((
        "Conditions".into(),
        format!(
            "{supports} support{}, {} load{} ({} N){}",
            if supports == 1 { "" } else { "s" },
            loads.len(),
            if loads.len() == 1 { "" } else { "s" },
            fmt_n(total),
            match kept {
                0 => String::new(),
                1 => ", 1 kept surface".into(),
                n => format!(", {n} kept surfaces"),
            }
        ),
    ));
    let m = &study.material;
    r.push(("Material".into(), m.label()));
    if let Some(res) = &study.result {
        r.push((
            "Last result".into(),
            format!("{:.1} g, {:.0} % of the body", res.mass_g, res.volume_fraction * 100.0),
        ));
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::study::FaceRef;

    #[test]
    fn guides_through_the_steps() {
        let mut s = Study::default();
        assert!(next_step(&s).starts_with("1."));
        s.body = Some("B".into());
        assert!(next_step(&s).contains("support"));
        let f = FaceRef { body: "B".into(), index: 1, point: None, label: Some("−X face".into()) };
        s.add(Condition::Support { id: 0, faces: vec![f.clone()], dirs: [true; 3] }, None);
        assert!(next_step(&s).contains("load"));
        s.add(Condition::Load { id: 0, faces: vec![f], mode: "vector".into(), force: [0.0, 0.0, -2000.0], magnitude: 0.0, case: 1 }, None);
        assert!(next_step(&s).starts_with("3."));
        let body = BodyInfo { lo: [0.0; 3], hi: [80.0, 40.0, 20.0], volume: 64_000.0 };
        let rows = rows(&s, Some(&body));
        let text: String = rows.iter().map(|(k, v)| format!("{k}: {v}
")).collect();
        assert!(text.contains("80×40×20 mm"), "{text}");
        assert!(text.contains("(2\u{2009}000 N)"), "{text}");
        assert_eq!(s.conditions[1].summary().1, "2\u{2009}000 N, −Z");
        assert_eq!(s.conditions[0].where_(), "−X face");
    }
}
