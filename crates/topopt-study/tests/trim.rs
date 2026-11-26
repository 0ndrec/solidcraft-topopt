//! Trimming dead weight: a loose spike goes, the load path stays.

use topopt_core::{Grid, Problem, Settings};
use topopt_study::trim::{TrimSettings, trim_dead_ends};

#[test]
fn a_spike_on_a_beam_is_trimmed_and_the_beam_kept() {
    // A 24-cell cantilever beam, 3×3 cells across, clamped at x = 0, loaded at its free end; a
    // spike one cell thick stands on it near the free end, where the bending moment is almost
    // zero, and carries nothing. (On the middle of the beam the same spike stiffens the top flange
    // locally: removing it there costs 0.24 % of the stiffness, and it rightly stays.)
    let g = Grid::new(24, 8, 8, 1.0);
    let mut p = Problem::new(g);
    p.fix_box([0.0, 0.0, 0.0], [0.0, 8.0, 8.0], [true; 3]);
    p.load_box([24.0, 2.0, 0.0], [24.0, 5.0, 3.0], [0.0, 0.0, -1.0]);
    let beam = |i: usize, j: usize, k: usize| (2..5).contains(&j) && k < 3 && i < 24;
    let spike = |i: usize, j: usize, k: usize| i == 22 && j == 3 && (3..8).contains(&k);
    let x: Vec<f64> = (0..g.nel())
        .map(|e| {
            let [i, j, k] = g.elem_ijk(e);
            if beam(i, j, k) || spike(i, j, k) { 1.0 } else { 0.0 }
        })
        .collect();
    let (out, rep) = trim_dead_ends(&p, &Settings::default(), &x, 0.5, &TrimSettings::default()).unwrap();
    println!("{rep:?}");
    for e in 0..g.nel() {
        let [i, j, k] = g.elem_ijk(e);
        if spike(i, j, k) {
            assert_eq!(out[e], 0.0, "the spike cell {i},{j},{k} is trimmed");
        } else {
            assert_eq!(out[e], x[e], "everything else is unchanged ({i},{j},{k})");
        }
    }
    assert_eq!(rep.removed_elements, 5);
    assert!((rep.compliance_after - rep.compliance_before).abs() <= 1e-3 * rep.compliance_before, "{rep:?}");
}
