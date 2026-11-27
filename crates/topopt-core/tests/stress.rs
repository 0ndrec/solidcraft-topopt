//! Displacement and stress of a solid cantilever against beam theory.

use topopt_core::{Grid, Optimizer, Problem, Settings};

#[test]
fn a_solid_cantilever_bends_as_beam_theory_says() {
    // 60 × 6 × 6 mm PETG beam, 0.5 mm elements, clamped at x = 0, 10 N down at the free end.
    let (l, b, h, cell) = (60.0, 6.0, 6.0, 0.5);
    let g = Grid::new((l / cell) as usize, (b / cell) as usize, (h / cell) as usize, cell);
    let mut pr = Problem::new(g);
    pr.material.e = 2100.0;
    pr.material.nu = 0.38;
    pr.fix_box([0.0, 0.0, 0.0], [0.0, b, h], [true; 3]);
    let p = 10.0;
    pr.load_box([l, 0.0, 0.0], [l, b, h], [0.0, 0.0, -p]);
    let settings = Settings { volfrac: 0.5, pcg_tol: 1e-10, pcg_max_iter: 2000, ..Default::default() };
    let mut opt = Optimizer::new(pr.clone(), settings).unwrap();
    let check = opt.check(&vec![1.0; g.nel()]).unwrap();
    let c = &check[0];

    // Tip deflection: Euler-Bernoulli plus Timoshenko shear.
    let (e, nu) = (2100.0, 0.38);
    let i = b * h.powi(3) / 12.0;
    let shear = p * l / (5.0 / 6.0 * e / (2.0 * (1.0 + nu)) * b * h);
    let tip = p * l.powi(3) / (3.0 * e * i) + shear;
    // Bending stress at mid-span, at the outermost element centres (z = ±(h/2 − cell/2)).
    let m = p * l / 2.0;
    let sigma = m * (h / 2.0 - cell / 2.0) / i;
    let mid = (0..g.nel())
        .filter(|&el| {
            let [x, _, z] = g.elem_center(el);
            (x - l / 2.0).abs() < cell && (z - h / 2.0).abs() > h / 2.0 - cell
        })
        .map(|el| c.von_mises[el])
        .fold(0.0, f64::max);
    println!("tip {:.4} mm (beam theory {tip:.4}); mid-span outer fibre {mid:.3} MPa (beam theory {sigma:.3})", c.max_displacement);
    // Fully integrated Hex8 is a little stiff in bending; 12 elements through the depth keep it
    // within a few percent.
    assert!((c.max_displacement / tip - 1.0).abs() < 0.05, "tip {} vs {tip}", c.max_displacement);
    assert!((mid / sigma - 1.0).abs() < 0.05, "stress {mid} vs {sigma}");
    assert!((c.compliance - p * c.max_displacement).abs() / c.compliance < 0.05);
}
