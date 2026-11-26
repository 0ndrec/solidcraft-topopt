//! From an optimized density field to a closed, 2-manifold, outward-oriented triangle mesh that
//! SolidCraft accepts as a closed mesh body and can convert to a solid.
//!
//! [`extract`] runs the whole pipeline:
//! 1. **field** ([`field`]): optional blur, then islands removed (solid components without an
//!    anchor such as a support or load pad) and closed cavities filled, using the exact
//!    connectivity of the contour; values nudged off the iso value;
//! 2. **contour** ([`mt`]): marching tetrahedra, watertight, manifold and self-intersection-free
//!    by construction;
//! 3. **smooth** ([`smooth`]): Taubin λ|μ; vertices near the original part surface are projected
//!    back onto it every step, so mounting faces and holes keep their CAD geometry. Retried with
//!    fewer iterations if it ever makes the surface cross itself;
//! 4. **decimate** ([`decimate`]): QEM with manifold and flip guards, to a triangle budget;
//! 5. **validate** ([`mesh::Mesh::validate`]): every guarantee checked, not assumed.

pub mod decimate;
pub mod field;
mod intersect;
pub mod io;
pub mod mesh;
pub mod mt;
pub mod smooth;

pub use decimate::{DecimateOptions, decimate};
pub use field::{Cleanup, Field, clean_components};
pub use io::{read_ply, read_stl};
pub use mesh::{Mesh, MeshReport};
pub use mt::marching_tetrahedra;
pub use smooth::{Projector, taubin};

#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub iso: f64,
    /// Samples closer than this to `iso` are moved off it (also the smallest edge fraction).
    pub nudge: f64,
    /// Blur passes before contouring. Off by default: blurring moves the surface inwards on thin
    /// convex members (−7 % volume on a torus of minor radius 3 cells, against −1.8 % without).
    pub blur_passes: usize,
    pub remove_islands: bool,
    pub fill_cavities: bool,
    pub smooth_iterations: usize,
    pub lambda: f64,
    pub mu: f64,
    /// Vertices within this many sample spacings of the original surface are kept on it.
    pub pin_band: f64,
    /// Decimate to this many triangles (`None`: keep all).
    pub target_triangles: Option<usize>,
    /// Largest RMS surface change decimation may make, in sample spacings.
    pub max_decimation_error: f64,
    pub check_self_intersections: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            iso: 0.5,
            nudge: 0.02,
            blur_passes: 0,
            remove_islands: true,
            fill_cavities: true,
            smooth_iterations: 10,
            lambda: 0.5,
            mu: -0.53,
            pin_band: 1.0,
            target_triangles: Some(40_000),
            max_decimation_error: 0.25,
            check_self_intersections: true,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Report {
    pub cleanup: Cleanup,
    pub contour_triangles: usize,
    pub contour_volume: f64,
    pub smoothed_volume: f64,
    pub pinned_vertices: usize,
    /// Smoothing iterations actually kept (fewer than asked if smoothing caused crossings).
    pub smooth_iterations: usize,
    pub decimated_from: usize,
    pub collapses: usize,
    pub mesh: MeshReport,
}

/// Runs the pipeline on a field (consumed). `anchors`: samples whose solid component must be kept
/// (one flag per sample). `surface`: the original part, for pinning (see [`Options::pin_band`]).
pub fn extract(mut field: Field, anchors: Option<&[bool]>, surface: Option<&Projector>, opts: &Options) -> Result<(Mesh, Report), String> {
    let mut report = Report::default();
    if opts.blur_passes > 0 {
        field.blur(opts.blur_passes);
    }
    field.zero_border();
    field.nudge(opts.iso, opts.nudge);
    report.cleanup = clean_components(&mut field, opts.iso, anchors, opts.remove_islands, opts.fill_cavities);
    let contour = marching_tetrahedra(&field, opts.iso);
    if contour.triangles.is_empty() {
        return Err("the density field has no material above the iso value".into());
    }
    report.contour_triangles = contour.triangles.len();
    report.contour_volume = contour.volume();
    let h = field.h;

    // Pinned vertices: near the original surface, and facing it. A vertex whose projection lands
    // on an edge or corner of the original surface (its direction disagrees with the contour's
    // normal there) is left free: pinning both sides of an edge folds the triangles across it.
    // Normals from a freely smoothed copy: the raw contour's staircase normals point along the
    // axes, not along the surface they approximate.
    let normals = {
        let mut m = contour.clone();
        taubin(&mut m, 5, opts.lambda, opts.mu, &|_, _| None);
        mesh::vertex_normals(&m)
    };
    let mut pinned: Vec<bool> = match surface {
        Some(s) => contour
            .positions
            .iter()
            .zip(&normals)
            .map(|(&p, n)| {
                let (q, d) = s.project(p);
                if d >= opts.pin_band * h {
                    return false;
                }
                if d < 1e-9 * h {
                    return true;
                }
                let dir = [(p[0] - q[0]) / d, (p[1] - q[1]) / d, (p[2] - q[2]) / d];
                (dir[0] * n[0] + dir[1] * n[1] + dir[2] * n[2]).abs() > 0.6
            })
            .collect(),
        None => vec![false; contour.positions.len()],
    };

    // Smooth. If the result crosses itself or has degenerate triangles, free the vertices of the
    // offending triangles and retry; after a few rounds, halve the iterations instead.
    let mut iterations = opts.smooth_iterations;
    let mut mesh;
    let mut rounds = 0;
    loop {
        mesh = contour.clone();
        let constrain = |v: usize, p: [f64; 3]| -> Option<[f64; 3]> {
            match surface {
                Some(s) if pinned[v] => Some(s.project(p).0),
                _ => None,
            }
        };
        if iterations == 0 {
            for v in 0..mesh.positions.len() {
                if let Some(q) = constrain(v, mesh.positions[v]) {
                    mesh.positions[v] = q;
                }
            }
        } else {
            taubin(&mut mesh, iterations, opts.lambda, opts.mu, &constrain);
        }
        let mut bad: Vec<u32> = mesh::degenerate_triangles(&mesh);
        if opts.check_self_intersections {
            for (a, b) in intersect::find_self_intersections(&mesh) {
                bad.push(a);
                bad.push(b);
            }
        }
        if bad.is_empty() {
            break;
        }
        rounds += 1;
        let mut freed = 0;
        for &t in &bad {
            for &v in &mesh.triangles[t as usize] {
                // Free the vertex and its ring, so the fold has room to relax.
                if pinned[v as usize] {
                    pinned[v as usize] = false;
                    freed += 1;
                }
            }
        }
        if freed == 0 || rounds > 6 {
            if iterations == 0 && !pinned.iter().any(|&p| p) {
                return Err("smoothing could not avoid self-intersections".into());
            }
            if iterations == 0 {
                pinned.iter_mut().for_each(|p| *p = false);
            }
            iterations /= 2;
            rounds = 0;
        }
    }
    report.pinned_vertices = pinned.iter().filter(|&&p| p).count();
    report.smooth_iterations = iterations;
    report.smoothed_volume = mesh.volume();

    // Decimate; back off to a larger budget if the result fails validation.
    if let Some(target) = opts.target_triangles {
        if mesh.triangles.len() > target {
            report.decimated_from = mesh.triangles.len();
            let mut budget = target;
            for _ in 0..4 {
                let mut m = mesh.clone();
                let d = DecimateOptions { target_triangles: budget, max_error: opts.max_decimation_error * h, ..Default::default() };
                let collapses = decimate(&mut m, &pinned, &d);
                if m.validate(opts.check_self_intersections).ok() {
                    report.collapses = collapses;
                    mesh = m;
                    break;
                }
                budget = budget * 3 / 2;
            }
        }
    }
    mesh.compact();
    report.mesh = mesh.validate(opts.check_self_intersections);
    if !report.mesh.ok() {
        return Err(format!("the mesh failed validation: {:?}", report.mesh));
    }
    Ok((mesh, report))
}
