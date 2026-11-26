//! From an optimized density field to a closed, 2-manifold, outward-oriented triangle mesh that
//! SolidCraft accepts as a closed mesh body and can convert to a solid.
//!
//! [`extract`] runs the whole pipeline:
//! 1. **field** ([`field`]): optional blur, then islands removed (solid components without an
//!    anchor such as a support or load pad) and closed cavities filled, using the exact
//!    connectivity of the contour; values nudged off the iso value;
//! 2. **contour** ([`mt`]): marching tetrahedra, watertight, manifold and self-intersection-free
//!    by construction;
//! 3. **pin** vertices near the original part surface to it: corners onto its corners, crease
//!    vertices onto its sharp edges, the rest onto its faces, so mounting faces, edges and holes
//!    keep their CAD geometry;
//! 4. **smooth** ([`smooth`]): Taubin, HC or bilateral, pins re-applied after every pass; then the
//!    free vertices are offset along their normals to restore the volume smoothing took. Retried
//!    with fewer pins or iterations if the surface ever crosses itself;
//! 5. **decimate** ([`decimate`]): QEM with manifold and flip guards, to a triangle budget, keeping
//!    corner and edge pins over face pins;
//! 6. **validate** ([`mesh::Mesh::validate`]): every guarantee checked, not assumed.

pub mod blend;
pub mod decimate;
pub mod field;
mod intersect;
pub mod io;
pub mod mesh;
pub mod mt;
pub mod smooth;

pub use blend::{BlendReport, round_joints, round_union, signed_distance};
pub use decimate::{DecimateOptions, decimate, decimate_ranked};
pub use field::{Cleanup, Field, clean_components, upsample_cells, upsample_mask};
pub use io::{read_ply, read_stl};
pub use mesh::{Mesh, MeshReport};
pub use mt::marching_tetrahedra;
pub use smooth::{Projector, Smoother, bilateral, flip_to_creases, hc_laplacian, insert_crease_vertices, restore_volume, taubin};

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
    pub smoother: Smoother,
    pub smooth_iterations: usize,
    /// Taubin factors (also the staircase passes of [`Smoother::Bilateral`]).
    pub lambda: f64,
    pub mu: f64,
    /// HC factors: pull towards the original (α) and towards the previous positions (β).
    pub hc_alpha: f64,
    pub hc_beta: f64,
    /// Bilateral scales: distance between triangles, in sample spacings, and difference of their
    /// unit normals (0.35: normals more than about 30° apart barely mix).
    pub bilateral_sigma: f64,
    pub bilateral_normal_sigma: f64,
    /// After smoothing, offset the free vertices along their normals so the enclosed volume is
    /// what it was before smoothing (with the pins applied).
    pub preserve_volume: bool,
    /// Vertices within this many sample spacings of the original surface are kept on it.
    pub pin_band: f64,
    /// Pin crease vertices onto the sharp edges, and the nearest vertex onto each corner, of the
    /// original surface (see [`Projector::with_feature_angle`]), instead of onto one of its faces.
    pub sharp_edges: bool,
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
            smoother: Smoother::Taubin,
            smooth_iterations: 10,
            lambda: 0.5,
            mu: -0.53,
            hc_alpha: 0.0,
            hc_beta: 0.6,
            bilateral_sigma: 1.0,
            bilateral_normal_sigma: 0.35,
            preserve_volume: true,
            pin_band: 1.0,
            sharp_edges: true,
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
    /// Volume with the pins applied, before smoothing: the target of volume preservation.
    pub pinned_volume: f64,
    pub smoothed_volume: f64,
    /// Vertices kept on the original surface: all of them, and those on its edges and corners.
    pub pinned_vertices: usize,
    pub edge_vertices: usize,
    pub corner_vertices: usize,
    /// Sharp edges were kept (asked for, a surface given, and the mesh did not need the fallback).
    pub sharp_edges: bool,
    /// Vertices inserted on the sharp edges of the original surface.
    pub crease_vertices: usize,
    /// Edge flips that turned a crease of the original surface into mesh edges.
    pub crease_flips: usize,
    /// Smoothing iterations actually kept (fewer than asked if smoothing caused crossings).
    pub smooth_iterations: usize,
    pub decimated_from: usize,
    pub collapses: usize,
    pub mesh: MeshReport,
}

/// How a vertex is held during smoothing (ordered: a higher kind wins in decimation).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Pin {
    Free = 0,
    Face = 1,
    Edge = 2,
    Corner = 3,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Classifies the contour's vertices against the original surface. `normals`: of a freely
/// smoothed copy (the raw contour's staircase normals point along the axes).
fn classify(contour: &Mesh, normals: &[[f64; 3]], s: &Projector, h: f64, opts: &Options) -> (Vec<Pin>, Vec<[f64; 3]>) {
    let n = contour.positions.len();
    let mut pins = vec![Pin::Free; n];
    let mut corner_at = vec![[0.0; 3]; n];
    for v in 0..n {
        let (p, nv) = (contour.positions[v], normals[v]);
        let (q, d, face_normal) = s.project_with_normal(p);
        if d >= opts.pin_band * h {
            continue;
        }
        // On a sharp edge: the vertices inserted there, and contour vertices that fall on it.
        if opts.sharp_edges && s.nearest_edge(p).is_some_and(|e| e.distance < 1e-6 * h) {
            pins[v] = Pin::Edge;
            continue;
        }
        // A face pin only where the vertex is that face of the material: its outward normal
        // agrees with the face's. The inner side of a thin plate lying along a face points the
        // other way: pinning it too would flatten the plate to nothing, its two sides on top of
        // each other.
        if dot(nv, face_normal) <= 0.6 {
            continue;
        }
        if d < 1e-9 * h {
            pins[v] = Pin::Face;
            continue;
        }
        // And only where the vertex faces the surface: a vertex whose projection lands on an
        // edge or corner of the original surface (its direction disagrees with the contour's
        // normal there) is left free, as pinning both sides of an edge folds the triangles.
        let dir = [(p[0] - q[0]) / d, (p[1] - q[1]) / d, (p[2] - q[2]) / d];
        if dot(dir, nv).abs() > 0.6 {
            pins[v] = Pin::Face;
        }
    }
    if opts.sharp_edges {
        // Each corner takes the one nearest vertex within reach (several on one point would
        // collapse triangles).
        for &c in s.corners() {
            let best = (0..n)
                .filter(|&v| pins[v] != Pin::Corner)
                .map(|v| {
                    let p = contour.positions[v];
                    (v, ((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)).sqrt())
                })
                .filter(|&(_, d)| d < opts.pin_band * h)
                .min_by(|a, b| a.1.total_cmp(&b.1));
            if let Some((v, _)) = best {
                pins[v] = Pin::Corner;
                corner_at[v] = c;
            }
        }
    }
    (pins, corner_at)
}

/// Runs the pipeline on a field (consumed). `anchors`: samples whose solid component must be kept
/// (one flag per sample). `surface`: the original part, for pinning (see [`Options::pin_band`]).
/// If keeping the sharp edges makes the surface cross itself beyond repair, the mesh is made again
/// without them (`Report::sharp_edges` then says false).
pub fn extract(field: Field, anchors: Option<&[bool]>, surface: Option<&Projector>, opts: &Options) -> Result<(Mesh, Report), String> {
    if opts.sharp_edges && surface.is_some() {
        if let Ok(done) = extract_once(field.clone(), anchors, surface, opts) {
            return Ok(done);
        }
        return extract_once(field, anchors, surface, &Options { sharp_edges: false, ..*opts });
    }
    extract_once(field, anchors, surface, opts)
}

fn extract_once(mut field: Field, anchors: Option<&[bool]>, surface: Option<&Projector>, opts: &Options) -> Result<(Mesh, Report), String> {
    let mut report = Report { sharp_edges: opts.sharp_edges && surface.is_some(), ..Default::default() };
    if opts.blur_passes > 0 {
        field.blur(opts.blur_passes);
    }
    field.zero_border();
    field.nudge(opts.iso, opts.nudge);
    report.cleanup = clean_components(&mut field, opts.iso, anchors, opts.remove_islands, opts.fill_cavities);
    let mut contour = marching_tetrahedra(&field, opts.iso);
    if contour.triangles.is_empty() {
        return Err("the density field has no material above the iso value".into());
    }
    report.contour_triangles = contour.triangles.len();
    report.contour_volume = contour.volume();
    let h = field.h;

    // Vertices on the creases of the original surface, where the contour rounds them off.
    if let (Some(s), true) = (surface, opts.sharp_edges) {
        report.crease_vertices = insert_crease_vertices(&mut contour, s, 0.75 * h).len();
    }

    let normals = {
        let mut m = contour.clone();
        taubin(&mut m, 5, opts.lambda, opts.mu, &|_, _| None);
        mesh::vertex_normals(&m)
    };
    let (mut pins, corner_at) = match surface {
        Some(s) => classify(&contour, &normals, s, h, opts),
        None => (vec![Pin::Free; contour.positions.len()], Vec::new()),
    };

    // Smooth. If the result crosses itself or has degenerate triangles, free the vertices of the
    // offending triangles and retry; after a few rounds, halve the iterations instead.
    let mut iterations = opts.smooth_iterations;
    let mut mesh;
    let mut rounds = 0;
    loop {
        mesh = contour.clone();
        let constrain = |v: usize, p: [f64; 3]| -> Option<[f64; 3]> {
            let s = surface?;
            match pins[v] {
                Pin::Free => None,
                Pin::Face => Some(s.project(p).0),
                Pin::Edge => s.nearest_edge(p).map(|e| e.point),
                Pin::Corner => Some(corner_at[v]),
            }
        };
        for v in 0..mesh.positions.len() {
            if let Some(q) = constrain(v, mesh.positions[v]) {
                mesh.positions[v] = q;
            }
        }
        report.pinned_volume = mesh.volume();
        if iterations > 0 {
            match opts.smoother {
                Smoother::Taubin => taubin(&mut mesh, iterations, opts.lambda, opts.mu, &constrain),
                Smoother::Hc => hc_laplacian(&mut mesh, iterations, opts.hc_alpha, opts.hc_beta, &constrain),
                Smoother::Bilateral => {
                    let staircase = iterations.div_ceil(2);
                    taubin(&mut mesh, staircase, opts.lambda, opts.mu, &constrain);
                    bilateral(&mut mesh, iterations - staircase, opts.bilateral_sigma * h, opts.bilateral_normal_sigma, &constrain);
                }
            }
            if opts.preserve_volume {
                // The free vertices carry the correction; their share of the surface limits it,
                // so a surface pinned almost everywhere does not get bumps in its few free spots.
                let free: Vec<bool> = pins.iter().map(|&p| p == Pin::Free).collect();
                let share = free.iter().filter(|&&f| f).count() as f64 / free.len().max(1) as f64;
                restore_volume(&mut mesh, report.pinned_volume, &free, 0.25 * h * share);
            }
        }
        // Make the creases mesh edges: the contour's triangles often straddle them diagonally.
        report.crease_flips = 0;
        if let (Some(s), true) = (surface, opts.sharp_edges) {
            let on_crease: Vec<bool> = pins.iter().map(|&p| p >= Pin::Edge).collect();
            report.crease_flips = flip_to_creases(&mut mesh, &on_crease, s, 0.02 * h);
        }
        let mut bad: Vec<u32> = mesh::degenerate_triangles(&mesh);
        bad.extend(mesh::folded_triangles(&mesh));
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
                // Free the vertex, so the fold has room to relax.
                if pins[v as usize] != Pin::Free {
                    pins[v as usize] = Pin::Free;
                    freed += 1;
                }
            }
        }
        if freed == 0 || rounds > 6 {
            if iterations == 0 && pins.iter().all(|&p| p == Pin::Free) {
                return Err("smoothing could not avoid self-intersections".into());
            }
            if iterations == 0 {
                pins.iter_mut().for_each(|p| *p = Pin::Free);
            }
            iterations /= 2;
            rounds = 0;
        }
    }
    report.pinned_vertices = pins.iter().filter(|&&p| p != Pin::Free).count();
    report.edge_vertices = pins.iter().filter(|&&p| p == Pin::Edge).count();
    report.corner_vertices = pins.iter().filter(|&&p| p == Pin::Corner).count();
    report.smooth_iterations = iterations;
    report.smoothed_volume = mesh.volume();
    let rank: Vec<u8> = pins.iter().map(|&p| p as u8).collect();

    // Decimate. A result that crosses itself is decimated again from the start with the crossing
    // triangles' vertices locked; if that does not settle it, the budget grows.
    if let Some(target) = opts.target_triangles {
        if mesh.triangles.len() > target {
            report.decimated_from = mesh.triangles.len();
            let mut rank = rank;
            let mut budget = target;
            'budgets: for _ in 0..4 {
                for _ in 0..6 {
                    let mut m = mesh.clone();
                    let d = DecimateOptions { target_triangles: budget, max_error: opts.max_decimation_error * h, ..Default::default() };
                    let collapses = decimate_ranked(&mut m, &rank, &d);
                    let vr = m.validate(opts.check_self_intersections);
                    if vr.ok() {
                        report.collapses = collapses;
                        mesh = m;
                        break 'budgets;
                    }
                    let crossings = if opts.check_self_intersections { intersect::find_self_intersections(&m) } else { Vec::new() };
                    if crossings.is_empty() {
                        break;
                    }
                    // Lock the crossing triangles' corners and their neighbours in the mesh before
                    // decimation (the same vertex ids).
                    let mut lock: Vec<u32> = crossings.iter().flat_map(|&(a, b)| [a, b]).flat_map(|t| m.triangles[t as usize]).collect();
                    lock.sort_unstable();
                    lock.dedup();
                    let is_locked: std::collections::HashSet<u32> = lock.iter().copied().collect();
                    for tri in &mesh.triangles {
                        if tri.iter().any(|v| is_locked.contains(v)) {
                            for &v in tri {
                                rank[v as usize] = decimate::LOCKED;
                            }
                        }
                    }
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
