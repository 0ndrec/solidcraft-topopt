//! Voxelization of a CAD body by ray casting, independent of how rays are cast.
//!
//! For each axis, one ray runs through every row of cell centres. A ray reports only the first
//! surface it meets, so each row is walked hit by hit: after a hit the ray starts again just past
//! it. Entering hits (normal against the ray) open an "inside" interval, exiting hits close it.
//! A cell is inside when at least two of the three axes say so, which forgives a ray grazing an
//! edge. Every hit is kept, with its face: those samples place supports, loads and pads.

use topopt_core::Grid;

/// One ray/surface intersection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    pub point: [f64; 3],
    /// Outward unit normal.
    pub normal: [f64; 3],
    /// Face index, as `model.faces` and selection references number them.
    pub face: i64,
}

/// Casts rays at one body: for each `(origin, unit direction)`, the nearest hit on that body
/// ahead of the origin, if any.
pub trait Caster {
    fn cast(&mut self, rays: &[([f64; 3], [f64; 3])]) -> Result<Vec<Option<Hit>>, String>;
}

pub struct Voxels {
    pub inside: Vec<bool>,
    pub hits: Vec<Hit>,
    /// Rays that found an odd number of crossings (open or self-touching geometry).
    pub unbalanced_rays: usize,
    pub calls: usize,
}

/// The grid over a box: element size `h` (or from `max_elements`), every dimension a multiple
/// of 4 (of 2 when small) so multigrid can coarsen, centred on the box.
pub fn plan_grid(lo: [f64; 3], hi: [f64; 3], element_size: f64, max_elements: usize) -> (Grid, Vec<String>) {
    let mut warnings = Vec::new();
    let size = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
    let round = |n: usize| if n >= 8 { n.div_ceil(4) * 4 } else { n.max(2).div_ceil(2) * 2 };
    let dims_for = |h: f64| -> [usize; 3] { std::array::from_fn(|d| round((size[d] / h).ceil().max(1.0) as usize)) };
    let count = |n: [usize; 3]| n[0] * n[1] * n[2];
    let mut h = if element_size > 0.0 { element_size } else { (size[0] * size[1] * size[2] / max_elements as f64).cbrt() };
    let mut n = dims_for(h);
    // Rounding up to multiples of 4 adds cells; grow h until the count fits.
    while count(n) > max_elements.max(8) {
        if element_size > 0.0 && h == element_size {
            warnings.push(format!(
                "element size {element_size} mm needs {} elements; raised to stay under {max_elements}",
                count(n)
            ));
        }
        h *= 1.03;
        n = dims_for(h);
    }
    let thinnest = n.iter().copied().min().unwrap();
    if thinnest < 8 {
        warnings.push(format!("the body is only {thinnest} elements thick somewhere: use a smaller element size for thin parts"));
    }
    let origin: [f64; 3] = std::array::from_fn(|d| lo[d] - (n[d] as f64 * h - size[d]) / 2.0);
    let mut g = Grid::new(n[0], n[1], n[2], h);
    g.origin = origin;
    (g, warnings)
}

/// Inside/outside of every cell of `grid` for the body the caster sees.
pub fn voxelize(caster: &mut dyn Caster, grid: &Grid, mut progress: impl FnMut(f32)) -> Result<Voxels, String> {
    let n = [grid.nx, grid.ny, grid.nz];
    let h = grid.h;
    let mut votes = vec![0u8; grid.nel()];
    let mut hits = Vec::new();
    let mut unbalanced = 0;
    let mut calls = 0;
    for axis in 0..3 {
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        let mut dir = [0.0; 3];
        dir[axis] = 1.0;
        let start = grid.origin[axis] - h;
        let end = grid.origin[axis] + (n[axis] as f64 + 1.0) * h;
        // One ray per row of cell centres.
        let rows: Vec<[f64; 3]> = (0..n[u] * n[v])
            .map(|r| {
                let (a, b) = (r % n[u], r / n[u]);
                let mut o = [0.0; 3];
                o[u] = grid.origin[u] + (a as f64 + 0.5) * h;
                o[v] = grid.origin[v] + (b as f64 + 0.5) * h;
                o[axis] = start;
                o
            })
            .collect();
        // Crossings per row: (position along the axis, entering?).
        let mut crossings: Vec<Vec<(f64, bool)>> = vec![Vec::new(); rows.len()];
        let mut active: Vec<usize> = (0..rows.len()).collect();
        let mut origin: Vec<[f64; 3]> = rows.clone();
        let eps = 1e-4 * h;
        for _round in 0..256 {
            if active.is_empty() {
                break;
            }
            let mut next = Vec::new();
            for chunk in active.chunks(1024) {
                let rays: Vec<([f64; 3], [f64; 3])> = chunk.iter().map(|&r| (origin[r], dir)).collect();
                let res = caster.cast(&rays)?;
                calls += 1;
                for (&r, hit) in chunk.iter().zip(res) {
                    let Some(hit) = hit else { continue };
                    let s = hit.point[axis];
                    if s > end || s < origin[r][axis] - eps {
                        continue;
                    }
                    let d = hit.normal[axis];
                    if d.abs() > 1e-9 {
                        crossings[r].push((s, d < 0.0));
                    }
                    hits.push(hit);
                    origin[r][axis] = s + eps;
                    next.push(r);
                }
            }
            active = next;
        }
        for (r, list) in crossings.iter().enumerate() {
            // Inside intervals from enter/exit hits.
            let mut intervals = Vec::new();
            let mut open: Option<f64> = None;
            for &(s, entering) in list {
                match (entering, open) {
                    (true, None) => open = Some(s),
                    (false, Some(a)) => {
                        intervals.push((a, s));
                        open = None;
                    }
                    _ => {}
                }
            }
            if open.is_some() || list.len() % 2 == 1 {
                unbalanced += 1;
            }
            let (a, b) = (r % n[u], r / n[u]);
            for c in 0..n[axis] {
                let x = grid.origin[axis] + (c as f64 + 0.5) * h;
                if intervals.iter().any(|&(lo, hi)| x > lo && x < hi) {
                    let mut ijk = [0usize; 3];
                    ijk[axis] = c;
                    ijk[u] = a;
                    ijk[v] = b;
                    votes[grid.elem(ijk[0], ijk[1], ijk[2])] += 1;
                }
            }
        }
        progress((axis + 1) as f32 / 3.0);
    }
    Ok(Voxels { inside: votes.iter().map(|&v| v >= 2).collect(), hits, unbalanced_rays: unbalanced, calls })
}
