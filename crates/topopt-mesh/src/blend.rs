//! Rounded unions on cell fields: fillets where the optimized material meets the parts of the
//! design that are kept as they are (pads under supports and loads, keep-in regions).
//!
//! The two bodies are described by signed distance fields `a` and `b` (negative inside), computed
//! exactly from their cell masks ([`signed_distance`]). Their rounded union ([`round_union`]) is
//!
//! ```text
//! u = max(r − a, 0),   v = max(r − b, 0),   union = max(r, min(a, b)) − √(u² + v²)
//! ```
//!
//! which equals the plain union `min(a, b)` except where a point is within `r` of both bodies, and
//! there gives a circular fillet of radius `r` (exactly, for faces meeting at a right angle). The
//! polynomial smooth minimum was tried first: its fillet is only about 0.6 × its parameter, and it
//! adds a quarter of it along the whole contact, so a "radius" setting did not mean millimetres. The field is only ever raised, and only in that fillet zone:
//! elsewhere the optimized surface stays exactly as it was (the distance fields are exact only to
//! half a cell, while the density field places the surface between samples). The two bodies must be
//! different parts: the rounded union of a body with itself would grow its whole surface by `k/4`.

/// One-dimensional squared Euclidean distance transform of `f` (Felzenszwalb & Huttenlocher
/// 2012): `d[q] = min_p (q − p)² + f[p]`, in O(n).
fn edt_1d(f: &[f64], d: &mut [f64], v: &mut [usize], z: &mut [f64]) {
    let n = f.len();
    let mut k = 0usize;
    v[0] = 0;
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;
    for q in 1..n {
        if !f[q].is_finite() {
            continue;
        }
        loop {
            let p = v[k];
            if !f[p].is_finite() {
                // An empty start: replace it.
                v[k] = q;
                z[k] = f64::NEG_INFINITY;
                z[k + 1] = f64::INFINITY;
                break;
            }
            let s = ((f[q] + (q * q) as f64) - (f[p] + (p * p) as f64)) / (2.0 * (q as f64 - p as f64));
            if s <= z[k] && k > 0 {
                k -= 1;
                continue;
            }
            if s <= z[k] {
                v[k] = q;
                z[k] = f64::NEG_INFINITY;
                z[k + 1] = f64::INFINITY;
                break;
            }
            k += 1;
            v[k] = q;
            z[k] = s;
            z[k + 1] = f64::INFINITY;
            break;
        }
    }
    k = 0;
    for q in 0..n {
        while z[k + 1] < q as f64 {
            k += 1;
        }
        let p = v[k];
        d[q] = if f[p].is_finite() { (q as f64 - p as f64).powi(2) + f[p] } else { f64::INFINITY };
    }
}

/// Squared distance (in cells) from every cell to the nearest cell with `seed` set.
fn squared_distance(dims: [usize; 3], seed: &[bool]) -> Vec<f64> {
    let [nx, ny, nz] = dims;
    let mut g: Vec<f64> = seed.iter().map(|&s| if s { 0.0 } else { f64::INFINITY }).collect();
    let n = nx.max(ny).max(nz);
    let (mut line, mut out, mut v, mut z) = (vec![0.0; n], vec![0.0; n], vec![0usize; n], vec![0.0; n + 1]);
    let idx = |i: usize, j: usize, k: usize| i + nx * (j + ny * k);
    // Along x, then y, then z.
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                line[i] = g[idx(i, j, k)];
            }
            edt_1d(&line[..nx], &mut out[..nx], &mut v, &mut z);
            for i in 0..nx {
                g[idx(i, j, k)] = out[i];
            }
        }
    }
    for k in 0..nz {
        for i in 0..nx {
            for j in 0..ny {
                line[j] = g[idx(i, j, k)];
            }
            edt_1d(&line[..ny], &mut out[..ny], &mut v, &mut z);
            for j in 0..ny {
                g[idx(i, j, k)] = out[j];
            }
        }
    }
    for j in 0..ny {
        for i in 0..nx {
            for k in 0..nz {
                line[k] = g[idx(i, j, k)];
            }
            edt_1d(&line[..nz], &mut out[..nz], &mut v, &mut z);
            for k in 0..nz {
                g[idx(i, j, k)] = out[k];
            }
        }
    }
    g
}

/// Signed distance (in units of the cell size `h`; negative inside) from the cell centres to the
/// surface of the body made of the cells in `inside`, which is taken half a cell beyond its cells'
/// centres. Exact to about half a cell.
pub fn signed_distance(dims: [usize; 3], h: f64, inside: &[bool]) -> Vec<f64> {
    let outside: Vec<bool> = inside.iter().map(|&i| !i).collect();
    let to_in = squared_distance(dims, inside);
    let to_out = squared_distance(dims, &outside);
    (0..inside.len())
        .map(|c| {
            if inside[c] {
                -(to_out[c].sqrt() - 0.5) * h
            } else {
                (to_in[c].sqrt() - 0.5) * h
            }
        })
        .collect()
}

/// Rounded union of two signed distances with fillet radius `r` (see the module docs).
pub fn round_union(a: f64, b: f64, r: f64) -> f64 {
    if r <= 0.0 {
        return a.min(b);
    }
    let (u, v) = ((r - a).max(0.0), (r - b).max(0.0));
    r.max(a.min(b)) - (u * u + v * v).sqrt()
}

/// What [`round_joints`] changed.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct BlendReport {
    /// Cells whose value was raised (the fillets).
    pub raised_cells: usize,
}

/// Adds fillets of radius `radius` where body `a` (the optimized material) meets body `b` (the
/// kept parts), on a cell field `values` whose surface is the `iso` level. Cells outside `domain`
/// are never raised. `h` is the cell size. Only cells in the fillet zone change, and only upwards.
pub fn round_joints(dims: [usize; 3], h: f64, values: &mut [f64], iso: f64, a: &[bool], b: &[bool], domain: &[bool], radius: f64) -> BlendReport {
    let mut report = BlendReport::default();
    if radius <= 0.0 || !a.iter().any(|&x| x) || !b.iter().any(|&x| x) {
        return report;
    }
    let (da, db) = (signed_distance(dims, h, a), signed_distance(dims, h, b));
    // The blended field: the fillet ramp in the fillet zone, the field itself elsewhere.
    let mut target = values.to_vec();
    let mut zone = vec![false; values.len()];
    for c in 0..values.len() {
        if !domain[c] {
            continue;
        }
        let (x, y) = (da[c], db[c]);
        // Only within the fillet radius of both bodies.
        if x >= radius || y >= radius {
            continue;
        }
        let s = round_union(x, y, radius);
        if s >= x.min(y) {
            continue;
        }
        // A density ramp across two cells, `iso` on the blended surface.
        target[c] = values[c].max((iso - s / (2.0 * h)).clamp(0.0, 1.0));
        zone[c] = true;
    }
    // The distances are exact to half a cell only: smooth the blended field so the fillet carries
    // no steps (a [¼ ½ ¼] pass along each axis), and take it in the zone and the cells next to it.
    blur3(dims, &mut target);
    let near = grow(dims, &zone);
    for c in 0..values.len() {
        if near[c] && domain[c] && target[c] > values[c] + 1e-9 {
            values[c] = target[c].min(1.0);
            report.raised_cells += 1;
        }
    }
    report
}

/// The mask grown by one cell along each axis.
fn grow(dims: [usize; 3], m: &[bool]) -> Vec<bool> {
    let strides = [1, dims[0], dims[0] * dims[1]];
    let mut out = m.to_vec();
    for c in 0..m.len() {
        if !m[c] {
            continue;
        }
        for axis in 0..3 {
            let (n, st) = (dims[axis], strides[axis]);
            let i = (c / st) % n;
            if i > 0 {
                out[c - st] = true;
            }
            if i + 1 < n {
                out[c + st] = true;
            }
        }
    }
    out
}

/// One separable [¼ ½ ¼] pass along each axis (edges clamped).
fn blur3(dims: [usize; 3], v: &mut [f64]) {
    let strides = [1, dims[0], dims[0] * dims[1]];
    let mut tmp = v.to_vec();
    for axis in 0..3 {
        let (n, st) = (dims[axis], strides[axis]);
        for c in 0..v.len() {
            let i = (c / st) % n;
            let lo = if i > 0 { v[c - st] } else { v[c] };
            let hi = if i + 1 < n { v[c + st] } else { v[c] };
            tmp[c] = 0.25 * lo + 0.5 * v[c] + 0.25 * hi;
        }
        v.copy_from_slice(&tmp);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distances_are_exact_on_a_single_cell() {
        let dims = [7, 7, 7];
        let mut inside = vec![false; 343];
        inside[3 + 7 * (3 + 7 * 3)] = true;
        let d = signed_distance(dims, 2.0, &inside);
        // The cell (3,3,3) is a cube of side 2: centre −1·h/2·… = −(1 − 0.5)·2.
        assert!((d[3 + 7 * (3 + 7 * 3)] + 1.0).abs() < 1e-12);
        // Cell (6,3,3): 3 cells away, surface half a cell nearer: (3 − 0.5)·2 = 5.
        assert!((d[6 + 7 * (3 + 7 * 3)] - 5.0).abs() < 1e-12);
        // Cell (5,5,3): √8 cells.
        assert!((d[5 + 7 * (5 + 7 * 3)] - (8f64.sqrt() - 0.5) * 2.0).abs() < 1e-12);
    }

    #[test]
    fn round_union_makes_a_circular_fillet_of_radius_r() {
        // Two half-spaces x < 0 and y < 0 meeting at a right angle: the fillet is the circle of
        // radius r centred at (r, r), so on the bisector the surface is at t = r (1 − 1/√2).
        let r = 3.0;
        let t = r * (1.0 - 1.0 / 2f64.sqrt());
        assert!(round_union(t, t, r).abs() < 1e-12);
        // Away from the corner it is the plain union.
        assert_eq!(round_union(5.0, -1.0, r), -1.0);
        assert_eq!(round_union(4.0, 3.5, r), 3.5);
        assert_eq!(round_union(1.0, 2.0, 0.0), 1.0);
    }

    #[test]
    fn distance_transform_matches_brute_force() {
        let dims = [9, 8, 7];
        let n = 9 * 8 * 7;
        let mut seed = 12345u64;
        for trial in 0..20 {
            let mask: Vec<bool> = (0..n)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    (seed >> 33) % (3 + trial as u64) == 0
                })
                .collect();
            if !mask.iter().any(|&m| m) {
                continue;
            }
            let d = squared_distance(dims, &mask);
            let pos = |c: usize| [(c % 9) as f64, ((c / 9) % 8) as f64, (c / 72) as f64];
            for c in 0..n {
                let p = pos(c);
                let best = (0..n)
                    .filter(|&s| mask[s])
                    .map(|s| {
                        let q = pos(s);
                        (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)
                    })
                    .fold(f64::INFINITY, f64::min);
                assert_eq!(d[c], best, "trial {trial}, cell {c}");
            }
        }
    }

    #[test]
    fn a_bar_on_a_plate_gets_a_fillet_and_nothing_else_changes() {
        // A plate z < 5 (kept part) and a bar 4×4 cells standing on it (optimized material).
        let dims = [20, 20, 20];
        let n = 8000;
        let idx = |i: usize, j: usize, k: usize| i + 20 * (j + 20 * k);
        let mut plate = vec![false; n];
        let mut bar = vec![false; n];
        for k in 0..20 {
            for j in 0..20 {
                for i in 0..20 {
                    if k < 5 {
                        plate[idx(i, j, k)] = true;
                    } else if (8..12).contains(&i) && (8..12).contains(&j) {
                        bar[idx(i, j, k)] = true;
                    }
                }
            }
        }
        let mut values: Vec<f64> = (0..n).map(|c| if plate[c] || bar[c] { 1.0 } else { 0.0 }).collect();
        let before = values.clone();
        let rep = round_joints(dims, 1.0, &mut values, 0.5, &bar, &plate, &vec![true; n], 3.0);
        assert!(rep.raised_cells > 0);
        // In the corner next to the bar's foot: filled.
        assert!(values[idx(12, 10, 5)] >= 0.5, "{}", values[idx(12, 10, 5)]);
        // Far along the plate, and high up the bar: untouched.
        assert_eq!(values[idx(18, 10, 5)], before[idx(18, 10, 5)]);
        assert_eq!(values[idx(13, 10, 15)], before[idx(13, 10, 15)]);
        // Only raised, never lowered.
        assert!(values.iter().zip(&before).all(|(v, b)| v >= b));
    }
}
