//! A scalar field sampled on a regular grid, and the clean-up done before contouring.
//!
//! Element densities are cell-centred: sample `(i, j, k)` of a field built by
//! [`Field::from_cells`] sits at the centre of element `(i − pad, j − pad, k − pad)`, with `pad`
//! rings of zero samples around the grid so the surface always closes.
//!
//! Connectivity matches the contouring exactly. Marching tetrahedra on the Freudenthal
//! decomposition (see [`crate::mt`]) joins two samples on the same side of the iso value only
//! along tetrahedron edges, which are the 14 [`NEIGHBOURS`] below. Components counted with them are
//! exactly the components of the solid and of the void that the surface bounds, for solid and
//! void alike.

/// Tetrahedron-edge neighbours of a sample in the Freudenthal decomposition used by [`crate::mt`].
pub const NEIGHBOURS: [[i32; 3]; 14] = [
    [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
    [1, 1, 0], [-1, -1, 0], [1, 0, 1], [-1, 0, -1], [0, 1, 1], [0, -1, -1],
    [1, 1, 1], [-1, -1, -1],
];

#[derive(Clone, Debug)]
pub struct Field {
    /// Samples per axis.
    pub n: [usize; 3],
    /// Position of sample (0, 0, 0).
    pub origin: [f64; 3],
    /// Sample spacing.
    pub h: f64,
    pub values: Vec<f64>,
}

/// What [`clean_components`] changed.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Cleanup {
    pub solid_components: usize,
    pub removed_components: usize,
    pub removed_samples: usize,
    pub cavities_filled: usize,
    pub filled_samples: usize,
}

impl Field {
    #[inline]
    pub fn idx(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.n[0] * (j + self.n[1] * k)
    }

    #[inline]
    pub fn ijk(&self, s: usize) -> [usize; 3] {
        [s % self.n[0], (s / self.n[0]) % self.n[1], s / (self.n[0] * self.n[1])]
    }

    pub fn pos(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + i as f64 * self.h,
            self.origin[1] + j as f64 * self.h,
            self.origin[2] + k as f64 * self.h,
        ]
    }

    /// Samples at the centres of an `nx × ny × nz` grid of cells of size `h` whose corner (0,0,0)
    /// is at `grid_origin`, with `pad ≥ 1` rings of zeros around it.
    pub fn from_cells(dims: [usize; 3], h: f64, grid_origin: [f64; 3], cells: &[f64], pad: usize) -> Field {
        assert!(pad >= 1, "the field needs at least one ring of zeros to close the surface");
        assert_eq!(cells.len(), dims[0] * dims[1] * dims[2]);
        let n = [dims[0] + 2 * pad, dims[1] + 2 * pad, dims[2] + 2 * pad];
        let origin = [
            grid_origin[0] + (0.5 - pad as f64) * h,
            grid_origin[1] + (0.5 - pad as f64) * h,
            grid_origin[2] + (0.5 - pad as f64) * h,
        ];
        let mut f = Field { n, origin, h, values: vec![0.0; n[0] * n[1] * n[2]] };
        for k in 0..dims[2] {
            for j in 0..dims[1] {
                for i in 0..dims[0] {
                    let s = f.idx(i + pad, j + pad, k + pad);
                    f.values[s] = cells[i + dims[0] * (j + dims[1] * k)];
                }
            }
        }
        f
    }

    /// A per-cell flag moved onto the padded samples of a field made by [`Field::from_cells`].
    pub fn cell_mask(&self, dims: [usize; 3], pad: usize, cells: &[bool]) -> Vec<bool> {
        let mut m = vec![false; self.values.len()];
        for k in 0..dims[2] {
            for j in 0..dims[1] {
                for i in 0..dims[0] {
                    m[self.idx(i + pad, j + pad, k + pad)] = cells[i + dims[0] * (j + dims[1] * k)];
                }
            }
        }
        m
    }

    /// Separable [¼ ½ ¼] smoothing, `passes` times. The outermost ring is reset to zero.
    pub fn blur(&mut self, passes: usize) {
        let n = self.n;
        let mut tmp = vec![0.0; self.values.len()];
        for _ in 0..passes {
            for axis in 0..3 {
                let stride = [1, n[0], n[0] * n[1]][axis];
                for s in 0..self.values.len() {
                    let c = self.ijk(s)[axis];
                    let prev = if c > 0 { self.values[s - stride] } else { 0.0 };
                    let next = if c + 1 < n[axis] { self.values[s + stride] } else { 0.0 };
                    tmp[s] = 0.25 * prev + 0.5 * self.values[s] + 0.25 * next;
                }
                std::mem::swap(&mut self.values, &mut tmp);
            }
        }
        self.zero_border();
    }

    /// Sets the outermost ring of samples to 0 so the surface closes inside the field.
    pub fn zero_border(&mut self) {
        for s in 0..self.values.len() {
            let [i, j, k] = self.ijk(s);
            if i == 0 || j == 0 || k == 0 || i + 1 == self.n[0] || j + 1 == self.n[1] || k + 1 == self.n[2] {
                self.values[s] = 0.0;
            }
        }
    }

    /// Moves samples within `eps` of `iso` away from it. With no sample exactly at the iso value
    /// the contour is a 2-manifold and every edge crossing has `t ∈ [eps, 1 − eps]`.
    pub fn nudge(&mut self, iso: f64, eps: f64) {
        for v in &mut self.values {
            let d = *v - iso;
            if d.abs() < eps {
                *v = if d >= 0.0 { iso + eps } else { iso - eps };
            }
        }
    }

    /// Labels the components of `{v > iso}` (`solid = true`) or `{v ≤ iso}`; returns per-sample
    /// labels (`u32::MAX` for samples on the other side) and the component sizes.
    pub fn components(&self, iso: f64, solid: bool) -> (Vec<u32>, Vec<usize>) {
        let inside = |v: f64| (v > iso) == solid;
        let mut label = vec![u32::MAX; self.values.len()];
        let mut sizes = Vec::new();
        let mut stack = Vec::new();
        for start in 0..self.values.len() {
            if label[start] != u32::MAX || !inside(self.values[start]) {
                continue;
            }
            let id = sizes.len() as u32;
            let mut size = 0;
            label[start] = id;
            stack.push(start);
            while let Some(s) = stack.pop() {
                size += 1;
                let [i, j, k] = self.ijk(s);
                for d in NEIGHBOURS {
                    let (a, b, c) = (i as i32 + d[0], j as i32 + d[1], k as i32 + d[2]);
                    if a < 0 || b < 0 || c < 0 || a >= self.n[0] as i32 || b >= self.n[1] as i32 || c >= self.n[2] as i32 {
                        continue;
                    }
                    let t = self.idx(a as usize, b as usize, c as usize);
                    if label[t] == u32::MAX && inside(self.values[t]) {
                        label[t] = id;
                        stack.push(t);
                    }
                }
            }
            sizes.push(size);
        }
        (label, sizes)
    }
}

/// Removes solid components that contain no anchor sample (all but the largest when there are no
/// anchors), then fills void components that do not reach the field's border (closed cavities).
pub fn clean_components(field: &mut Field, iso: f64, anchors: Option<&[bool]>, remove_islands: bool, fill_cavities: bool) -> Cleanup {
    let mut report = Cleanup::default();
    let (label, sizes) = field.components(iso, true);
    report.solid_components = sizes.len();
    if remove_islands && sizes.len() > 1 {
        let mut keep = vec![false; sizes.len()];
        if let Some(a) = anchors {
            for (s, &is_anchor) in a.iter().enumerate() {
                if is_anchor && label[s] != u32::MAX {
                    keep[label[s] as usize] = true;
                }
            }
        }
        if !keep.iter().any(|&k| k) {
            let largest = (0..sizes.len()).max_by_key(|&c| sizes[c]).unwrap();
            keep[largest] = true;
        }
        for (s, &l) in label.iter().enumerate() {
            if l != u32::MAX && !keep[l as usize] {
                field.values[s] = 0.0;
                report.removed_samples += 1;
            }
        }
        report.removed_components = keep.iter().filter(|&&k| !k).count();
    }
    if fill_cavities {
        let (label, sizes) = field.components(iso, false);
        // Sample 0 is a corner of the zero border: its component is the outside.
        let outside = label[0];
        for (s, &l) in label.iter().enumerate() {
            if l != u32::MAX && l != outside {
                field.values[s] = 1.0;
                report.filled_samples += 1;
            }
        }
        report.cavities_filled = sizes.len().saturating_sub(1);
    }
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn islands_removed_cavities_filled() {
        let dims = [10, 10, 10];
        let mut cells = vec![0.0; 1000];
        let c = |i: usize, j: usize, k: usize| i + 10 * (j + 10 * k);
        // A hollow box (shell 1 thick, cavity inside) and a separate speck.
        for k in 1..8 {
            for j in 1..8 {
                for i in 1..8 {
                    cells[c(i, j, k)] = 1.0;
                }
            }
        }
        for k in 3..6 {
            for j in 3..6 {
                for i in 3..6 {
                    cells[c(i, j, k)] = 0.0;
                }
            }
        }
        cells[c(9, 9, 9)] = 1.0;
        let mut f = Field::from_cells(dims, 1.0, [0.0; 3], &cells, 1);
        let r = clean_components(&mut f, 0.5, None, true, true);
        assert_eq!(r.solid_components, 2);
        assert_eq!(r.removed_components, 1);
        assert_eq!(r.removed_samples, 1);
        assert_eq!(r.cavities_filled, 1);
        assert_eq!(r.filled_samples, 27);
        assert_eq!(f.components(0.5, true).1.len(), 1);
        assert_eq!(f.components(0.5, false).1.len(), 1);
    }
}
