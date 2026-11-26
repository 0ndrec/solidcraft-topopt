//! Structured hexahedral grid: `nx × ny × nz` cubic elements of side `h` (mm).
//!
//! Node `(i, j, k)`, `0 ≤ i ≤ nx`, has index `i + (nx+1)·(j + (ny+1)·k)` and degrees of freedom
//! `3·node + d` for `d ∈ {x, y, z}`. Element `(i, j, k)` has index `i + nx·(j + ny·k)`; its eight
//! corners are in the order of [`CORNERS`], which [`crate::ke`] follows.

use serde::{Deserialize, Serialize};

/// Corner offsets of an element in local-node order: bottom face counter-clockwise, then top.
pub const CORNERS: [[usize; 3]; 8] =
    [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Grid {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    /// Element edge length, mm.
    pub h: f64,
    /// World position of node (0, 0, 0), mm.
    #[serde(default)]
    pub origin: [f64; 3],
}

impl Grid {
    pub fn new(nx: usize, ny: usize, nz: usize, h: f64) -> Self {
        assert!(nx > 0 && ny > 0 && nz > 0 && h > 0.0, "grid must be non-empty with h > 0");
        Grid { nx, ny, nz, h, origin: [0.0; 3] }
    }

    pub fn nel(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    pub fn nnodes(&self) -> usize {
        (self.nx + 1) * (self.ny + 1) * (self.nz + 1)
    }

    pub fn ndof(&self) -> usize {
        3 * self.nnodes()
    }

    #[inline]
    pub fn node(&self, i: usize, j: usize, k: usize) -> usize {
        i + (self.nx + 1) * (j + (self.ny + 1) * k)
    }

    #[inline]
    pub fn node_ijk(&self, n: usize) -> [usize; 3] {
        let (sx, sy) = (self.nx + 1, self.ny + 1);
        [n % sx, (n / sx) % sy, n / (sx * sy)]
    }

    #[inline]
    pub fn elem(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.nx * (j + self.ny * k)
    }

    #[inline]
    pub fn elem_ijk(&self, e: usize) -> [usize; 3] {
        [e % self.nx, (e / self.nx) % self.ny, e / (self.nx * self.ny)]
    }

    /// The element's eight node indices in [`CORNERS`] order.
    #[inline]
    pub fn elem_nodes(&self, i: usize, j: usize, k: usize) -> [usize; 8] {
        let n0 = self.node(i, j, k);
        let sx = self.nx + 1;
        let sxy = sx * (self.ny + 1);
        [n0, n0 + 1, n0 + 1 + sx, n0 + sx, n0 + sxy, n0 + 1 + sxy, n0 + 1 + sx + sxy, n0 + sx + sxy]
    }

    /// The element's 24 dof indices: node-major, x/y/z within a node.
    #[inline]
    pub fn elem_dofs(&self, i: usize, j: usize, k: usize) -> [usize; 24] {
        let nodes = self.elem_nodes(i, j, k);
        let mut d = [0usize; 24];
        for (a, &n) in nodes.iter().enumerate() {
            d[3 * a] = 3 * n;
            d[3 * a + 1] = 3 * n + 1;
            d[3 * a + 2] = 3 * n + 2;
        }
        d
    }

    /// World position of a node, mm.
    pub fn node_pos(&self, n: usize) -> [f64; 3] {
        let [i, j, k] = self.node_ijk(n);
        [self.origin[0] + i as f64 * self.h, self.origin[1] + j as f64 * self.h, self.origin[2] + k as f64 * self.h]
    }

    /// World position of an element's centre, mm.
    pub fn elem_center(&self, e: usize) -> [f64; 3] {
        let [i, j, k] = self.elem_ijk(e);
        [
            self.origin[0] + (i as f64 + 0.5) * self.h,
            self.origin[1] + (j as f64 + 0.5) * self.h,
            self.origin[2] + (k as f64 + 0.5) * self.h,
        ]
    }

    /// The grid with twice the element size, when every dimension is even.
    pub fn coarsen(&self) -> Option<Grid> {
        if self.nx % 2 == 0 && self.ny % 2 == 0 && self.nz % 2 == 0 {
            Some(Grid { nx: self.nx / 2, ny: self.ny / 2, nz: self.nz / 2, h: self.h * 2.0, origin: self.origin })
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexing_round_trips() {
        let g = Grid::new(4, 3, 2, 1.0);
        assert_eq!(g.nel(), 24);
        assert_eq!(g.nnodes(), 5 * 4 * 3);
        for e in 0..g.nel() {
            let [i, j, k] = g.elem_ijk(e);
            assert_eq!(g.elem(i, j, k), e);
            let nodes = g.elem_nodes(i, j, k);
            for (a, c) in CORNERS.iter().enumerate() {
                assert_eq!(nodes[a], g.node(i + c[0], j + c[1], k + c[2]));
            }
        }
        for n in 0..g.nnodes() {
            let [i, j, k] = g.node_ijk(n);
            assert_eq!(g.node(i, j, k), n);
        }
        assert_eq!(g.coarsen(), None);
        assert_eq!(Grid::new(4, 2, 2, 1.0).coarsen().unwrap().nel(), 2);
    }
}
