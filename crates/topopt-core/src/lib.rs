//! Density-based topology optimization (SIMP) on a structured 3D Hex8 grid.
//!
//! - [`grid`]: indexing of the structured grid.
//! - [`ke`]: the Hex8 element stiffness.
//! - [`problem`]: what to optimize (supports, loads, passive regions, material).
//! - [`fem`]: matrix-free stiffness operators.
//! - [`mg`] + [`pcg`]: the multigrid-preconditioned conjugate-gradient solver.
//! - [`filter`]: density filter and Heaviside projection.
//! - [`opt`]: the optimization loop (optimality criteria).
//!
//! Pure Rust, no I/O: it runs natively and inside the SolidCraft plugin sandbox
//! (wasm32-unknown-unknown). Single-threaded by default; the `parallel` feature spreads the heavy
//! loops over rayon's threads ([`par`]).

pub mod fem;
pub mod filter;
pub mod grid;
pub mod ke;
pub mod mg;
pub mod opt;
pub mod par;
pub mod pcg;
pub mod problem;

pub use grid::Grid;
pub use opt::{Iteration, Optimizer, Outcome, Settings};
pub use problem::{ElementKind, LoadCase, Material, NodeFix, NodeForce, Problem};
