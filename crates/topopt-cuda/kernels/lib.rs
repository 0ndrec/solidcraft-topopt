//! TopOpt GPU kernels: the matrix-free Hex8 stiffness operator, the multigrid V-cycle pieces and
//! the vector operations of conjugate gradients, on the structured grid of `topopt-core`.
//!
//! Built for `nvptx64-nvidia-cuda` with a nightly toolchain by `scripts/build-ptx.sh`; the PTX is
//! checked in (`../src/kernels.ptx`) and loaded by the CUDA driver at run time, so running the
//! solver needs no toolkit and no nightly.
//!
//! Layout (as `topopt_core::Grid`): node `(i, j, k)` has index `i + sx·(j + sy·k)` with
//! `sx = nx + 1`, `sy = ny + 1`, and dofs `3·node + d`; element `(i, j, k)` has index
//! `i + nx·(j + ny·k)`; an element's corners follow `CORNERS`. Element matrices are 24×24 with
//! entry (row `r`, column `c`) at `c·24 + r`, as `topopt-core` reads them. Supports: `fix[node]`
//! has bit `d` set when dof `d` is fixed; operators take fixed dofs as zero on input and are the
//! identity on them.
//!
//! One thread per node: it gathers its three rows from the (up to) eight elements around it, so
//! nothing is written twice and no atomics are needed. Reductions use a global scratch area per
//! block (no shared memory) and a fixed number of blocks, so results are deterministic.

#![no_std]
#![feature(abi_ptx, stdarch_nvptx, core_float_math, asm_experimental_arch)]
#![allow(internal_features, clippy::missing_safety_doc, clippy::too_many_arguments)]

use core::arch::asm;
use core::arch::nvptx::{_block_dim_x, _block_idx_x, _grid_dim_x, _syncthreads, _thread_idx_x};
use core::ptr::{read_volatile, write_volatile};

#[panic_handler]
fn panic(_: &core::panic::PanicInfo) -> ! {
    unsafe { core::arch::nvptx::trap() }
}

/// Corner offsets in `CORNERS` order: bottom face counter-clockwise, then top.
const CX: [u32; 8] = [0, 1, 1, 0, 0, 1, 1, 0];
const CY: [u32; 8] = [0, 0, 1, 1, 0, 0, 1, 1];
const CZ: [u32; 8] = [0, 0, 0, 0, 1, 1, 1, 1];

/// Runs `$body` eight times with `$c` a constant 0…7: the loops over an element's corners are
/// unrolled, so every matrix index is a compile-time constant.
macro_rules! unroll8 {
    ($c:ident, $body:block) => {
        { const $c: usize = 0; $body }
        { const $c: usize = 1; $body }
        { const $c: usize = 2; $body }
        { const $c: usize = 3; $body }
        { const $c: usize = 4; $body }
        { const $c: usize = 5; $body }
        { const $c: usize = 6; $body }
        { const $c: usize = 7; $body }
    };
}

/// Runs `$body` three times with `$c` a constant 0…2.
macro_rules! unroll3 {
    ($c:ident, $body:block) => {
        { const $c: usize = 0; $body }
        { const $c: usize = 1; $body }
        { const $c: usize = 2; $body }
    };
}

/// Loads through the read-only data cache (`ld.global.nc`, CUDA's `__ldg`): on Pascal, plain
/// global loads are cached in L2 only, and the operator reads every matrix entry and neighbour
/// value several times.
macro_rules! ldg {
    ($ty:literal, $class:ident, $t:ty, $p:expr) => {{
        let v: $t;
        let p = $p;
        unsafe {
            asm!(
                concat!("{{ .reg .b64 %g; cvta.to.global.u64 %g, {p}; ld.global.nc.", $ty, " {v}, [%g]; }}"),
                v = out($class) v,
                p = in(reg64) p,
                options(pure, readonly, nostack, preserves_flags)
            )
        };
        v
    }};
}

trait Real: Copy + core::ops::Add<Output = Self> + core::ops::Mul<Output = Self> {
    const ZERO: Self;
    fn fma(a: Self, b: Self, c: Self) -> Self;
    /// `*p` through the read-only cache; `p` must point to device memory nothing writes during
    /// the kernel.
    unsafe fn ldg(p: *const Self) -> Self;
}

impl Real for f32 {
    const ZERO: f32 = 0.0;
    #[inline(always)]
    fn fma(a: f32, b: f32, c: f32) -> f32 {
        core::f32::math::mul_add(a, b, c)
    }
    #[inline(always)]
    unsafe fn ldg(p: *const f32) -> f32 {
        ldg!("f32", reg32, f32, p)
    }
}

impl Real for f64 {
    const ZERO: f64 = 0.0;
    #[inline(always)]
    fn fma(a: f64, b: f64, c: f64) -> f64 {
        core::f64::math::mul_add(a, b, c)
    }
    #[inline(always)]
    unsafe fn ldg(p: *const f64) -> f64 {
        ldg!("f64", reg64, f64, p)
    }
}

#[inline(always)]
unsafe fn ldg_u8(p: *const u8) -> u32 {
    ldg!("u8", reg32, u32, p)
}

#[inline(always)]
fn global_id() -> u32 {
    unsafe { (_block_idx_x() as u32) * (_block_dim_x() as u32) + _thread_idx_x() as u32 }
}

/// Node coordinates of node `n` on an `nx × ny × nz`-element grid, or `None` past the last node.
#[inline(always)]
fn node_ijk(n: u32, nx: u32, ny: u32, nz: u32) -> Option<(u32, u32, u32)> {
    let (sx, sy) = (nx + 1, ny + 1);
    if n >= sx * sy * (nz + 1) {
        return None;
    }
    let t = n / sx;
    Some((n - t * sx, t % sy, t / sy))
}

/// Element `(i, j, k)` of which this node is corner `a`, if it exists.
#[inline(always)]
fn element_at(i: u32, j: u32, k: u32, a: usize, nx: u32, ny: u32, nz: u32) -> Option<(u32, u32, u32)> {
    if i < CX[a] || j < CY[a] || k < CZ[a] {
        return None;
    }
    let (ei, ej, ek) = (i - CX[a], j - CY[a], k - CZ[a]);
    if ei >= nx || ej >= ny || ek >= nz { None } else { Some((ei, ej, ek)) }
}

/// A node's three values, fixed dofs read as zero.
#[inline(always)]
unsafe fn load_node<T: Real>(x: *const T, fix: *const u8, n: usize) -> [T; 3] {
    unsafe {
        let m = ldg_u8(fix.add(n));
        [
            if m & 1 != 0 { T::ZERO } else { T::ldg(x.add(3 * n)) },
            if m & 2 != 0 { T::ZERO } else { T::ldg(x.add(3 * n + 1)) },
            if m & 4 != 0 { T::ZERO } else { T::ldg(x.add(3 * n + 2)) },
        ]
    }
}

#[inline(always)]
unsafe fn store_node<T: Real>(x: *const T, fix: *const u8, n: usize, acc: [T; 3], y: *mut T) {
    unsafe {
        let m = ldg_u8(fix.add(n));
        for d in 0..3 {
            *y.add(3 * n + d) = if (m >> d) & 1 != 0 { *x.add(3 * n + d) } else { acc[d] };
        }
    }
}

/// The 27 nodes around node `n = (i, j, k)`, at index `dx + 3·(dy + 3·dz)` for offsets
/// `dx − 1, dy − 1, dz − 1`: each is loaded once, though up to eight elements share it. Fixed dofs
/// and nodes off the grid read as zero.
#[inline(always)]
unsafe fn neighbourhood<T: Real>(x: *const T, fix: *const u8, n: usize, ijk: (u32, u32, u32), dims: (u32, u32, u32), sx: usize, sxy: usize) -> [[T; 3]; 27] {
    let ((i, j, k), (nx, ny, nz)) = (ijk, dims);
    let mut v = [[T::ZERO; 3]; 27];
    unroll3!(DZ, {
        unroll3!(DY, {
            unroll3!(DX, {
                let inside = (DX != 0 || i > 0) && (DX != 2 || i < nx) && (DY != 0 || j > 0) && (DY != 2 || j < ny) && (DZ != 0 || k > 0) && (DZ != 2 || k < nz);
                if inside {
                    let m = n + DX + sx * DY + sxy * DZ - (1 + sx + sxy);
                    v[DX + 3 * (DY + 3 * DZ)] = unsafe { load_node(x, fix, m) };
                }
            });
        });
    });
    v
}

/// `t += M[rows of corner A] · x_e` for the element of which this node is corner `A`: `m(q)` is
/// entry `q` of the element's matrix, (row r, column c) at `q = c·24 + r`; `node(o)` the values
/// of the node at offset `o` in the [`neighbourhood`] numbering.
#[inline(always)]
fn element_rows<T: Real, const A: usize>(m: impl Fn(usize) -> T, node: impl Fn(usize) -> [T; 3], t: &mut [T; 3]) {
    unroll8!(B, {
        let v = node((CX[B] + 1 - CX[A]) as usize + 3 * ((CY[B] + 1 - CY[A]) as usize + 3 * (CZ[B] + 1 - CZ[A]) as usize));
        for d in 0..3 {
            let q = 3 * B * 24 + 3 * A + d;
            t[d] = T::fma(m(q), v[0], t[d]);
            t[d] = T::fma(m(q + 24), v[1], t[d]);
            t[d] = T::fma(m(q + 48), v[2], t[d]);
        }
    });
}

/// The node at offset `o` (as in [`neighbourhood`]) from node `n`, loaded where it is needed.
#[inline(always)]
unsafe fn neighbour<T: Real>(x: *const T, fix: *const u8, n: usize, o: usize, sx: usize, sxy: usize) -> [T; 3] {
    unsafe { load_node(x, fix, n + o % 3 + sx * ((o / 3) % 3) + sxy * (o / 9) - (1 + sx + sxy)) }
}

/// `y = K·x` on the fine level: element `e`'s matrix is `emod[e] · ke`, with `ke(q)` its
/// entries. With `GATHER`, the 27 neighbours are loaded once up front (fewer loads, more
/// registers: it pays in single precision only).
#[inline(always)]
unsafe fn apply_fine<T: Real, const GATHER: bool>(nx: u32, ny: u32, nz: u32, ke: impl Fn(usize) -> T + Copy, emod: *const T, fix: *const u8, x: *const T, y: *mut T) {
    let n = global_id();
    let Some((i, j, k)) = node_ijk(n, nx, ny, nz) else { return };
    let (sx, sxy) = ((nx + 1) as usize, ((nx + 1) * (ny + 1)) as usize);
    let nb = if GATHER { unsafe { neighbourhood(x, fix, n as usize, (i, j, k), (nx, ny, nz), sx, sxy) } } else { [[T::ZERO; 3]; 27] };
    let mut acc = [T::ZERO; 3];
    unroll8!(A, {
        if let Some((ei, ej, ek)) = element_at(i, j, k, A, nx, ny, nz) {
            let mut t = [T::ZERO; 3];
            if GATHER {
                element_rows::<T, A>(ke, |o| nb[o], &mut t);
            } else {
                element_rows::<T, A>(ke, |o| unsafe { neighbour(x, fix, n as usize, o, sx, sxy) }, &mut t);
            }
            let s = unsafe { T::ldg(emod.add((ei + nx * (ej + ny * ek)) as usize)) };
            for d in 0..3 {
                acc[d] = T::fma(s, t[d], acc[d]);
            }
        }
    });
    unsafe { store_node(x, fix, n as usize, acc, y) };
}

/// An element matrix passed by value: kernel parameters live in the constant bank, so the
/// multiply-adds take its entries as operands, with no loads at all.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct Ke32(pub [f32; 576]);

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn apply_fine_f32(nx: u32, ny: u32, nz: u32, ke: Ke32, emod: *const f32, fix: *const u8, x: *const f32, y: *mut f32) {
    unsafe { apply_fine::<f32, true>(nx, ny, nz, |q| ke.0[q], emod, fix, x, y) }
}

/// In double precision the matrix (4.5 KB) exceeds the 4 KB of kernel parameters; it is read
/// through the read-only cache instead.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn apply_fine_f64(nx: u32, ny: u32, nz: u32, ke: *const f64, emod: *const f64, fix: *const u8, x: *const f64, y: *mut f64) {
    unsafe { apply_fine::<f64, false>(nx, ny, nz, |q| f64::ldg(ke.add(q)), emod, fix, x, y) }
}

/// `y = K·x` on a coarse level: one 24×24 matrix per element, stored entry-major
/// (`mats[q·nel + e]`) so that neighbouring threads read neighbouring addresses.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn apply_own_f32(nx: u32, ny: u32, nz: u32, mats: *const f32, fix: *const u8, x: *const f32, y: *mut f32) {
    let n = global_id();
    let Some((i, j, k)) = node_ijk(n, nx, ny, nz) else { return };
    let (sx, sxy) = ((nx + 1) as usize, ((nx + 1) * (ny + 1)) as usize);
    let nel = (nx * ny * nz) as usize;
    let mut acc = [0.0f32; 3];
    unroll8!(A, {
        if let Some((ei, ej, ek)) = element_at(i, j, k, A, nx, ny, nz) {
            let m = unsafe { mats.add((ei + nx * (ej + ny * ek)) as usize) };
            element_rows::<f32, A>(|q| unsafe { f32::ldg(m.add(q * nel)) }, |o| unsafe { neighbour(x, fix, n as usize, o, sx, sxy) }, &mut acc);
        }
    });
    unsafe { store_node(x, fix, n as usize, acc, y) };
}

/// Writes `1 / diagonal` (1 on fixed dofs and where the diagonal is 0).
#[inline(always)]
unsafe fn write_inv_diag(fix: *const u8, n: usize, diag: [f64; 3], out: *mut f32) {
    unsafe {
        let m = *fix.add(n);
        for d in 0..3 {
            *out.add(3 * n + d) = if (m >> d) & 1 != 0 || diag[d] == 0.0 { 1.0 } else { (1.0 / diag[d]) as f32 };
        }
    }
}

/// Inverse diagonal of the fine level (`emod[e] · ke`).
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn inv_diag_fine(nx: u32, ny: u32, nz: u32, ke: *const f32, emod: *const f32, fix: *const u8, out: *mut f32) {
    let n = global_id();
    let Some((i, j, k)) = node_ijk(n, nx, ny, nz) else { return };
    let mut diag = [0.0f64; 3];
    for a in 0..8 {
        let Some((ei, ej, ek)) = element_at(i, j, k, a, nx, ny, nz) else { continue };
        let s = unsafe { *emod.add((ei + nx * (ej + ny * ek)) as usize) } as f64;
        for d in 0..3 {
            diag[d] += s * unsafe { *ke.add((3 * a + d) * 25) } as f64;
        }
    }
    unsafe { write_inv_diag(fix, n as usize, diag, out) };
}

/// Inverse diagonal of a coarse level (own element matrices, entry-major).
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn inv_diag_own(nx: u32, ny: u32, nz: u32, mats: *const f32, fix: *const u8, out: *mut f32) {
    let n = global_id();
    let Some((i, j, k)) = node_ijk(n, nx, ny, nz) else { return };
    let nel = (nx * ny * nz) as usize;
    let mut diag = [0.0f64; 3];
    for a in 0..8 {
        let Some((ei, ej, ek)) = element_at(i, j, k, a, nx, ny, nz) else { continue };
        let m = unsafe { mats.add((ei + nx * (ej + ny * ek)) as usize) };
        for d in 0..3 {
            diag[d] += unsafe { *m.add((3 * a + d) * 25 * nel) } as f64;
        }
    }
    unsafe { write_inv_diag(fix, n as usize, diag, out) };
}

/// Damped Jacobi: `x += ω·D⁻¹·(b − ax)`; with `first != 0`, `x = ω·D⁻¹·b` (x was zero).
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn jacobi(n: u32, omega: f32, invd: *const f32, b: *const f32, ax: *const f32, x: *mut f32, first: u32) {
    let i = global_id();
    if i >= n {
        return;
    }
    let i = i as usize;
    unsafe {
        let w = omega * *invd.add(i);
        *x.add(i) = if first != 0 { w * *b.add(i) } else { f32::fma(w, *b.add(i) - *ax.add(i), *x.add(i)) };
    }
}

/// `coarse = Pᵀ·(b − ax)` with `P` masked by both fixed sets: coarse node `I` collects fine node
/// `2I` with weight 1 and `2I ± 1` with weight ½, per axis.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn restrict_residual(
    cnx: u32,
    cny: u32,
    cnz: u32,
    cfix: *const u8,
    fnx: u32,
    fny: u32,
    fnz: u32,
    ffix: *const u8,
    b: *const f32,
    ax: *const f32,
    coarse: *mut f32,
) {
    let cn = global_id();
    let Some((ic, jc, kc)) = node_ijk(cn, cnx, cny, cnz) else { return };
    let (fsx, fsxy) = (fnx + 1, (fnx + 1) * (fny + 1));
    let mut v = [0.0f32; 3];
    for dz in 0..3u32 {
        let fk = 2 * kc + dz;
        if fk == 0 || fk > fnz + 1 {
            continue;
        }
        let fk = fk - 1;
        let wz = if dz == 1 { 1.0 } else { 0.5 };
        for dy in 0..3u32 {
            let fj = 2 * jc + dy;
            if fj == 0 || fj > fny + 1 {
                continue;
            }
            let fj = fj - 1;
            let wy = if dy == 1 { 1.0 } else { 0.5 };
            for dx in 0..3u32 {
                let fi = 2 * ic + dx;
                if fi == 0 || fi > fnx + 1 {
                    continue;
                }
                let fi = fi - 1;
                let w = wz * wy * if dx == 1 { 1.0 } else { 0.5 };
                let fnode = (fi + fsx * fj + fsxy * fk) as usize;
                let m = unsafe { *ffix.add(fnode) };
                for d in 0..3 {
                    if (m >> d) & 1 == 0 {
                        v[d] = f32::fma(w, unsafe { *b.add(3 * fnode + d) - *ax.add(3 * fnode + d) }, v[d]);
                    }
                }
            }
        }
    }
    let cn = cn as usize;
    let m = unsafe { *cfix.add(cn) };
    for d in 0..3 {
        unsafe { *coarse.add(3 * cn + d) = if (m >> d) & 1 != 0 { 0.0 } else { v[d] } };
    }
}

/// `fine += P·coarse` (masked by both fixed sets): an even fine coordinate takes coarse `i/2`, an
/// odd one `(i−1)/2` and `(i+1)/2` with weight ½ each.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn prolong_add(
    cnx: u32,
    cny: u32,
    _cnz: u32,
    cfix: *const u8,
    fnx: u32,
    fny: u32,
    fnz: u32,
    ffix: *const u8,
    coarse: *const f32,
    fine: *mut f32,
) {
    let fnode = global_id();
    let Some((i, j, k)) = node_ijk(fnode, fnx, fny, fnz) else { return };
    let (csx, csxy) = (cnx + 1, (cnx + 1) * (cny + 1));
    let (nxs, nys, nzs) = ((i & 1) + 1, (j & 1) + 1, (k & 1) + 1);
    let (wx, wy, wz) = (1.0 / nxs as f32, 1.0 / nys as f32, 1.0 / nzs as f32);
    let mut v = [0.0f32; 3];
    for c in 0..nzs {
        for bb in 0..nys {
            for aa in 0..nxs {
                let cn = ((i >> 1) + aa + csx * ((j >> 1) + bb) + csxy * ((k >> 1) + c)) as usize;
                let m = unsafe { *cfix.add(cn) };
                let w = wx * wy * wz;
                for d in 0..3 {
                    if (m >> d) & 1 == 0 {
                        v[d] = f32::fma(w, unsafe { *coarse.add(3 * cn + d) }, v[d]);
                    }
                }
            }
        }
    }
    let fnode = fnode as usize;
    let m = unsafe { *ffix.add(fnode) };
    for d in 0..3 {
        if (m >> d) & 1 == 0 {
            unsafe { *fine.add(3 * fnode + d) += v[d] };
        }
    }
}

// ---- Reductions: dot products accumulated in double; each block reduces through its slice of a
// global scratch area (`blockDim` doubles per block) and writes one partial sum; one block then
// sums the partials into `out[slot]`.

#[inline(always)]
unsafe fn block_sum(v: f64, scratch: *mut f64) -> f64 {
    unsafe {
        let t = _thread_idx_x() as usize;
        let s = scratch.add(_block_idx_x() as usize * _block_dim_x() as usize);
        write_volatile(s.add(t), v);
        _syncthreads();
        let mut w = _block_dim_x() as usize / 2;
        while w > 0 {
            if t < w {
                write_volatile(s.add(t), read_volatile(s.add(t)) + read_volatile(s.add(t + w)));
            }
            _syncthreads();
            w /= 2;
        }
        read_volatile(s)
    }
}

#[inline(always)]
unsafe fn dot_partial<T: Copy + Into<f64>>(n: u32, a: *const T, b: *const T, scratch: *mut f64, partial: *mut f64) {
    unsafe {
        let stride = (_grid_dim_x() as u32) * (_block_dim_x() as u32);
        let mut i = global_id();
        let mut acc = 0.0f64;
        while i < n {
            acc += (*a.add(i as usize)).into() * (*b.add(i as usize)).into();
            i += stride;
        }
        let s = block_sum(acc, scratch);
        if _thread_idx_x() == 0 {
            *partial.add(_block_idx_x() as usize) = s;
        }
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn dot_f32(n: u32, a: *const f32, b: *const f32, scratch: *mut f64, partial: *mut f64) {
    unsafe { dot_partial(n, a, b, scratch, partial) }
}

#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn dot_f64(n: u32, a: *const f64, b: *const f64, scratch: *mut f64, partial: *mut f64) {
    unsafe { dot_partial(n, a, b, scratch, partial) }
}

/// `out[slot] = Σ partial[0..m]` (one block).
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn sum_partials(m: u32, partial: *const f64, scratch: *mut f64, out: *mut f64, slot: u32) {
    unsafe {
        let mut i = _thread_idx_x() as u32;
        let mut acc = 0.0f64;
        while i < m {
            acc += *partial.add(i as usize);
            i += _block_dim_x() as u32;
        }
        let s = block_sum(acc, scratch);
        if _thread_idx_x() == 0 {
            *out.add(slot as usize) = s;
        }
    }
}

// ---- Conjugate gradients with the scalars kept on the device.

/// `num / den`, or 0 when `den` is not positive (CG's breakdown; the host sees it and stops).
#[inline(always)]
fn ratio(num: f64, den: f64) -> f32 {
    if den > 0.0 { (num / den) as f32 } else { 0.0 }
}

/// `α = scal[num] / scal[den]`; `x += α·p`; `r −= α·q`.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn cg_xr(n: u32, scal: *const f64, num: u32, den: u32, x: *mut f32, r: *mut f32, p: *const f32, q: *const f32) {
    let i = global_id();
    if i >= n {
        return;
    }
    let i = i as usize;
    unsafe {
        let alpha = ratio(*scal.add(num as usize), *scal.add(den as usize));
        *x.add(i) = f32::fma(alpha, *p.add(i), *x.add(i));
        *r.add(i) = f32::fma(-alpha, *q.add(i), *r.add(i));
    }
}

/// `β = scal[num] / scal[den]`; `p = z + β·p`.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn cg_p(n: u32, scal: *const f64, num: u32, den: u32, p: *mut f32, z: *const f32) {
    let i = global_id();
    if i >= n {
        return;
    }
    let i = i as usize;
    unsafe {
        let beta = ratio(*scal.add(num as usize), *scal.add(den as usize));
        *p.add(i) = f32::fma(beta, *p.add(i), *z.add(i));
    }
}

// ---- Mixed precision: the residual in double, the correction in single.

/// `r = f − ax`.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn residual_f64(n: u32, f: *const f64, ax: *const f64, r: *mut f64) {
    let i = global_id();
    if i >= n {
        return;
    }
    let i = i as usize;
    unsafe { *r.add(i) = *f.add(i) - *ax.add(i) };
}

/// `out = r · scale`, rounded to single.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn to_f32(n: u32, r: *const f64, scale: f64, out: *mut f32) {
    let i = global_id();
    if i >= n {
        return;
    }
    let i = i as usize;
    unsafe { *out.add(i) = (*r.add(i) * scale) as f32 };
}

/// `u += scale · d`.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn add_f32_to_f64(n: u32, scale: f64, d: *const f32, u: *mut f64) {
    let i = global_id();
    if i >= n {
        return;
    }
    let i = i as usize;
    unsafe { *u.add(i) = f64::fma(scale, *d.add(i) as f64, *u.add(i)) };
}

// ---- Power iteration for the largest eigenvalue of D⁻¹·A.

/// `v *= s`.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn scale_f32(n: u32, s: f32, v: *mut f32) {
    let i = global_id();
    if i >= n {
        return;
    }
    unsafe { *v.add(i as usize) *= s };
}

/// `av = fixed ? 0 : av · invd`, per dof.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn scale_inv_diag(nnodes: u32, fix: *const u8, invd: *const f32, av: *mut f32) {
    let n = global_id();
    if n >= nnodes {
        return;
    }
    let n = n as usize;
    unsafe {
        let m = *fix.add(n);
        for d in 0..3 {
            let i = 3 * n + d;
            *av.add(i) = if (m >> d) & 1 != 0 { 0.0 } else { *av.add(i) * *invd.add(i) };
        }
    }
}

// ---- Coarse operators.

/// The first coarse level's element matrices, entry-major (`out[q·nel + e]`): coarse element
/// `e` is the sum over its eight children `c` of `emod[child] · g[c]`, with `g[c] = P_cᵀ·ke·P_c`
/// (`g[576·c + q]`), exactly the Galerkin product the CPU forms.
#[unsafe(no_mangle)]
pub unsafe extern "ptx-kernel" fn galerkin_first(cnx: u32, cny: u32, cnz: u32, g: *const f32, emod: *const f32, out: *mut f32) {
    let e = global_id();
    let nel = cnx * cny * cnz;
    if e >= nel {
        return;
    }
    let t = e / cnx;
    let (ci, cj, ck) = (e - t * cnx, t % cny, t / cny);
    let (fnx, fny) = (2 * cnx, 2 * cny);
    let mut s = [0.0f32; 8];
    unroll8!(C, {
        s[C] = unsafe { f32::ldg(emod.add((2 * ci + CX[C] + fnx * (2 * cj + CY[C] + fny * (2 * ck + CZ[C]))) as usize)) };
    });
    let (e, nel) = (e as usize, nel as usize);
    for q in 0..576 {
        let mut v = 0.0f32;
        unroll8!(C, {
            v = f32::fma(s[C], unsafe { f32::ldg(g.add(576 * C + q)) }, v);
        });
        unsafe { *out.add(q * nel + e) = v };
    }
}
