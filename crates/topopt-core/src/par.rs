//! Data-parallel loops: rayon with the `parallel` feature (the native solver), plain loops
//! without it (the plugin sandbox has one thread).

/// Calls `f(index, chunk)` for consecutive `len`-sized chunks of `out`, in parallel when enabled.
pub(crate) fn chunks<T: Send>(out: &mut [T], len: usize, f: impl Fn(usize, &mut [T]) + Sync + Send) {
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        out.par_chunks_mut(len.max(1)).enumerate().for_each(|(i, c)| f(i, c));
    }
    #[cfg(not(feature = "parallel"))]
    {
        out.chunks_mut(len.max(1)).enumerate().for_each(|(i, c)| f(i, c));
    }
}

/// Calls `f(i)` for `i` in `0..n`, in parallel when enabled.
#[cfg_attr(not(feature = "parallel"), allow(dead_code))]
pub(crate) fn each(n: usize, f: impl Fn(usize) + Sync + Send) {
    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        (0..n).into_par_iter().for_each(f);
    }
    #[cfg(not(feature = "parallel"))]
    {
        (0..n).for_each(f);
    }
}

/// Number of worker threads (1 without the feature).
pub fn threads() -> usize {
    #[cfg(feature = "parallel")]
    {
        rayon::current_num_threads()
    }
    #[cfg(not(feature = "parallel"))]
    {
        1
    }
}

/// A raw pointer that may cross threads; callers guarantee disjoint writes.
#[derive(Clone, Copy)]
pub(crate) struct SharedMut(pub *mut f64);
unsafe impl Send for SharedMut {}
unsafe impl Sync for SharedMut {}

/// Vectors shorter than this are processed on one thread (spawning costs more than it saves).
const MIN_PARALLEL: usize = 1 << 14;

/// `out[i] = f(i)` for every `i`, in parallel blocks when enabled and worth it.
pub(crate) fn fill(out: &mut [f64], f: impl Fn(usize) -> f64 + Sync + Send) {
    if threads() > 1 && out.len() >= MIN_PARALLEL {
        let block = 4096;
        chunks(out, block, |b, c| {
            let base = b * block;
            for (i, v) in c.iter_mut().enumerate() {
                *v = f(base + i);
            }
        });
    } else {
        for (i, v) in out.iter_mut().enumerate() {
            *v = f(i);
        }
    }
}

/// `Σ a[i]·b[i]`, in parallel when enabled and worth it.
pub(crate) fn dot(a: &[f64], b: &[f64]) -> f64 {
    #[cfg(feature = "parallel")]
    if a.len() >= MIN_PARALLEL {
        use rayon::prelude::*;
        return a.par_chunks(4096).zip(b.par_chunks(4096)).map(|(x, y)| x.iter().zip(y).map(|(p, q)| p * q).sum::<f64>()).sum();
    }
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// `out[i] = f(i, out[i])` for every `i`, in place, in parallel blocks when enabled and worth it.
pub(crate) fn update(out: &mut [f64], f: impl Fn(usize, f64) -> f64 + Sync + Send) {
    if threads() > 1 && out.len() >= MIN_PARALLEL {
        let block = 4096;
        chunks(out, block, |b, c| {
            let base = b * block;
            for (i, v) in c.iter_mut().enumerate() {
                *v = f(base + i, *v);
            }
        });
    } else {
        for (i, v) in out.iter_mut().enumerate() {
            *v = f(i, *v);
        }
    }
}
