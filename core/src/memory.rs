//! Blocked Array

use std::ops::{Index, IndexMut};

/// Stores 2D arrays in a blocked memory layout.
#[derive(Clone)]
pub struct BlockedArray<T> {
    /// Array data.
    data: Vec<T>,

    /// Number of blocks in u-dimension.
    u_blocks: usize,

    /// Size in u-dimension.
    u_res: usize,

    /// Size in v-dimension.
    v_res: usize,
}

impl<T> BlockedArray<T>
where
    T: Copy + Default,
{
    /// Create a new `BlockedArray<T>` with default values.
    ///
    /// * `u_res`          - Size in u-dimension.
    /// * `v_res`          - Size in v-dimension.
    /// * `log_block_size` - Logarithm base 2 value used to ensure block sizes
    ///                      will be powers of 2.
    pub fn new(u_res: usize, v_res: usize) -> Self {
        let u_res_rounded = round_up(u_res);
        let v_res_rounded = round_up(v_res);

        let u_blocks = u_res_rounded >> LOG_BLOCK_SIZE;
        let n_alloc = u_res_rounded * v_res_rounded;

        Self {
            data: vec![T::default(); n_alloc],
            u_blocks,
            u_res,
            v_res,
        }
    }

    /// Create a new `BlockedArray<T>` from slice.
    ///
    /// * `u_res`          - Size in u-dimension.
    /// * `v_res`          - Size in v-dimension.
    /// * `log_block_size` - Logarithm base 2 value used to ensure block sizes
    ///                      will be powers of 2.
    pub fn from_slice(u_res: usize, v_res: usize, data: &[T]) -> Self {
        let mut ret = Self::new(u_res, v_res);
        for v in 0..v_res {
            for u in 0..u_res {
                ret[(u, v)] = data[v * u_res + u];
            }
        }
        ret
    }

    /// Returns the size in u-dimension.
    pub fn u_size(&self) -> usize {
        self.u_res
    }

    /// Returns the size in v-dimension.
    pub fn v_size(&self) -> usize {
        self.v_res
    }

    /// Returns a linear `Vec<T>`.
    pub fn linear_vec(&self) -> Vec<T> {
        let mut a = Vec::with_capacity(self.u_res * self.v_res);
        let mut i = 0;
        for v in 0..self.v_res {
            for u in 0..self.u_res {
                a[i] = self[(u, v)];
                i += 1;
            }
        }
        a
    }

    /// Returns the offset into `data` for a given (x, y) index.
    ///
    /// * `index` -  A tuple containing (x, y).
    fn offset(&self, index: (usize, usize)) -> usize {
        let (u, v) = index;
        let bu = block(u);
        let bv = block(v);
        let ou = offset(u);
        let ov = offset(v);
        BLOCK_SIZE * BLOCK_SIZE * (self.u_blocks * bv + bu) + BLOCK_SIZE * ov + ou
    }
}

impl<T> Index<(usize, usize)> for BlockedArray<T>
where
    T: Copy + Default,
{
    type Output = T;

    /// Index the `BlockedArray<T>` by a tuple of (x, y) indexes.
    ///
    /// * `index` -  A tuple containing (x, y).
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let offset = self.offset(index);
        &self.data[offset]
    }
}

impl<T> IndexMut<(usize, usize)> for BlockedArray<T>
where
    T: Copy + Default,
{
    /// Index the `BlockedArray<T>` by a tuple of (x, y) indexes.
    ///
    /// * `index` -  A tuple containing (x, y).
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let offset = self.offset(index);
        &mut self.data[offset]
    }
}

/// Use log base 2 to ensure blocks are powers of 2.
const LOG_BLOCK_SIZE: usize = 2;

/// Block size.
const BLOCK_SIZE: usize = 1 << LOG_BLOCK_SIZE;

/// Rounds up to power of 2.
#[inline]
fn round_up(x: usize) -> usize {
    (x + BLOCK_SIZE - 1) & !(BLOCK_SIZE - 1)
}

/// Return the block for an index.
///
/// * `a` - The index.
#[inline]
pub fn block(a: usize) -> usize {
    a >> LOG_BLOCK_SIZE
}

/// Return the offset for an index.
///
/// * `a` - The index.
#[inline]
pub fn offset(a: usize) -> usize {
    a & (BLOCK_SIZE - 1)
}
