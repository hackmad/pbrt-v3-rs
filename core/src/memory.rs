//! Blocked Array

use std::ops::{Index, IndexMut};

/// Stores 2D arrays in a blocked memory layout.
///
/// * `T`              - Type of value to store.
/// * `LOG_BLOCK_SIZE` - Use log base to ensure blocks are powers of this base. Typically powers of 2 are used.
#[derive(Clone)]
pub struct BlockedArray<T, const LOG_BLOCK_SIZE: usize> {
    /// Array data.
    data: Vec<T>,

    /// Number of blocks in u-dimension.
    u_blocks: usize,

    /// Size in u-dimension.
    u_res: usize,

    /// Size in v-dimension.
    v_res: usize,

    /// Block size.
    block_size: usize,
}

impl<T, const LOG_BLOCK_SIZE: usize> BlockedArray<T, LOG_BLOCK_SIZE>
where
    T: Copy + Default,
{
    /// Create a new `BlockedArray<T>` with default values.
    ///
    /// * `u_res` - Size in u-dimension.
    /// * `v_res` - Size in v-dimension.
    pub fn new(u_res: usize, v_res: usize) -> Self {
        let block_size = 1 << LOG_BLOCK_SIZE;

        let u_res_rounded = round_up(u_res, block_size);
        let v_res_rounded = round_up(v_res, block_size);

        let u_blocks = u_res_rounded >> LOG_BLOCK_SIZE;
        let n_alloc = u_res_rounded * v_res_rounded;

        let data: Vec<T> = vec![T::default(); n_alloc];

        Self {
            data,
            block_size,
            u_blocks,
            u_res,
            v_res,
        }
    }

    /// Create a new `BlockedArray<T>` from slice.
    ///
    /// * `u_res` - Size in u-dimension.
    /// * `v_res` - Size in v-dimension.
    /// * `data`  - Slice containing the items.
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
    #[inline(always)]
    pub fn u_size(&self) -> usize {
        self.u_res
    }

    /// Returns the size in v-dimension.
    #[inline(always)]
    pub fn v_size(&self) -> usize {
        self.v_res
    }

    /// Returns a linear `Vec<T>`.
    pub fn linear_vec(&self) -> Vec<T> {
        let mut a = Vec::with_capacity(self.u_res * self.v_res);
        for v in 0..self.v_res {
            for u in 0..self.u_res {
                a.push(self[(u, v)]);
            }
        }
        a
    }

    /// Returns the offset into `data` for a given (x, y) coordinate.
    ///
    /// * `index` -  A tuple containing (x, y).
    #[inline(always)]
    fn coord_offset(&self, coord: (usize, usize)) -> usize {
        let (u, v) = coord;
        let bu = self.block(u);
        let bv = self.block(v);
        let ou = self.offset(u);
        let ov = self.offset(v);
        self.block_size * self.block_size * (self.u_blocks * bv + bu) + self.block_size * ov + ou
    }

    /// Return the block for an index.
    ///
    /// * `a` - The index.
    #[inline(always)]
    fn block(&self, a: usize) -> usize {
        a >> LOG_BLOCK_SIZE
    }

    /// Return the offset for an index.
    ///
    /// * `a` - The index.
    #[inline(always)]
    fn offset(&self, a: usize) -> usize {
        a & (self.block_size - 1)
    }
}

impl<T, const LOG_BLOCK_SIZE: usize> Index<(usize, usize)> for BlockedArray<T, LOG_BLOCK_SIZE>
where
    T: Copy + Default,
{
    type Output = T;

    /// Index the `BlockedArray<T>` by a tuple of (x, y) indexes.
    ///
    /// * `coord` -  A tuple containing (x, y).
    fn index(&self, coord: (usize, usize)) -> &Self::Output {
        let offset = self.coord_offset(coord);
        &self.data[offset]
    }
}

impl<T, const LOG_BLOCK_SIZE: usize> IndexMut<(usize, usize)> for BlockedArray<T, LOG_BLOCK_SIZE>
where
    T: Copy + Default,
{
    /// Index the `BlockedArray<T>` by a tuple of (x, y) indexes.
    ///
    /// * `coord` -  A tuple containing (x, y).
    fn index_mut(&mut self, coord: (usize, usize)) -> &mut Self::Output {
        let offset = self.coord_offset(coord);
        &mut self.data[offset]
    }
}

/// Rounds up to next power.
///
/// * `x`          - Value to round up.
/// * `block_size` - Block size based on LOG_BLOCK_SIZE.
#[inline(always)]
fn round_up(x: usize, block_size: usize) -> usize {
    (x + block_size - 1) & !(block_size - 1)
}
