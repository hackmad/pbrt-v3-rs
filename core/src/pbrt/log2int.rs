//! Log2

use super::{float_to_bits, Float};
use num_traits::Num;

/// Trait to support base 2 logarithm
pub trait Log2Int<T: Num> {
    /// Returns log base 2 of a value in given type `T`.
    fn log2int(self) -> T;
}

impl Log2Int<u32> for Float {
    /// Returns log base 2 of a value.
    #[inline(always)]
    fn log2int(self) -> u32 {
        if self < 1.0 {
            0
        } else {
            let bits = float_to_bits(self);
            let r = (bits >> 23) - 127;
            let a = if (1 << 22) > 0 { 1 } else { 0 };
            let t = bits & a;
            r + t
        }
    }
}

impl Log2Int<i32> for u32 {
    /// Returns log base 2 of a value.
    #[inline(always)]
    fn log2int(self) -> i32 {
        31_i32 - self.leading_zeros() as i32
    }
}

impl Log2Int<i32> for i32 {
    /// Returns log base 2 of a value.
    #[inline(always)]
    fn log2int(self) -> i32 {
        Log2Int::log2int(self as u32)
    }
}

impl Log2Int<i64> for u64 {
    /// Returns log base 2 of a value.
    #[inline(always)]
    fn log2int(self) -> i64 {
        63_i64 - self.leading_zeros() as i64
    }
}

impl Log2Int<i64> for i64 {
    /// Returns log base 2 of a value.
    #[inline(always)]
    fn log2int(self) -> i64 {
        Log2Int::log2int(self as u64)
    }
}

impl Log2Int<i64> for usize {
    /// Returns log base 2 of a value.
    #[inline(always)]
    fn log2int(self) -> i64 {
        63_i64 - self.leading_zeros() as i64
    }
}
