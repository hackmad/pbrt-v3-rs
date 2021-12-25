//! Random Number Generator.

use crate::pbrt::*;

/// 32-bit precision value for 1 - epsilon.
pub const FLOAT_ONE_MINUS_EPSILON: f32 = hexf32!("0x1.fffffep-1"); // 0.99999994

/// 64-bit precision value for 1 - epsilon.
pub const DOUBLE_ONE_MINUS_EPSILON: f64 = hexf64!("0x1.fffffffffffffp-1"); // 0.99999999999999989

/// 1 - epsilon in the precision we've selected for `Float`.
pub const ONE_MINUS_EPSILON: Float = FLOAT_ONE_MINUS_EPSILON;

const PCG32_DEFAULT_STATE: u64 = 0x853c49e6748fea9b;
const PCG32_DEFAULT_STREAM: u64 = 0xda3e39cb94b95bdb;
const PCG32_MULT: u64 = 0x5851f42d4c957f2d;

/// Implements the pseudo-random number generator.
#[derive(Clone)]
pub struct RNG {
    state: u64,
    inc: u64,
}

impl Default for RNG {
    /// Return a new instance of `RNG` with default state and stream.
    fn default() -> Self {
        Self {
            state: PCG32_DEFAULT_STATE,
            inc: PCG32_DEFAULT_STREAM,
        }
    }
}

impl RNG {
    /// Create a new `RNG` by seeding it with the given starting sequence.
    ///
    /// * `sequence_index` - The starting sequence to seed with.
    pub fn new(sequence_index: u64) -> Self {
        let mut ret = Self { state: 0, inc: 0 };
        ret.set_sequence(sequence_index);
        ret
    }

    /// Initialize the random number generator sequence.
    ///
    /// * `init_seq` - The starting sequence to seed with.
    #[inline(always)]
    fn set_sequence(&mut self, init_seq: u64) {
        self.state = 0;
        let (inc, _) = init_seq.overflowing_shl(1);
        self.inc = inc | 1;
        let _ = self.uniform_u32();

        let (state, _) = self.state.overflowing_add(PCG32_DEFAULT_STATE);
        self.state = state;
        let _ = self.uniform_u32();
    }

    /// Returns a uniformly distributed u32 value.
    #[inline(always)]
    pub fn uniform_u32(&mut self) -> u32 {
        let old_state = self.state;
        let (new_state, _) = old_state.overflowing_mul(PCG32_MULT);
        let (new_state, _) = new_state.overflowing_add(self.inc);
        self.state = new_state;

        let (xor_shifted, _) = old_state.overflowing_shr(18);
        let (xor_shifted, _) = (xor_shifted ^ old_state).overflowing_shr(27);
        let xor_shifted = xor_shifted as u32;

        let (rot, _) = old_state.overflowing_shr(59);
        let rot = rot as u32;

        let (r1, _) = xor_shifted.overflowing_shr(rot);
        let (bits, _) = (!rot).overflowing_add(1);
        let (r2, _) = xor_shifted.overflowing_shl(bits & 31);

        r1 | r2
    }

    /// Returns a uniformly distributed value over the closed interval containing
    /// the given bounds.
    ///
    /// * `lower_bound` - The upper bound.
    /// * `upper_bound` - The upper bound.
    pub fn bounded_uniform_u32(&mut self, lower_bound: u32, upper_bound: u32) -> u32 {
        let b = upper_bound - lower_bound;
        let threshold = (!b + 1) % b;
        loop {
            let r = self.uniform_u32();
            if r >= threshold {
                return lower_bound + r % b;
            }
        }
    }

    /// Returns a uniformly distributed value over the half open interval [0.0, 1.0).
    pub fn uniform_float(&mut self) -> Float {
        min(
            self.uniform_u32() as Float * hexf32!("0x1.0p-32") as Float,
            FLOAT_ONE_MINUS_EPSILON,
        )
    }

    /// Randomly permute a slice containing n-dimensional values in a linear
    /// structure.
    ///
    /// * `v`            - The slice to shuffle.
    /// * `count`        - Number n-dimensional values.
    /// * `n_dimensions` - Number of total dimensions.
    pub fn shuffle<T>(&mut self, v: &mut [T], count: usize, n_dimensions: usize) {
        debug_assert!(count * n_dimensions <= v.len());

        for i in 0..count {
            let other = i + self.bounded_uniform_u32(0, (count - i) as u32) as usize;
            for j in 0..n_dimensions {
                v.swap(n_dimensions * i + j, n_dimensions * other + j);
            }
        }
    }
}
