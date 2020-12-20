//! Random Number Generator.

#![allow(dead_code)]
use super::pbrt::{min, Float};
use hexf::*;
use rand::distributions::uniform::SampleUniform;
use rand::distributions::{Distribution, Standard, Uniform};
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg32;

/// 32-bit precision value for 1 - epsilon.
pub const FLOAT_ONE_MINUS_EPSILON: f32 = hexf32!("0x1.fffffep-1"); // 0.99999994

/// 64-bit precision value for 1 - epsilon.
pub const DOUBLE_ONE_MINUS_EPSILON: f64 = hexf64!("0x1.fffffffffffffp-1"); // 0.99999999999999989

/// 1 - epsilon in the precision we've selected for `Float`.
pub const ONE_MINUS_EPSILON: Float = FLOAT_ONE_MINUS_EPSILON;

const PCG32_DEFAULT_STATE: u64 = 0x853c49e6748fea9b;
const PCG32_DEFAULT_STREAM: u64 = 0xda3e39cb94b95bdb;

/// Interface for generating uniform samples.
pub trait UniformRandom<T>
where
    Standard: Distribution<T>,
{
    /// Returns a uniformly distributed value.
    fn uniform(&mut self) -> T;
}

/// Implements the pseudo-random number generator.
#[derive(Clone)]
pub struct RNG {
    r: Pcg32,
}

impl Default for RNG {
    /// Return a new instance of `RNG` with default state and stream.
    fn default() -> Self {
        Self {
            r: Pcg32::new(PCG32_DEFAULT_STATE, PCG32_DEFAULT_STREAM),
        }
    }
}

impl RNG {
    /// Create a new `RNG` by seeding it with the given starting sequence.
    ///
    /// * `sequence_index` - The starting sequence to seed with.
    pub fn new(sequence_index: u64) -> Self {
        Self {
            r: Pcg32::seed_from_u64(sequence_index),
        }
    }

    /// Returns a uniformly distributed value over the closed interval containing
    /// the given bounds.
    ///
    /// * `lower_bound` - The upper bound.
    /// * `upper_bound` - The upper bound.
    pub fn bounded_uniform<T>(&mut self, lower_bound: T, upper_bound: T) -> T
    where
        T: SampleUniform,
    {
        let between = Uniform::from(lower_bound..upper_bound);
        between.sample(&mut self.r)
    }

    /// Randomly permute a slice containing n-dimensional values in a linear
    /// structure.
    ///
    /// * `v`            - The slice to shuffle.
    /// * `count`        - Number n-dimensional values.
    /// * `n_dimensions` - Number of total dimensions.
    pub fn shuffle<T>(&mut self, v: &mut [T], count: usize, n_dimensions: usize) {
        debug_assert!(count * n_dimensions == v.len());

        for i in 0..count {
            let other = i + self.bounded_uniform(0, count - i);
            for j in 0..n_dimensions {
                v.swap(n_dimensions * i + j, n_dimensions * other + j);
            }
        }
    }
}

/// Use default implementation for `UniformRandom` that wraps `Rng::gen<T>()`.
macro_rules! uniform_rand {
    ($t: ty) => {
        impl UniformRandom<$t> for RNG {
            /// Returns a uniformly distributed value over the closed interval
            /// containing the minimum value to maximum value of `$t`.
            fn uniform(&mut self) -> $t {
                let r = &mut self.r;
                r.gen::<$t>()
            }
        }
    };
}
uniform_rand!(u8);
uniform_rand!(u16);
uniform_rand!(u32);
uniform_rand!(u64);

impl UniformRandom<Float> for RNG {
    /// Returns a uniformly distributed value over the half open interval [0.0, 1.0).
    fn uniform(&mut self) -> Float {
        let r = &mut self.r;
        min(ONE_MINUS_EPSILON, r.gen::<Float>())
    }
}
