//! Random Number Generator.

#![allow(dead_code)]
use super::pbrt::Float;
use hexf::*;

/// 32-bit precision value for 1 - epsilon.
const FLOAT_ONE_MINUS_EPSILON: f32 = hexf32!("0x1.fffffep-1"); // 0.99999994

/// 64-bit precision value for 1 - epsilon.
const DOUBLE_ONE_MINUS_EPSILON: f64 = hexf64!("0x1.fffffffffffffp-1"); // 0.99999999999999989

/// 1 - epsilon in the precision we've selected for `Float`.
pub const ONE_MINUS_EPSILON: Float = FLOAT_ONE_MINUS_EPSILON;
