//! Axis

use std::ops::Add;

/// Axis enumeration
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Axis {
    X = 0,
    Y = 1,
    Z = 2,
}
impl From<usize> for Axis {
    fn from(i: usize) -> Self {
        match i {
            0 => Axis::X,
            1 => Axis::Y,
            2 => Axis::Z,
            _ => panic!("invalid axis value"),
        }
    }
}
impl From<u8> for Axis {
    fn from(i: u8) -> Self {
        match i {
            0 => Axis::X,
            1 => Axis::Y,
            2 => Axis::Z,
            _ => panic!("invalid axis value"),
        }
    }
}
impl From<Axis> for u8 {
    fn from(axis: Axis) -> Self {
        match axis {
            Axis::X => 0_u8,
            Axis::Y => 1_u8,
            Axis::Z => 2_u8,
        }
    }
}
impl From<Axis> for usize {
    fn from(axis: Axis) -> usize {
        match axis {
            Axis::X => 0_usize,
            Axis::Y => 1_usize,
            Axis::Z => 2_usize,
        }
    }
}
impl Add<usize> for Axis {
    type Output = Axis;
    fn add(self, i: usize) -> Self::Output {
        Axis::from((self as usize + i) % 3)
    }
}
impl Default for Axis {
    fn default() -> Self {
        Axis::X
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
use proptest::prelude::*;

#[cfg(test)]
pub fn axis_2d_strategy() -> impl Strategy<Value = Axis> {
    prop_oneof![Just(Axis::X), Just(Axis::Y)]
}

#[cfg(test)]
pub fn axis_3d_strategy() -> impl Strategy<Value = Axis> {
    prop_oneof![Just(Axis::X), Just(Axis::Y), Just(Axis::Z)]
}
