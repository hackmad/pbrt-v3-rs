//! 3-D Coordinate System

#![allow(dead_code)]
use super::abs;
use super::vector3::{vector3, Vector3};
use num_traits::Float;

/// A coordinate system containing 3 orthogonal vectors.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct CoordinateSystem<T> {
    /// The first unit vector.
    pub v1: Vector3<T>,

    /// The second unit vector.
    pub v2: Vector3<T>,

    /// The third unit vector.
    pub v3: Vector3<T>,
}

impl<T: Float> From<Vector3<T>> for CoordinateSystem<T> {
    /// Create a new coordinate system from a single unit vector.
    ///
    /// A second vector is constructing from the first by zeroing one of the
    /// coordinates and swapping the remaining 2 and negating one of them. This
    /// vector is also normalized.
    ///
    /// The third vector is the cross product of the give vector and the second
    /// vector. Since both these are normalized, the third vector will be a unit
    /// vector.
    ///
    /// * `v1` - The first unit vector to form part of the coordinate system.
    fn from(v1: Vector3<T>) -> Self {
        let (v2, v3) = get_coordinate_system_vectors(&v1);
        CoordinateSystem { v1, v2, v3 }
    }
}

/// Create a new coordinate system from a single unit vector and return
/// the new vectors.
///
/// A second vector is constructing from the first by zeroing one of the
/// coordinates and swapping the remaining 2 and negating one of them. This
/// vector is also normalized.
///
/// The third vector is the cross product of the give vector and the second
/// vector. Since both these are normalized, the third vector will be a unit
/// vector.
///
/// * `v1` - The first unit vector to form part of the coordinate system.
pub fn get_coordinate_system_vectors<T: Float>(v1: &Vector3<T>) -> (Vector3<T>, Vector3<T>) {
    let v2 = if abs(v1.x) > abs(v1.y) {
        vector3(-v1.z, T::zero(), v1.x) / (v1.x * v1.x + v1.z * v1.z).sqrt()
    } else {
        vector3(T::zero(), v1.z, -v1.y) / (v1.y * v1.y + v1.z * v1.z).sqrt()
    };

    let v3 = v1.cross(&v2);

    (v2, v3)
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
#[macro_use]
mod tests {
    use super::super::common::Dot;
    use super::*;

    #[test]
    fn from_unit_x_axis() {
        let v1 = vector3(1.0, 0.0, 0.0);
        let cs = CoordinateSystem::from(v1);
        assert!(cs.v1 == v1);
        assert!(cs.v2 == vector3(0.0, 0.0, 1.0));
        assert!(cs.v3 == vector3(0.0, -1.0, 0.0));
    }

    #[test]
    fn from_x_axis() {
        let v1 = vector3(2.0, 0.0, 0.0);
        let cs = CoordinateSystem::from(v1);
        assert!(cs.v1 == v1);
        assert!(cs.v2 == vector3(0.0, 0.0, 1.0));
        assert!(cs.v3 == vector3(0.0, -2.0, 0.0));
    }

    #[test]
    fn from_vector_x_greater_than_y() {
        let v1 = vector3(0.5, 0.2, 0.5);
        let cs = CoordinateSystem::from(v1);
        assert!(cs.v1 == v1);
        assert!(cs.v1.dot(&cs.v2) == 0.0);
        assert!(cs.v1.dot(&cs.v3) == 0.0);
        assert!(cs.v2.dot(&cs.v3) == 0.0);
    }

    #[test]
    fn from_vector_x_less_than_y() {
        let v1 = vector3(0.2, 0.5, 0.5);
        let cs = CoordinateSystem::from(v1);
        assert!(cs.v1 == v1);
        assert!(cs.v1.dot(&cs.v2) == 0.0);
        assert!(cs.v1.dot(&cs.v3) == 0.0);
        assert!(cs.v2.dot(&cs.v3) == 0.0);
    }
}
