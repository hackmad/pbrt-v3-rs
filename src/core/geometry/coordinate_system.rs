//! 3-D Coordinate System

#![allow(dead_code)]
use crate::core::geometry::Vector3;
use crate::core::pbrt::abs;
use num_traits::Float;

/// Create a new coordinate system from a single unit vector and return
/// the new vectors.
///
/// * `v1` - The first unit vector to form part of the coordinate system.
/// * `v2` - Constructed from the first by zeroing one of the coordinates and
///          swapping the remaining 2 and negating one of them. This vector is
///          normalized.
/// * `v3` - Cross product of `v1 and `v2`. Since both these are normalized,
///          this will be a unit vector.
pub fn coordinate_system<T: Float>(v1: &Vector3<T>, v2: &mut Vector3<T>, v3: &mut Vector3<T>) {
    if abs(v1.x) > abs(v1.y) {
        *v2 = Vector3::new(-v1.z, T::zero(), v1.x) / (v1.x * v1.x + v1.z * v1.z).sqrt();
    } else {
        *v2 = Vector3::new(T::zero(), v1.z, -v1.y) / (v1.y * v1.y + v1.z * v1.z).sqrt();
    };
    *v3 = v1.cross(v2);
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
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        let mut v2 = Vector3::new(0.0, 0.0, 0.0);
        let mut v3 = Vector3::new(0.0, 0.0, 0.0);
        coordinate_system(&v1, &mut v2, &mut v3);
        assert!(v2 == Vector3::new(0.0, 0.0, 1.0));
        assert!(v3 == Vector3::new(0.0, -1.0, 0.0));
    }

    #[test]
    fn from_x_axis() {
        let v1 = Vector3::new(2.0, 0.0, 0.0);
        let mut v2 = Vector3::new(0.0, 0.0, 0.0);
        let mut v3 = Vector3::new(0.0, 0.0, 0.0);
        coordinate_system(&v1, &mut v2, &mut v3);
        assert!(v2 == Vector3::new(0.0, 0.0, 1.0));
        assert!(v3 == Vector3::new(0.0, -2.0, 0.0));
    }

    #[test]
    fn from_vector_x_greater_than_y() {
        let v1 = Vector3::new(0.5, 0.2, 0.5);
        let mut v2 = Vector3::new(0.0, 0.0, 0.0);
        let mut v3 = Vector3::new(0.0, 0.0, 0.0);
        coordinate_system(&v1, &mut v2, &mut v3);
        assert!(v1.dot(&v2) == 0.0);
        assert!(v1.dot(&v3) == 0.0);
        assert!(v2.dot(&v3) == 0.0);
    }

    #[test]
    fn from_vector_x_less_than_y() {
        let v1 = Vector3::new(0.2, 0.5, 0.5);
        let mut v2 = Vector3::new(0.0, 0.0, 0.0);
        let mut v3 = Vector3::new(0.0, 0.0, 0.0);
        coordinate_system(&v1, &mut v2, &mut v3);
        assert!(v1.dot(&v2) == 0.0);
        assert!(v1.dot(&v3) == 0.0);
        assert!(v2.dot(&v3) == 0.0);
    }
}
