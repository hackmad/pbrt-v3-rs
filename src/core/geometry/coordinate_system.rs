//! 3-D Coordinate System

#![allow(dead_code)]
use crate::core::geometry::Vector3;
use crate::core::pbrt::abs;
use num_traits::Float;

/// Create a new coordinate system from a single unit vector `v1` and return
/// the new vectors:
/// * `v2` constructed from the first by zeroing one of the coordinates and
/// swapping the remaining 2 and negating one of them. This vector is normalized.
/// * `v3` Cross product of `v1 and `v2`. Since both these are normalized, this
/// will be a unit vector.
///
/// * `v1` - The first unit vector to form part of the coordinate system.
pub fn coordinate_system<T: Float>(v1: &Vector3<T>) -> (Vector3<T>, Vector3<T>) {
    let v2 = if abs(v1.x) > abs(v1.y) {
        Vector3::new(-v1.z, T::zero(), v1.x) / (v1.x * v1.x + v1.z * v1.z).sqrt()
    } else {
        Vector3::new(T::zero(), v1.z, -v1.y) / (v1.y * v1.y + v1.z * v1.z).sqrt()
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
        let v1 = Vector3::new(1.0, 0.0, 0.0);
        let (v2, v3) = coordinate_system(&v1);
        assert!(v2 == Vector3::new(0.0, 0.0, 1.0));
        assert!(v3 == Vector3::new(0.0, -1.0, 0.0));
    }

    #[test]
    fn from_x_axis() {
        let v1 = Vector3::new(2.0, 0.0, 0.0);
        let (v2, v3) = coordinate_system(&v1);
        assert!(v2 == Vector3::new(0.0, 0.0, 1.0));
        assert!(v3 == Vector3::new(0.0, -2.0, 0.0));
    }

    #[test]
    fn from_vector_x_greater_than_y() {
        let v1 = Vector3::new(0.5, 0.2, 0.5);
        let (v2, v3) = coordinate_system(&v1);
        assert!(v1.dot(&v2) == 0.0);
        assert!(v1.dot(&v3) == 0.0);
        assert!(v2.dot(&v3) == 0.0);
    }

    #[test]
    fn from_vector_x_less_than_y() {
        let v1 = Vector3::new(0.2, 0.5, 0.5);
        let (v2, v3) = coordinate_system(&v1);
        assert!(v1.dot(&v2) == 0.0);
        assert!(v1.dot(&v3) == 0.0);
        assert!(v2.dot(&v3) == 0.0);
    }
}
