//! Rays

#![allow(dead_code)]
use super::{Float, Point3f, Vector3f};

/// A Ray
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Ray {
    /// Origin.
    pub o: Point3f,

    /// Direction.
    pub d: Vector3f,

    /// Maximum extent of the ray.
    pub t_max: Float,

    /// Time value.
    pub time: Float,

    /// Auxilliary rays offset by one sample in x and y direction.
    pub differentials: Option<RayDifferential>,
}

impl Ray {
    /// Returns true if either coordinate is NaN.
    pub fn has_nans(&self) -> bool {
        self.o.has_nans() || self.d.has_nans() || self.t_max.is_nan()
    }

    /// Get position along the ray at given parameter.
    ///
    /// * `t` - Parameter to evaluate.
    pub fn at(&self, t: Float) -> Point3f {
        self.o + self.d * t
    }

    /// Scale the differential rays to account for spacing between samples on
    /// the film plane.
    ///
    /// * `s` - The weight used to scale the differential rays.
    pub fn scale_differentials(&mut self, s: Float) {
        if let Some(d) = self.differentials {
            self.differentials = Some(RayDifferential {
                rx_origin: self.o + (d.rx_origin - self.o) * s,
                ry_origin: self.o + (d.ry_origin - self.o) * s,
                rx_direction: self.d + (d.rx_direction - self.d) * s,
                ry_direction: self.d + (d.ry_direction - self.d) * s,
            });
        }
    }
}

/// Returns a ray with no differentials.
///
/// * `o`     - Origin.
/// * `d`     - Direction.
/// * `t_max` - Maximum extent of the ray.
/// * `time`  - Time value.
pub fn ray(o: Point3f, d: Vector3f, t_max: Float, time: Float) -> Ray {
    Ray {
        o,
        d,
        t_max,
        time,
        differentials: None::<RayDifferential>,
    }
}

/// Returns a ray with differential.
///
/// * `o`             - Origin.
/// * `d`             - Direction.
/// * `t_max`         - Maximum extent of the ray.
/// * `time`          - Time value.
/// * `differentials` - Auxilliary rays offset by one sample in x and y direction.
pub fn ray_with_differentials(
    o: Point3f,
    d: Vector3f,
    t_max: Float,
    time: Float,
    differentials: RayDifferential,
) -> Ray {
    Ray {
        o,
        d,
        t_max,
        time,
        differentials: Some(differentials),
    }
}

/// A ray differential is offset by one sample in the x and y direction of a ray.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct RayDifferential {
    /// Origin of ray offset in x-direction.
    pub rx_origin: Point3f,

    /// Origin of ray offset in y-direction.
    pub ry_origin: Point3f,

    /// Direction of ray offset in x-direction.
    pub rx_direction: Vector3f,

    /// Direction of ray offset in y-direction.
    pub ry_direction: Vector3f,
}

/// Returns a ray differential.
///
/// * `xo` - Origin for x-direction differential.
/// * `xd` - Direction for x-direction differential.
/// * `yo` - Origin for y-direction differential.
/// * `yd` - Direction for y-direction differential.
pub fn ray_differential(xo: Point3f, yo: Point3f, xd: Vector3f, yd: Vector3f) -> RayDifferential {
    RayDifferential {
        rx_origin: xo,
        ry_origin: yo,
        rx_direction: xd,
        ry_direction: yd,
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
#[macro_use]
mod tests {
    use super::super::{point3, vector3, Point3, Vector3, INFINITY};
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn has_nans() {
        let nan_point = point3(f32::NAN, f32::NAN, f32::NAN);
        let nan_vector = vector3(f32::NAN, f32::NAN, f32::NAN);
        let nan_t_max = f32::NAN;
        let point = point3(0.0, 0.0, 0.0);
        let vector = vector3(1.0, 0.0, 0.0);

        assert!(ray(nan_point, vector, 0.0, 0.0).has_nans());
        assert!(ray(point, nan_vector, 0.0, 0.0).has_nans());
        assert!(ray(point, vector, nan_t_max, 0.0).has_nans());
        assert!(ray(nan_point, nan_vector, nan_t_max, 0.0).has_nans());
        assert!(!ray(point, vector, 0.0, 0.0).has_nans());
    }

    #[test]
    fn at() {
        let o = point3(0.0, 0.0, 0.0);
        let d = vector3(1.0, 1.0, 1.0);
        let r = ray(o, d, INFINITY, 0.0);
        assert!(r.at(0.0) == o);
        assert!(r.at(1.0) == Point3::from(d));
    }

    #[test]
    fn scale_differentials_none() {
        let o = point3(0.0, 0.0, 0.0);
        let d = vector3(1.0, 1.0, 1.0);

        let mut r = ray(o, d, INFINITY, 0.0);
        assert!(r.differentials.is_none());

        r.scale_differentials(2.0);
        assert!(r.differentials.is_none());
    }

    #[test]
    fn scale_differentials_some() {
        let o = point3(0.0, 0.0, 0.0);
        let d = vector3(1.0, 1.0, 1.0);
        let xo = point3(1.0, 0.0, 0.0);
        let yo = point3(0.0, 1.0, 0.0);
        let xd = vector3(1.0, 0.0, 0.0);
        let yd = vector3(0.0, 1.0, 0.0);

        let rd = ray_differential(xo, yo, xd, yd);

        let mut r = ray_with_differentials(o, d, INFINITY, 0.0, rd);
        assert!(!r.differentials.is_none());

        r.scale_differentials(2.0);
        assert!(!r.differentials.is_none());

        let scaled_differentials = r.differentials.unwrap();
        assert!(scaled_differentials.rx_origin == o + 2.0 * (xo - o));
        assert!(scaled_differentials.ry_origin == o + 2.0 * (yo - o));
        assert!(scaled_differentials.rx_direction == d + 2.0 * (xd - d));
        assert!(scaled_differentials.ry_direction == d + 2.0 * (yd - d));
    }

    // Define some properties for tests.
    prop_range!(range_f32, f32, -100.0..100.0f32);

    prop_point3!(
        point3_f32,
        f32,
        -100.0..100.0f32,
        -100.0..100.0f32,
        -100.0..100.0f32
    );

    prop_vector3!(
        vector3_f32,
        f32,
        -100.0..100.0f32,
        -100.0..100.0f32,
        -100.0..100.0f32
    );

    proptest! {
        #[test]
        fn at_f32(o in point3_f32(), d in vector3_f32(), t in range_f32()) {
            let r = ray(o, d, INFINITY, 0.0);
            prop_assert_eq!(r.at(t), o + t * d);
        }
    }
}
