//! Rays

#![allow(dead_code)]
use super::{
    next_float_down, next_float_up, ArcMedium, Dot, Float, Normal3f, Point3f, Vector3, Vector3f,
    INFINITY,
};
use std::fmt::{Debug, Formatter, Result};

/// A Ray
#[derive(Clone)]
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

    /// Medium containing the origin.
    pub medium: Option<ArcMedium>,
}

impl Ray {
    /// Returns a ray with no differentials.
    ///
    /// * `o`      - Origin.
    /// * `d`      - Direction.
    /// * `t_max`  - Maximum extent of the ray.
    /// * `time`   - Time value.
    /// * `medium` - Medium containing origin `o`.
    pub fn new(
        o: Point3f,
        d: Vector3f,
        t_max: Float,
        time: Float,
        medium: Option<ArcMedium>,
    ) -> Self {
        Self {
            o,
            d,
            t_max,
            time,
            differentials: None::<RayDifferential>,
            medium,
        }
    }

    /// Returns a ray with differential.
    ///
    /// * `o`             - Origin.
    /// * `d`             - Direction.
    /// * `t_max`         - Maximum extent of the ray.
    /// * `time`          - Time value.
    /// * `differentials` - Auxilliary rays offset by one sample in x and y direction.
    /// * `medium` - Medium containing origin `o`.
    pub fn new_with_differentials(
        o: Point3f,
        d: Vector3f,
        t_max: Float,
        time: Float,
        differentials: RayDifferential,
        medium: Option<ArcMedium>,
    ) -> Self {
        Self {
            o,
            d,
            t_max,
            time,
            differentials: Some(differentials),
            medium,
        }
    }

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

    /// Offset ray origin along the surface normal.
    ///
    /// `p`       - Intersection point.
    /// `p_error` - Floating point error for intersection points.
    /// `w`       - The direction.
    /// `n`       - Surface normal at the point `p`.
    pub fn offset_origin(p: &Point3f, p_error: &Vector3f, n: &Normal3f, w: &Vector3f) -> Point3f {
        let d = n.abs().dot(p_error);

        let mut offset = d * Vector3::from(*n);
        if w.dot(n) < 0.0 {
            offset = -offset;
        }

        let mut po = *p + offset;

        // Round offset point po away from p.
        for axis in 0..3 {
            if offset[axis] > 0.0 {
                po[axis] = next_float_up(po[axis]);
            } else if offset[axis] < 0.0 {
                po[axis] = next_float_down(po[axis]);
            }
        }

        po
    }
}

impl Debug for Ray {
    /// Display the ray parameters.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.debug_struct("Ray")
            .field("o", &self.o)
            .field("d", &self.d)
            .field("t_max", &self.t_max)
            .field("time", &self.time)
            .field("differentials", &self.differentials)
            .finish()
    }
}

impl Default for Ray {
    /// Returns a default value for `Ray`.
    fn default() -> Self {
        Self {
            o: Point3f::default(),
            d: Vector3f::default(),
            t_max: INFINITY,
            time: 0.0,
            differentials: None,
            medium: None,
        }
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

impl RayDifferential {
    /// Returns a ray differential.
    ///
    /// * `xo` - Origin for x-direction differential.
    /// * `xd` - Direction for x-direction differential.
    /// * `yo` - Origin for y-direction differential.
    /// * `yd` - Direction for y-direction differential.
    pub fn new(xo: Point3f, yo: Point3f, xd: Vector3f, yd: Vector3f) -> Self {
        Self {
            rx_origin: xo,
            ry_origin: yo,
            rx_direction: xd,
            ry_direction: yd,
        }
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
#[macro_use]
mod tests {
    use super::super::{Point3, Vector3, INFINITY};
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn has_nans() {
        let nan_point = Point3::new(f32::NAN, f32::NAN, f32::NAN);
        let nan_vector = Vector3::new(f32::NAN, f32::NAN, f32::NAN);
        let nan_t_max = f32::NAN;
        let point = Point3::new(0.0, 0.0, 0.0);
        let vector = Vector3::new(1.0, 0.0, 0.0);

        assert!(Ray::new(nan_point, vector, 0.0, 0.0, None).has_nans());
        assert!(Ray::new(point, nan_vector, 0.0, 0.0, None).has_nans());
        assert!(Ray::new(point, vector, nan_t_max, 0.0, None).has_nans());
        assert!(Ray::new(nan_point, nan_vector, nan_t_max, 0.0, None).has_nans());
        assert!(!Ray::new(point, vector, 0.0, 0.0, None).has_nans());
    }

    #[test]
    fn at() {
        let o = Point3::new(0.0, 0.0, 0.0);
        let d = Vector3::new(1.0, 1.0, 1.0);
        let r = Ray::new(o, d, INFINITY, 0.0, None);
        assert!(r.at(0.0) == o);
        assert!(r.at(1.0) == Point3::from(d));
    }

    #[test]
    fn scale_differentials_none() {
        let o = Point3::new(0.0, 0.0, 0.0);
        let d = Vector3::new(1.0, 1.0, 1.0);

        let mut r = Ray::new(o, d, INFINITY, 0.0, None);
        assert!(r.differentials.is_none());

        r.scale_differentials(2.0);
        assert!(r.differentials.is_none());
    }

    #[test]
    fn scale_differentials_some() {
        let o = Point3::new(0.0, 0.0, 0.0);
        let d = Vector3::new(1.0, 1.0, 1.0);
        let xo = Point3::new(1.0, 0.0, 0.0);
        let yo = Point3::new(0.0, 1.0, 0.0);
        let xd = Vector3::new(1.0, 0.0, 0.0);
        let yd = Vector3::new(0.0, 1.0, 0.0);

        let rd = RayDifferential::new(xo, yo, xd, yd);

        let mut r = Ray::new_with_differentials(o, d, INFINITY, 0.0, rd, None);
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
            let r = Ray::new(o, d, INFINITY, 0.0, None);
            prop_assert_eq!(r.at(t), o + t * d);
        }
    }
}
