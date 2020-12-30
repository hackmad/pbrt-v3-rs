//! 3-D Axis Aligned Bounding Boxes.

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use num_traits::bounds::Bounded;
use num_traits::{Num, Zero};
use std::fmt;
use std::mem::swap;
use std::ops::{DivAssign, Index, Mul};

/// 3-D Axis Aligned Bounding Box.
#[derive(Copy, Clone, Default, PartialEq)]
pub struct Bounds3<T: Num> {
    /// Minimum bounds.
    pub p_min: Point3<T>,

    /// Maximum bounds.
    pub p_max: Point3<T>,
}

/// 3-D bounding box containing `Float` points.
pub type Bounds3f = Bounds3<Float>;

/// 3-D bounding box containing `Int` points.
pub type Bounds3i = Bounds3<Int>;

impl<T: Num + PartialOrd + Copy> From<Point3<T>> for Bounds3<T> {
    /// Use a 3-D point as minimum and maximum 3-D bounds.
    ///
    /// * `p` - 3-D point.
    fn from(p: Point3<T>) -> Self {
        Bounds3 { p_min: p, p_max: p }
    }
}

impl From<Bounds3i> for Bounds3f {
    /// Convert a `Bounds3i` to `Bounds3f`.
    ///
    /// * `b` - The `Bounds3i` to convert.
    fn from(b: Bounds3i) -> Self {
        Self {
            p_min: b.p_min.into(),
            p_max: b.p_max.into(),
        }
    }
}

impl From<Bounds3f> for Bounds3i {
    /// Convert a `Bounds3f` to `Bounds3i`.
    ///
    /// * `b` - The `Bounds3f` to convert.
    fn from(b: Bounds3f) -> Self {
        Self {
            p_min: b.p_min.into(),
            p_max: b.p_max.into(),
        }
    }
}

impl<T: Num + Copy> Bounds3<T> {
    /// Creates a new 3-D bounding box from 2 points. The minimum and maximum bounds
    /// are used for each coordinate axis.
    ///
    /// * `p1` - First point.
    /// * `p2` - Second point.
    pub fn new(p1: Point3<T>, p2: Point3<T>) -> Self
    where
        T: PartialOrd + Copy,
    {
        Self {
            p_min: Point3::new(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z)),
            p_max: Point3::new(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z)),
        }
    }

    /// Returns a 3-D bounding box where minimum and maximum bounds are maximum and
    /// minimum values respectively of the type's limits. This is so we can easily grow
    /// the bounding box from nothing iteratively.
    pub fn empty() -> Self
    where
        T: Bounded + PartialOrd + Copy,
    {
        // Don't call new() because it'll create the largest bounding box
        // by flipping p_min and p_max.
        Self {
            p_min: Point3::new(T::max_value(), T::max_value(), T::max_value()),
            p_max: Point3::new(T::min_value(), T::min_value(), T::min_value()),
        }
    }

    /// Returns true if the bounds describes an empty box where any the
    /// components of any p_max are less than p_max.
    pub fn is_empty(&self) -> bool
    where
        T: PartialOrd,
    {
        self.p_max.x < self.p_min.x || self.p_max.y < self.p_min.y || self.p_max.z < self.p_min.z
    }

    /// Returns the vector along the box diagonal from the minimum point to
    /// the maximum point.
    pub fn diagonal(&self) -> Vector3<T> {
        self.p_max - self.p_min
    }

    /// Returns the surface area of the bounding box.
    pub fn surface_area(&self) -> T
    where
        T: PartialOrd,
    {
        if self.is_empty() {
            T::zero()
        } else {
            let d = self.diagonal();
            let h = d.x * d.y + d.x * d.z + d.y * d.z;
            h + h
        }
    }

    /// Returns the volume of the bounding box.
    pub fn volume(&self) -> T
    where
        T: PartialOrd,
    {
        if self.is_empty() {
            T::zero()
        } else {
            let d = self.diagonal();
            d.x * d.y + d.z
        }
    }

    /// Returns the index of which of the axes is longest. This is useful, for
    /// example, when deciding which axis to subdivide when building some of
    /// the ray-intersection acceleration structures.
    pub fn maximum_extent(&self) -> Axis
    where
        T: PartialOrd,
    {
        let d = self.diagonal();
        if d.x > d.y && d.x > d.z {
            Axis::X
        } else if d.y > d.z {
            Axis::Y
        } else {
            Axis::Z
        }
    }

    /// Returns true if extents of another bounding box overlap with this one.
    ///
    /// * `other` - The other bounding box.
    pub fn overlaps(&self, other: &Self) -> bool
    where
        T: PartialOrd,
    {
        let x = (self.p_max.x >= other.p_min.x) && (self.p_min.x <= other.p_max.x);
        let y = (self.p_max.y >= other.p_min.y) && (self.p_min.y <= other.p_max.y);
        let z = (self.p_max.z >= other.p_min.z) && (self.p_min.z <= other.p_max.z);
        x && y && z
    }

    /// Returns the continuous position of a point relative to the corners of the
    /// box, where a point at the minimum corner has offset `(0, 0, 0)` and a
    /// point at the maximum corner has offset is `(1, 1, 1)`.
    ///
    /// * `p` - The point.
    pub fn offset(&self, p: &Point3<T>) -> Vector3<T>
    where
        T: num_traits::Float + DivAssign<T> + PartialOrd + Copy,
    {
        let mut o = *p - self.p_min;
        if self.p_max.x > self.p_min.x {
            o.x /= self.p_max.x - self.p_min.x;
        }
        if self.p_max.y > self.p_min.y {
            o.y /= self.p_max.y - self.p_min.y;
        }
        if self.p_max.z > self.p_min.z {
            o.z /= self.p_max.z - self.p_min.z;
        }
        o
    }

    /// Returns true if a point is inside the bounding box.
    ///
    /// * `p` - The point.
    pub fn contains(&self, p: &Point3<T>) -> bool
    where
        T: PartialOrd,
    {
        (p.x >= self.p_min.x && p.x <= self.p_max.x)
            && (p.y >= self.p_min.y && p.y <= self.p_max.y)
            && (p.z >= self.p_min.z && p.z <= self.p_max.z)
    }

    /// Returns true if a point is inside the bounding box. The upper boundary
    /// is considered out of bounds. This is useful for integer-typed bounds.
    ///
    /// * `p` - The point.
    pub fn contains_exclusive(&self, p: &Point3<T>) -> bool
    where
        T: PartialOrd,
    {
        (p.x >= self.p_min.x && p.x < self.p_max.x)
            && (p.y >= self.p_min.y && p.y < self.p_max.y)
            && (p.z >= self.p_min.z && p.z < self.p_max.z)
    }

    /// Return the center and radius of a sphere bounded on the corners of the
    /// bounding box.
    pub fn bounding_sphere(&self) -> (Point3<T>, T)
    where
        T: num_traits::Float + Zero,
        Float: Mul<Point3<T>, Output = Point3<T>>,
    {
        let center = lerp(0.5, self.p_min, self.p_max);
        let radius = if self.contains(&center) {
            center.distance(self.p_max)
        } else {
            T::zero()
        };
        (center, radius)
    }

    /// Linearly interpolates between the corners of the box by the given amount
    /// in each dimension.
    ///
    /// * `t` - The interpolation parameter in x and y directions.
    pub fn lerp(&self, t: &Point3f) -> Point3<T>
    where
        Float: Mul<T, Output = T>,
    {
        Point3::new(
            lerp::<T>(t.x, self.p_min.x, self.p_max.x),
            lerp::<T>(t.y, self.p_min.y, self.p_max.y),
            lerp::<T>(t.z, self.p_min.z, self.p_max.z),
        )
    }

    /// Pad the bounding box by a constant factor in all dimensions.
    ///
    /// * `delta` - Padding amount.
    pub fn expand(&self, delta: T) -> Self
    where
        T: PartialOrd,
    {
        // Don't call bounds2<T>() to prevent flipping p_min and p_max when
        // the input is empty box.
        Self {
            p_min: self.p_min - Vector3::new(delta, delta, delta),
            p_max: self.p_max + Vector3::new(delta, delta, delta),
        }
    }

    /// Returns the coordinates of one of the four corners.
    ///
    /// * `corner` -
    pub fn corner(&self, corner: u8) -> Point3<T>
    where
        T: Copy,
    {
        debug_assert!(corner <= 4);
        let x = corner & 1;
        let y = if corner & 2 == 0 { 0 } else { 1 };
        let z = if corner & 4 == 0 { 0 } else { 1 };
        Point3::new(self[x].x, self[y].y, self[z].z)
    }

    /// Returns the near and far ray parameters where it intersects the bounding
    /// box. If no intersection occurs `None` is returned.
    ///
    /// * `ray` - The ray
    pub fn intersect_p(&self, ray: &Ray) -> Option<(Float, Float)>
    where
        T: num_traits::Float + Copy + PartialOrd + Into<Float>,
    {
        let mut t0 = 0.0;
        let mut t1 = ray.t_max;

        for i in 0..3 {
            // Update interval for ith bounding box slab
            let inv_ray_dir = 1.0 / ray.d[i];

            let mut t_near = (self.p_min[i].into() - ray.o[i]) * inv_ray_dir;
            let mut t_far = (self.p_max[i].into() - ray.o[i]) * inv_ray_dir;

            // Update parametric interval from slab intersection values
            if t_near > t_far {
                swap(&mut t_near, &mut t_far);
            }

            // Update tFar to ensure robust ray–bounds intersection
            t0 = if t_near > t0 { t_near } else { t0 };
            t1 = if t_far < t1 { t_far } else { t1 };
            if t0 > t1 {
                return None;
            }
        }

        Some((t0, t1))
    }

    /// Uses the reciprocal of a rays direction and returns `true` if it
    /// intersects the bounding box; otherwise `false`.
    ///
    /// * `ray`        - The ray.
    /// * `inv_dir`    - Reciprocal of `ray`'s direction.
    /// * `dir_is_neg` - Ray direction is negative.
    #[rustfmt::skip]
    pub fn intersect_p_inv(&self, ray: &Ray, inv_dir: &Vector3f, dir_is_neg: [u8; 3]) -> bool
    where
        T: num_traits::Float + Copy + PartialOrd + Into<Float>,
    {
        let bounds = *self;

        // Check for ray intersection against and slabs
        let mut t_min   = (bounds[    dir_is_neg[0]].x.into() - ray.o.x) * inv_dir.x;
        let mut t_max   = (bounds[1 - dir_is_neg[0]].x.into() - ray.o.x) * inv_dir.x;
        let t_y_min     = (bounds[    dir_is_neg[1]].y.into() - ray.o.y) * inv_dir.y;
        let mut t_y_max = (bounds[1 - dir_is_neg[1]].y.into() - ray.o.y) * inv_dir.y;

        // Update t_max and t_y_max to ensure robust bounds intersection
        let gamma_3 = gamma(3);
        t_max   *= 1.0 + 2.0 * gamma_3;
        t_y_max *= 1.0 + 2.0 * gamma_3;

        if t_min > t_y_max || t_y_min > t_max { return false; }

        if t_y_min > t_min { t_min = t_y_min; }
        if t_y_max < t_max { t_max = t_y_max; }

        // Check for ray intersection against slab
        let t_z_min = (bounds[    dir_is_neg[2]].z.into() - ray.o.z) * inv_dir.z;
        let t_z_max = (bounds[1 - dir_is_neg[2]].z.into() - ray.o.z) * inv_dir.z;

        // Update t_z_max to ensure robust bounds intersection
        if t_min > t_z_max || t_z_min > t_max { return false; }

        if t_z_min > t_min { t_min = t_z_min; }
        if t_z_max < t_max { t_max = t_z_max; }

        t_min < ray.t_max && t_max > 0.0
    }
}

impl<T: Num> Index<u8> for Bounds3<T> {
    type Output = Point3<T>;

    /// Index the minimum and maximum bounds.
    ///
    /// * `i` - 0 for minimum and 1 for maximum.
    fn index(&self, index: u8) -> &Self::Output {
        match index {
            0 => &self.p_min,
            1 => &self.p_max,
            _ => panic!("Invalid index for std::Index on Bounds3<T>"),
        }
    }
}

impl<T: Num + PartialOrd + Copy> Union<Point3<T>> for Bounds3<T> {
    /// Return a bounding box containing the itself and a point.
    ///
    /// * `other` - The point.
    fn union(&self, other: &Point3<T>) -> Self {
        Bounds3 {
            p_min: Point3::new(
                min(self.p_min.x, other.x),
                min(self.p_min.y, other.y),
                min(self.p_min.z, other.z),
            ),
            p_max: Point3::new(
                max(self.p_max.x, other.x),
                max(self.p_max.y, other.y),
                max(self.p_max.z, other.z),
            ),
        }
    }
}

impl<T: Num + PartialOrd + Copy> Union<Bounds3<T>> for Bounds3<T> {
    /// Return a bounding box containing both bounding boxes.
    ///
    /// * `other` - The other bounding box.
    fn union(&self, other: &Bounds3<T>) -> Self {
        Bounds3 {
            p_min: Point3::new(
                min(self.p_min.x, other.p_min.x),
                min(self.p_min.y, other.p_min.y),
                min(self.p_min.z, other.p_min.z),
            ),
            p_max: Point3::new(
                max(self.p_max.x, other.p_max.x),
                max(self.p_max.y, other.p_max.y),
                max(self.p_max.z, other.p_max.z),
            ),
        }
    }
}

impl<T: Num + PartialOrd + Copy> Intersect<Bounds3<T>> for Bounds3<T> {
    /// Return a bounding box containing the intersection of both bounding boxes.
    ///
    /// * `other` - The other bounding box.
    fn intersect(&self, other: &Bounds3<T>) -> Self {
        Bounds3 {
            p_min: Point3::new(
                max(self.p_min.x, other.p_min.x),
                max(self.p_min.y, other.p_min.y),
                max(self.p_min.z, other.p_min.z),
            ),
            p_max: Point3::new(
                min(self.p_max.x, other.p_max.x),
                min(self.p_max.y, other.p_max.y),
                min(self.p_max.z, other.p_max.z),
            ),
        }
    }
}

impl<T: Num + fmt::Debug> fmt::Debug for Bounds3<T> {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Bounds2")
            .field("p_min", &self.p_min)
            .field("p_max", &self.p_max)
            .finish()
    }
}

impl<T: Num + fmt::Display> fmt::Display for Bounds3<T> {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{{}, {}}}", self.p_min, self.p_max)
    }
}
