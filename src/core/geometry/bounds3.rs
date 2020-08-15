//! 3-D Axis Aligned Bounding Boxes.

#![allow(dead_code)]
use super::{
    lerp, max, min, point3, vector3, Axis, Float, Int, Intersect, Point3, Point3f, Union, Vector3,
};
use num_traits::bounds::Bounded;
use num_traits::{Num, Zero};
use std::ops::{DivAssign, Index, Mul};

/// 3-D Axis Aligned Bounding Box.
#[derive(Copy, Clone, Debug, PartialEq)]
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

/// Creates a new 3-D bounding box from 2 points. The minimum and maximum bounds
/// are used for each coordinate axis.
///
/// * `p1` - First point.
/// * `p2` - Second point.
pub fn bounds3<T: Num + PartialOrd + Copy>(p1: Point3<T>, p2: Point3<T>) -> Bounds3<T> {
    Bounds3 {
        p_min: point3(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z)),
        p_max: point3(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z)),
    }
}

/// Returns a 3-D bounding box where minimum and maximum bounds are maximum and
/// minimum values respectively of the type's limits. This is so we can easily grow
/// the bounding box from nothing iteratively.
pub fn empty_bounds2<T: Num + Bounded + PartialOrd + Copy>() -> Bounds3<T> {
    // Don't call bounds2<T>() because it'll create the largest bounding box
    // by flipping p_min and p_max.
    Bounds3 {
        p_min: point3(T::max_value(), T::max_value(), T::max_value()),
        p_max: point3(T::min_value(), T::min_value(), T::min_value()),
    }
}

impl<T: Num + PartialOrd + Copy> From<Point3<T>> for Bounds3<T> {
    /// Use a 3-D point as minimum and maximum 3-D bounds.
    ///
    /// * `p` - 3-D point.
    fn from(p: Point3<T>) -> Self {
        Bounds3 { p_min: p, p_max: p }
    }
}

impl<T: Num + Copy> Bounds3<T> {
    /// Returns true if the bounds describes an empty box where any the
    /// components of any p_max are less than p_max.
    fn is_empty(&self) -> bool
    where
        T: PartialOrd,
    {
        self.p_max.x < self.p_min.x || self.p_max.y < self.p_min.y || self.p_max.z < self.p_min.z
    }

    /// Returns the vector along the box diagonal from the minimum point to
    /// the maximum point.
    fn diagonal(&self) -> Vector3<T> {
        self.p_max - self.p_min
    }

    /// Returns the surface area of the bounding box.
    fn surface_area(&self) -> T
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
    fn volume(&self) -> T
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
    fn maximum_extent(&self) -> Axis
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
    fn overlaps(&self, other: &Self) -> bool
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
    fn offset(&self, p: &Point3<T>) -> Vector3<T>
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
    fn contains(&self, p: &Point3<T>) -> bool
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
    fn contains_exclusive(&self, p: &Point3<T>) -> bool
    where
        T: PartialOrd,
    {
        (p.x >= self.p_min.x && p.x < self.p_max.x)
            && (p.y >= self.p_min.y && p.y < self.p_max.y)
            && (p.z >= self.p_min.z && p.z < self.p_max.z)
    }

    /// Return the center and radius of a sphere bounded on the corners of the
    /// bounding box.
    fn bounding_sphere(&self) -> (Point3<T>, T)
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
    fn lerp(&self, t: &Point3f) -> Point3<T>
    where
        Float: Mul<T, Output = T>,
    {
        point3(
            lerp::<T>(t.x, self.p_min.x, self.p_max.x),
            lerp::<T>(t.y, self.p_min.y, self.p_max.y),
            lerp::<T>(t.z, self.p_min.z, self.p_max.z),
        )
    }

    /// Pad the bounding box by a constant factor in all dimensions.
    ///
    /// * `delta` - Padding amount.
    fn expand(&self, delta: T) -> Bounds3<T>
    where
        T: PartialOrd,
    {
        // Don't call bounds2<T>() to prevent flipping p_min and p_max when
        // the input is empty box.
        Bounds3 {
            p_min: self.p_min - vector3(delta, delta, delta),
            p_max: self.p_max + vector3(delta, delta, delta),
        }
    }

    /// Returns the coordinates of one of the four corners.
    ///
    /// * `corner` -
    fn corner(&self, corner: u8) -> Point3<T>
    where
        T: Copy,
    {
        debug_assert!(corner <= 4);
        let x = corner & 1;
        let y = if corner & 2 == 0 { 0 } else { 1 };
        let z = if corner & 4 == 0 { 0 } else { 1 };
        point3(self[x].x, self[y].y, self[z].z)
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
            p_min: point3(
                min(self.p_min.x, other.x),
                min(self.p_min.y, other.y),
                min(self.p_min.z, other.z),
            ),
            p_max: point3(
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
            p_min: point3(
                min(self.p_min.x, other.p_min.x),
                min(self.p_min.y, other.p_min.y),
                min(self.p_min.z, other.p_min.z),
            ),
            p_max: point3(
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
            p_min: point3(
                max(self.p_min.x, other.p_min.x),
                max(self.p_min.y, other.p_min.y),
                max(self.p_min.z, other.p_min.z),
            ),
            p_max: point3(
                min(self.p_max.x, other.p_max.x),
                min(self.p_max.y, other.p_max.y),
                min(self.p_max.z, other.p_max.z),
            ),
        }
    }
}
