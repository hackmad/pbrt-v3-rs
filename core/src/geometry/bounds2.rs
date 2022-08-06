//! 2D Axis Aligned Bounding Boxes.

#![allow(dead_code)]
use crate::geometry::*;
use crate::pbrt::*;
use itertools::{iproduct, Product};
use num_traits::bounds::Bounded;
use num_traits::{Num, Zero};
use std::fmt;
use std::ops::{DivAssign, Index, Mul, Range};

/// 2D Axis Aligned Bounding Box.
#[derive(Copy, Clone, Default, PartialEq)]
pub struct Bounds2<T: Num> {
    /// Minimum bounds.
    pub p_min: Point2<T>,

    /// Maximum bounds.
    pub p_max: Point2<T>,
}

/// 2-D bounding box containing `Float` points.
pub type Bounds2f = Bounds2<Float>;
impl Bounds2f {
    /// 2-D bounding box where minimum and maximum bounds are maximum and minimum
    /// floating point values. This is so we can easily grow the bounding box
    /// from nothing iteratively.
    pub const EMPTY: Self = Self {
        p_min: Point2f::MAX,
        p_max: Point2f::MIN,
    };
}

/// 2-D bounding box containing `Int` points.
pub type Bounds2i = Bounds2<Int>;
impl Bounds2i {
    /// 2-D bounding box where minimum and maximum bounds are maximum and minimum
    /// floating point values. This is so we can easily grow the bounding box
    /// from nothing iteratively.
    pub const EMPTY: Self = Self {
        p_min: Point2i::MAX,
        p_max: Point2i::MIN,
    };
}

impl<T: Num + PartialOrd + Copy> From<Point2<T>> for Bounds2<T> {
    /// Use a 2D point as minimum and maximum 2D bounds.
    ///
    /// * `p` - 2D point.
    fn from(p: Point2<T>) -> Self {
        Bounds2 { p_min: p, p_max: p }
    }
}

impl From<Bounds2i> for Bounds2f {
    /// Convert a `Bounds2i` to `Bounds2f`.
    ///
    /// * `b` - The `Bounds2i` to convert.
    fn from(b: Bounds2i) -> Self {
        Self {
            p_min: b.p_min.into(),
            p_max: b.p_max.into(),
        }
    }
}

impl From<Bounds2f> for Bounds2i {
    /// Convert a `Bounds2f` to `Bounds2i`.
    ///
    /// * `b` - The `Bounds2f` to convert.
    fn from(b: Bounds2f) -> Self {
        Self {
            p_min: b.p_min.into(),
            p_max: b.p_max.into(),
        }
    }
}

impl<T: Num + Copy> Bounds2<T> {
    /// Creates a new 2D bounding box from 2 points. The minimum and maximum bounds
    /// are used for each coordinate axis.
    ///
    /// * `p1` - First point.
    /// * `p2` - Second point.
    pub fn new(p1: Point2<T>, p2: Point2<T>) -> Self
    where
        T: PartialOrd + Copy,
    {
        Self {
            p_min: Point2::new(min(p1.x, p2.x), min(p1.y, p2.y)),
            p_max: Point2::new(max(p1.x, p2.x), max(p1.y, p2.y)),
        }
    }

    /// Returns a 2D bounding box where minimum and maximum bounds are maximum and
    /// minimum values respectively of the type's limits. This is so we can easily
    /// grow the boundng box from nothing iteratively.
    pub fn empty() -> Self
    where
        T: Bounded + PartialOrd + Copy,
    {
        // Don't call new() because it'll create the largest bounding box
        // by flipping p_min and p_max.
        Self {
            p_min: Point2::new(T::max_value(), T::max_value()),
            p_max: Point2::new(T::min_value(), T::min_value()),
        }
    }

    /// Returns true if the bounding box describes an empty box where any the
    /// components of any p_max are less than p_max.
    pub fn is_empty(&self) -> bool
    where
        T: PartialOrd,
    {
        self.p_max.x < self.p_min.x || self.p_max.y < self.p_min.y
    }

    /// Returns the vector along the box diagonal from the minimum point to
    /// the maximum point.
    pub fn diagonal(&self) -> Vector2<T> {
        self.p_max - self.p_min
    }

    /// Returns the area of the bounding box.
    pub fn area(&self) -> T
    where
        T: PartialOrd,
    {
        if self.is_empty() {
            T::zero()
        } else {
            let d = self.diagonal();
            d.x * d.y
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
        if d.x > d.y {
            Axis::X
        } else {
            Axis::Y
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
        x && y
    }

    /// Returns the continuous position of a point relative to the corners of the
    /// box, where a point at the minimum corner has offset `(0, 0)` and a
    /// point at the maximum corner has offset is `(1, 1)`.
    ///
    /// * `p` - The point.
    pub fn offset(&self, p: &Point2<T>) -> Vector2<T>
    where
        T: num_traits::Float + DivAssign<T> + PartialOrd + Copy,
    {
        let mut o = p - self.p_min;
        if self.p_max.x > self.p_min.x {
            o.x /= self.p_max.x - self.p_min.x;
        }
        if self.p_max.y > self.p_min.y {
            o.y /= self.p_max.y - self.p_min.y;
        }
        o
    }

    /// Returns true if a point is inside the bounding box.
    ///
    /// * `p` - The point.
    pub fn contains(&self, p: &Point2<T>) -> bool
    where
        T: PartialOrd,
    {
        p.x >= self.p_min.x && p.x <= self.p_max.x && p.y >= self.p_min.y && p.y <= self.p_max.y
    }

    /// Returns true if a point is inside the bounding box. The upper boundary
    /// is considered out of bounds. This is useful for integer-typed bounds.
    ///
    /// * `p` - The point.
    pub fn contains_exclusive(&self, p: &Point2<T>) -> bool
    where
        T: PartialOrd,
    {
        p.x >= self.p_min.x && p.x < self.p_max.x && p.y >= self.p_min.y && p.y < self.p_max.y
    }

    /// Return the center and radius of a circle bounded on the corners of the
    /// bounding box.
    pub fn bounding_circle(&self) -> (Point2<T>, T)
    where
        T: num_traits::Float + Zero,
        Float: Mul<Point2<T>, Output = Point2<T>>,
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
    pub fn lerp(&self, t: &Point2f) -> Point2<T>
    where
        Float: Mul<T, Output = T>,
    {
        Point2::new(
            lerp::<T>(t.x, self.p_min.x, self.p_max.x),
            lerp::<T>(t.y, self.p_min.y, self.p_max.y),
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
            p_min: self.p_min - Vector2::new(delta, delta),
            p_max: self.p_max + Vector2::new(delta, delta),
        }
    }

    /// Returns the coordinates of one of the four corners.
    ///
    /// * `corner` -
    pub fn corner(&self, corner: u8) -> Point2<T>
    where
        T: Copy,
    {
        debug_assert!(corner <= 4);
        let x = corner & 1;
        let y = if corner & 2 == 0 { 0 } else { 1 };
        Point2::new(self[x].x, self[y].y)
    }
}

impl<T: Num> Index<u8> for Bounds2<T> {
    type Output = Point2<T>;

    /// Index the minimum and maximum bounds.
    ///
    /// * `i` - 0 for minimum and 1 for maximum.
    fn index(&self, index: u8) -> &Self::Output {
        match index {
            0 => &self.p_min,
            1 => &self.p_max,
            _ => panic!("Invalid index for std::Index on Bounds2<T>"),
        }
    }
}

impl<T: Num + PartialOrd + Copy> Union<Point2<T>> for Bounds2<T> {
    /// Return a bounding box containing the itself and a point.
    ///
    /// * `other` - The point.
    fn union(&self, other: &Point2<T>) -> Self {
        Bounds2 {
            p_min: Point2::new(min(self.p_min.x, other.x), min(self.p_min.y, other.y)),
            p_max: Point2::new(max(self.p_max.x, other.x), max(self.p_max.y, other.y)),
        }
    }
}

impl<T: Num + PartialOrd + Copy> Union<Bounds2<T>> for Bounds2<T> {
    /// Return a bounding box containing both bounding boxes.
    ///
    /// * `other` - The other bounding box.
    fn union(&self, other: &Bounds2<T>) -> Self {
        Bounds2 {
            p_min: Point2::new(
                min(self.p_min.x, other.p_min.x),
                min(self.p_min.y, other.p_min.y),
            ),
            p_max: Point2::new(
                max(self.p_max.x, other.p_max.x),
                max(self.p_max.y, other.p_max.y),
            ),
        }
    }
}

impl<T: Num + PartialOrd + Copy> Intersect<Bounds2<T>> for Bounds2<T> {
    /// Return a bounding box containing the intersection of both bounding boxes.
    ///
    /// * `other` - The other bounding box.
    fn intersect(&self, other: &Bounds2<T>) -> Self {
        Bounds2 {
            p_min: Point2::new(
                max(self.p_min.x, other.p_min.x),
                max(self.p_min.y, other.p_min.y),
            ),
            p_max: Point2::new(
                min(self.p_max.x, other.p_max.x),
                min(self.p_max.y, other.p_max.y),
            ),
        }
    }
}

impl<T: Num + fmt::Debug> fmt::Debug for Bounds2<T> {
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

impl<T: Num + fmt::Display> fmt::Display for Bounds2<T> {
    /// Formats the value using the given formatter.
    ///
    /// * `f` - Formatter.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{{}, {}}}", self.p_min, self.p_max)
    }
}

/// An iterator that can step through integer coordinates in a bounding box
/// in a left-to-right (x-axis) and top-to-bottom (y-axis) scan order.
pub struct Bounds2iIterator {
    /// The iterator that will supply points.
    p: Product<Range<Int>, Range<Int>>,
}

impl IntoIterator for Bounds2i {
    type Item = Point2i;
    type IntoIter = Bounds2iIterator;

    /// Create an iterator for `Bounds2i`.
    fn into_iter(self) -> Self::IntoIter {
        let max_x = if self.p_min.x == self.p_max.x {
            self.p_max.x + 1
        } else {
            self.p_max.x
        };
        let max_y = if self.p_min.y == self.p_max.y {
            self.p_max.y + 1
        } else {
            self.p_max.y
        };
        Bounds2iIterator {
            p: iproduct!(self.p_min.y..max_y, self.p_min.x..max_x),
        }
    }
}

impl Iterator for Bounds2iIterator {
    type Item = Point2i;

    /// Get the next point.
    fn next(&mut self) -> Option<Self::Item> {
        self.p.next().map(|(y, x)| Point2i::new(x, y))
    }
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
#[macro_use]
mod tests {
    use crate::geometry::*;
    use crate::pbrt::*;
    use float_cmp::*;
    use num_traits::Zero;
    use proptest::prelude::*;

    #[test]
    fn empty_bounds2f_returns_min_greater_than_max_components() {
        let b = Bounds2::<f32>::empty();
        assert_eq!(b.p_min, Point2::new(f32::MAX, f32::MAX));
        assert_eq!(b.p_max, Point2::new(f32::MIN, f32::MIN));
    }

    #[test]
    fn empty_bounds2i_returns_min_greater_than_max_components() {
        let b = Bounds2::<i32>::empty();
        assert_eq!(b.p_min, Point2::new(i32::MAX, i32::MAX));
        assert_eq!(b.p_max, Point2::new(i32::MIN, i32::MIN));
    }

    #[test]
    fn bounding_circle_of_empty_box_returns_origin_and_zero_radius() {
        let (center, radius) = Bounds2::<f32>::empty().bounding_circle();
        assert_eq!(center, Point2::<f32>::zero());
        assert_eq!(radius, 0.0);
    }

    #[test]
    fn corner_returns_points_with_left_to_right_in_x_top_to_bottom_in_y() {
        let p1 = Point2::new(-1, -1);
        let p2 = Point2::new(1, 1);
        let b = Bounds2::new(p1, p2);
        assert_eq!(b.corner(0), Point2::new(p1.x, p1.y));
        assert_eq!(b.corner(1), Point2::new(p2.x, p1.y));
        assert_eq!(b.corner(2), Point2::new(p1.x, p2.y));
        assert_eq!(b.corner(3), Point2::new(p2.x, p2.y));
    }

    #[test]
    #[should_panic]
    #[allow(unused)]
    fn invalid_index_panics() {
        let p = Bounds2::<i32>::empty()[2];
    }

    #[test]
    fn index_returns_p_min_at_0_and_p_max_at_1() {
        let p1 = Point2::new(-1, -1);
        let p2 = Point2::new(1, 1);
        let b = Bounds2::new(p1, p2);
        assert_eq!(b[0], p1);
        assert_eq!(b[1], p2);
    }

    #[test]
    fn area_of_empty_bounds2i_returns_zero() {
        let empty = Bounds2::<i32>::empty();
        assert_eq!(empty.area(), 0);
    }

    #[test]
    fn area_of_empty_bounds2f_returns_zero() {
        let empty = Bounds2::<f32>::empty();
        assert_eq!(empty.area(), 0.0);
    }

    #[test]
    fn union_of_two_empty_bounds2i_returns_empty() {
        let empty = Bounds2::<i32>::empty();
        assert!(empty.union(&empty).is_empty());
    }

    #[test]
    fn union_of_two_empty_bounds2f_returns_empty() {
        let empty = Bounds2::<f32>::empty();
        assert!(empty.union(&empty).is_empty());
    }

    #[test]
    fn intersection_of_two_empty_bounds2i_returns_empty() {
        let empty = Bounds2::<i32>::empty();
        assert!(empty.intersect(&empty).is_empty());
    }

    #[test]
    fn intersection_of_two_empty_bounds2f_returns_empty() {
        let empty = Bounds2::<f32>::empty();
        assert!(empty.intersect(&empty).is_empty());
    }

    #[test]
    fn iterating_empty_bounds2i_return_none() {
        let empty = Bounds2::<i32>::empty();
        let mut iter = empty.into_iter();
        assert!(iter.next().is_none());
    }

    #[test]
    fn iterate_point_bounds2i_returns_point_only() {
        let b = Bounds2::new(Point2::new(0, 0), Point2::new(0, 0));
        let mut iter = b.into_iter();
        assert_eq!(iter.next(), Some(Point2::new(0, 0)));
        assert!(iter.next().is_none());
    }

    // Define some properties for tests.
    prop_range!(range_i32, i32, -100..100i32);
    prop_range!(range_f32, f32, -100.0..100.0f32);

    prop_non_zero_range!(non_zero_i32, i32, -100..100i32);
    prop_non_zero_range!(non_zero_f32, f32, -100.0..100.0f32);

    prop_point2!(point2_i32, i32, -100..100i32, -100..100i32);
    prop_point2!(point2_f32, f32, -100.0..100.0f32, -100.0..100.0f32);

    proptest! {
        #[test]
        fn bounds2i_sorts_x_and_y_components(p1 in point2_i32(), p2 in point2_i32()) {
            let b1 = Bounds2::new(p1, p2);
            let b2 = Bounds2::new(p2, p1);
            prop_assert_eq!(b1, b2);
            prop_assert_eq!(b1.p_min, Point2::new(min(p1.x, p2.x), min(p1.y, p2.y)));
            prop_assert_eq!(b1.p_max, Point2::new(max(p1.x, p2.x), max(p1.y, p2.y)));
        }

        #[test]
        fn bounds2f_sorts_x_and_y_components(p1 in point2_f32(), p2 in point2_f32()) {
            let b1 = Bounds2::new(p1, p2);
            let b2 = Bounds2::new(p2, p1);
            prop_assert_eq!(b1, b2);
            prop_assert_eq!(b1.p_min, Point2::new(min(p1.x, p2.x), min(p1.y, p2.y)));
            prop_assert_eq!(b1.p_max, Point2::new(max(p1.x, p2.x), max(p1.y, p2.y)));
        }

        #[test]
        fn bounds2i_from_point_sets_min_max_to_given_point(p in point2_i32()) {
            let b = Bounds2::from(p);
            prop_assert_eq!(b.p_min, p);
            prop_assert_eq!(b.p_max, p);
        }

        #[test]
        fn bounds2f_from_point_sets_min_max_to_given_point(p in point2_f32()) {
            let b = Bounds2::from(p);
            prop_assert_eq!(b.p_min, p);
            prop_assert_eq!(b.p_max, p);
        }

        #[test]
        fn diagonal_of_bounds2i_returns_vector_from_min_to_max(
            p1 in point2_i32(), p2 in point2_i32(),
        ) {
            let b = Bounds2::new(p1, p2);
            prop_assert_eq!(b.diagonal(), b.p_max - b.p_min);
        }

        #[test]
        fn diagonal_of_bounds2f_returns_vector_from_min_to_max(
            p1 in point2_f32(), p2 in point2_f32(),
        ) {
            let b = Bounds2::new(p1, p2);
            prop_assert_eq!(b.diagonal(), b.p_max - b.p_min);
        }

        #[test]
        fn area_of_non_empty_bounds2i_returns_product_of_diagonal_components(
            p in point2_i32(), dx in -10..10i32, dy in -10..10i32,
        ) {
            let b = Bounds2::new(p, p + Vector2::new(dx, dy));
            prop_assert_eq!(b.area(), abs(dx * dy));
        }

        #[test]
        fn area_of_non_empty_bounds2f_returns_product_of_diagonal_components(
            p in point2_f32(), dx in -10.0..10.0f32, dy in -10.0..10.0f32,
        ) {
            let b = Bounds2::new(p, p + Vector2::new(dx, dy));
            prop_assert!(approx_eq!(f32, b.area(), abs(dx * dy), epsilon = 0.0001));
        }

        #[test]
        fn maximum_extent_of_non_empty_bounds2i_returns_axis_with_max_diagonal_component(
            p in point2_i32(), d in 0..10i32,
        ) {
            prop_assert_eq!(Bounds2::new(p, p + Vector2::new(d + 1, d)).maximum_extent(), Axis::X);
            prop_assert_eq!(Bounds2::new(p, p + Vector2::new(d, d + 1)).maximum_extent(), Axis::Y);
            prop_assert_eq!(Bounds2::new(p, p + Vector2::new(d, d)).maximum_extent(), Axis::Y);
        }

        #[test]
        fn maximum_extent_of_non_empty_bounds2f_returns_axis_with_max_diagonal_component(
            p in point2_f32(), d in 0.0..10.0f32,
        ) {
            prop_assert_eq!(Bounds2::new(p, p + Vector2::new(d + 0.001, d)).maximum_extent(), Axis::X);
            prop_assert_eq!(Bounds2::new(p, p + Vector2::new(d, d + 0.001)).maximum_extent(), Axis::Y);

            // Floating point errors due to addition with x/y components will cause
            // this assertion to fail.
            // prop_assert_eq!(Bounds2::new(p, p + Vector2::new(d, d)).maximum_extent(), Axis::Y);
        }

        #[test]
        fn maximum_extent_of_non_empty_bounds2f_returns_y_axis_edge_case(
            x in -10.0..10.0f32, d in 0.0..10.0f32,
        ) {
            prop_assert_eq!(Bounds2::new(Point2::new(x, x), Point2::new(x + d, x + d)).maximum_extent(), Axis::Y);
        }

        #[test]
        fn overlaps_returns_true_when_two_bounds2f_overlap(
            p in point2_f32(),
            dx in 0.1..1.0f32, dy in 0.1..1.0f32,
            sx in 0.0..2.0f32, sy in 0.0..2.0f32,
        ) {
            let b1 = Bounds2::new(p - Vector2::new(dx, dy), p + Vector2::new(dx, dy));
            let bounds = vec!(
                Bounds2::new(b1.p_min - Vector2::new( sx * dx,      0.0), b1.p_max - Vector2::new( sx * dx,      0.0)),
                Bounds2::new(b1.p_min + Vector2::new( sx * dx,      0.0), b1.p_max + Vector2::new( sx * dx,      0.0)),
                Bounds2::new(b1.p_min - Vector2::new(     0.0,  sy * dy), b1.p_max - Vector2::new(     0.0,  sy * dy)),
                Bounds2::new(b1.p_min + Vector2::new(     0.0,  sy * dy), b1.p_max + Vector2::new(     0.0,  sy * dy)),
                Bounds2::new(b1.p_min - Vector2::new( sx * dx,  sy * dy), b1.p_max - Vector2::new( sx * dx,  sy * dy)),
                Bounds2::new(b1.p_min + Vector2::new( sx * dx,  sy * dy), b1.p_max + Vector2::new( sx * dx,  sy * dy)),
                Bounds2::new(b1.p_min + Vector2::new( sx * dx, -sy * dy), b1.p_max + Vector2::new( sx * dx, -sy * dy)),
                Bounds2::new(b1.p_min + Vector2::new(-sx * dx,  sy * dy), b1.p_max + Vector2::new(-sx * dx,  sy * dy)),
                Bounds2::new(b1.p_min - Vector2::new( sx * dx,      0.0), b1.p_max + Vector2::new( sx * dx,      0.0)),
                Bounds2::new(b1.p_min - Vector2::new(     0.0,  sy * dy), b1.p_max + Vector2::new(     0.0,  sy * dy)),
                Bounds2::new(b1.p_min - Vector2::new( sx * dx,  sy * dy), b1.p_max + Vector2::new( sx * dx,  sy * dy)),
            );
            for b2 in bounds {
                prop_assert!(b1.overlaps(&b2));
            }
        }

        #[test]
        fn overlaps_returns_false_when_two_bounds2f_do_not_overlap(
            p in point2_f32(),
            dx in 1.0..2.0f32, dy in 1.0..2.0f32,
            sx in 2.001..3.0f32, sy in 2.001..3.0f32,
        ) {
            let b1 = Bounds2::new(p - Vector2::new(dx, dy), p + Vector2::new(dx, dy));
            let bounds = vec!(
                Bounds2::new(b1.p_min - Vector2::new( sx * dx,      0.0), b1.p_max - Vector2::new( sx * dx,      0.0)),
                Bounds2::new(b1.p_min + Vector2::new( sx * dx,      0.0), b1.p_max + Vector2::new( sx * dx,      0.0)),
                Bounds2::new(b1.p_min - Vector2::new(     0.0,  sy * dy), b1.p_max - Vector2::new(     0.0,  sy * dy)),
                Bounds2::new(b1.p_min + Vector2::new(     0.0,  sy * dy), b1.p_max + Vector2::new(     0.0,  sy * dy)),
                Bounds2::new(b1.p_min - Vector2::new( sx * dx,  sy * dy), b1.p_max - Vector2::new( sx * dx,  sy * dy)),
                Bounds2::new(b1.p_min + Vector2::new( sx * dx, -sy * dy), b1.p_max + Vector2::new( sx * dx, -sy * dy)),
                Bounds2::new(b1.p_min + Vector2::new(-sx * dx,  sy * dy), b1.p_max + Vector2::new(-sx * dx,  sy * dy)),
            );
            for b2 in bounds {
                prop_assert!(!b1.overlaps(&b2));
            }
        }

        #[test]
        fn offset_of_any_point_within_an_empty_bounds2f_returns_vector_towards_p_min(
            p in point2_f32()
        ) {
            let b = Bounds2::<f32>::empty();
            prop_assert_eq!(b.offset(&p), p - b.p_min);
        }

        #[test]
        fn offset_of_any_point_within_non_empty_bounds2f_returns_components_in_0_1(
            p1 in point2_f32(), p2 in point2_f32(), p in point2_f32(),
        ) {
            let b = Bounds2::new(p1, p2);
            let offset = b.offset(&p);
            let o = p - b.p_min;
            prop_assert_eq!(offset.x, o.x / (b.p_max.x - b.p_min.x));
            prop_assert_eq!(offset.y, o.y / (b.p_max.y - b.p_min.y));
        }

        #[test]
        fn contains_returns_false_for_any_point_when_empty_bounds2f(p in point2_f32()) {
            prop_assert!(!Bounds2::<f32>::empty().contains(&p));
        }

        #[test]
        fn contains_returns_true_for_point_bounds2f_when_p_min_p_max_is_same_point(
            p in point2_f32(), dx in 0.001..1.0f32, dy in 0.001..1.0f32,
        ) {
            let b = Bounds2::from(p);
            prop_assert!(b.contains(&p));
            prop_assert!(!b.contains(&(p + Vector2::new( dx, 0.0))));
            prop_assert!(!b.contains(&(p + Vector2::new(0.0,  dy))));
            prop_assert!(!b.contains(&(p + Vector2::new( dx,  dy))));
            prop_assert!(!b.contains(&(p + Vector2::new(-dx, 0.0))));
            prop_assert!(!b.contains(&(p + Vector2::new(0.0, -dy))));
            prop_assert!(!b.contains(&(p + Vector2::new(-dx, -dy))));
        }

        #[test]
        fn contains_returns_false_when_point_is_outside_bounds2f(
            p in point2_f32(),
            dx in 0.001..1.0f32, dy in 0.001..1.0f32,
            sx in 0.001..1.0f32, sy in 0.001..1.0f32,
        ) {
            let b = Bounds2::new(p - Vector2::new(dx, dy), p + Vector2::new(dx, dy));

            let points = vec!(
                Point2::new(b.p_min.x - sx, b.p_min.y     ),
                Point2::new(b.p_min.x,      b.p_min.y - sy),
                Point2::new(b.p_min.x - sx, b.p_min.y - sy),
                Point2::new(b.p_max.x + sx, b.p_max.y     ),
                Point2::new(b.p_max.x,      b.p_max.y + sy),
                Point2::new(b.p_max.x + sx, b.p_max.y + sy),
                Point2::new(b.p_min.x - sx, b.p_max.y     ),
                Point2::new(b.p_min.x,      b.p_max.y + sy),
                Point2::new(b.p_min.x - sx, b.p_max.y + sy),
                Point2::new(b.p_max.x + sx, b.p_min.y     ),
                Point2::new(b.p_max.x,      b.p_min.y - sy),
                Point2::new(b.p_max.x + sx, b.p_min.y - sy),
            );
            for p in points {
                prop_assert!(!b.contains(&p));
            }
        }

        #[test]
        fn contains_returns_true_when_point_is_inside_bounds2f(
            p in point2_f32(),
            dx in 0.001..1.0f32, dy in 0.001..1.0f32,
            tx in 0.0..1.0f32, ty in 0.0..1.0f32,
        ) {
            let b = Bounds2::new(p - Vector2::new(dx, dy), p + Vector2::new(dx, dy));
            let p = Point2::new(
                lerp(tx, b.p_min.x, b.p_max.x),
                lerp(ty, b.p_min.y, b.p_max.y),
            );
            prop_assert!(b.contains(&p));
        }

        #[test]
        fn contains_exclusive_returns_false_for_any_point_when_empty_bounds2i(
            p in point2_i32(),
        ) {
            let b = Bounds2::<i32>::empty();
            prop_assert!(!b.contains_exclusive(&p));
        }

        #[test]
        fn contains_exclusive_returns_true_for_point_bounds2i_when_p_min_p_max_is_same_point(
            p in point2_i32(), dx in 0..5i32, dy in 0..5i32,
        ) {
            let b = Bounds2::from(p);
            prop_assert!(!b.contains_exclusive(&p));
            prop_assert!(!b.contains_exclusive(&(p + Vector2::new( dx,   0))));
            prop_assert!(!b.contains_exclusive(&(p + Vector2::new(  0,  dy))));
            prop_assert!(!b.contains_exclusive(&(p + Vector2::new( dx,  dy))));
            prop_assert!(!b.contains_exclusive(&(p + Vector2::new(-dx,   0))));
            prop_assert!(!b.contains_exclusive(&(p + Vector2::new(  0, -dy))));
            prop_assert!(!b.contains_exclusive(&(p + Vector2::new(-dx, -dy))));
        }

        #[test]
        fn contains_exclusive_returns_false_when_point_is_outside_bounds2i(
            p in point2_i32(),
            dx in 1..10i32, dy in 1..10i32,
            sx in 1..10i32, sy in 1..10i32,
            kx in 0..10i32, ky in 0..10i32,
        ) {
            let b = Bounds2::new(p - Vector2::new(dx, dy), p + Vector2::new(dx, dy));

            let points = vec!(
                Point2::new(b.p_min.x - sx, b.p_min.y     ),
                Point2::new(b.p_min.x,      b.p_min.y - sy),
                Point2::new(b.p_min.x - sx, b.p_min.y - sy),
                Point2::new(b.p_max.x + kx, b.p_max.y     ),
                Point2::new(b.p_max.x,      b.p_max.y + ky),
                Point2::new(b.p_max.x + kx, b.p_max.y + ky),
                Point2::new(b.p_min.x - sx, b.p_max.y     ),
                Point2::new(b.p_min.x,      b.p_max.y + ky),
                Point2::new(b.p_min.x - sx, b.p_max.y + ky),
                Point2::new(b.p_max.x + kx, b.p_min.y     ),
                Point2::new(b.p_max.x,      b.p_min.y - sy),
                Point2::new(b.p_max.x + kx, b.p_min.y - sy),
            );
            for p in points {
                prop_assert!(!b.contains_exclusive(&p));
            }
        }

        #[test]
        fn contains_exclusive_returns_true_when_point_is_inside_bounds2i(
            p in point2_i32(),
            dx in 1..10i32, dy in 1..10i32,
            tx in 0.0..0.999f32, ty in 0.0..0.999f32,
        ) {
            let b = Bounds2::new(p - Vector2::new(dx, dy), p + Vector2::new(dx, dy));
            let px = lerp(tx, b.p_min.x as f32, b.p_max.x as f32);
            let py = lerp(ty, b.p_min.y as f32, b.p_max.y as f32);
            let p = Point2::new(px.floor() as i32, py.floor() as i32);
            prop_assert!(b.contains_exclusive(&p));
        }

        #[test]
        fn bounding_circle_of_point_bounds2f_returns_point_as_center_and_zero_radius(
            p in point2_f32(),
        ) {
            let (center, radius) = Bounds2::from(p).bounding_circle();
            prop_assert_eq!(center, p);
            prop_assert_eq!(radius, 0.0);
        }

        #[test]
        fn bounding_circle_returns_midpoint_as_center_and_distance_to_p_max_as_radius(
            p1 in point2_f32(), p2 in point2_f32(),
        ) {
            let b = Bounds2::new(p1, p2);
            let (center, radius) = b.bounding_circle();
            prop_assert_eq!(center, lerp(0.5, p1, p2));
            prop_assert_eq!(radius, center.distance(b.p_max));
        }

        #[test]
        fn lerp_returns_corners_of_bounds2f_at_0_and_1(
            p1 in point2_f32(), p2 in point2_f32(),
        ) {
            let b = Bounds2::new(p1, p2);
            prop_assert_eq!(b.lerp(&Point2::new(0.0, 0.0)), b.p_min);
            prop_assert_eq!(b.lerp(&Point2::new(1.0, 1.0)), b.p_max);
        }

        #[test]
        fn lerp_interpolates_and_extrapolates_across_corners_of_bounds2f(
            p1 in point2_f32(), p2 in point2_f32(),
            tx in -2.0..2.0f32, ty in -2.0..2.0f32,
        ) {
            let b = Bounds2::new(p1, p2);
            let l = b.lerp(&Point2::new(tx, ty));
            prop_assert_eq!(l.x, lerp(tx, b.p_min.x, b.p_max.x));
            prop_assert_eq!(l.y, lerp(ty, b.p_min.y, b.p_max.y));
        }

        #[test]
        fn expand_returns_empty_when_bounds2i_is_empty(delta in 0..2i32) {
            prop_assert!(Bounds2::<i32>::empty().expand(delta).is_empty());
        }

        #[test]
        fn expand_returns_non_empty_bounds2i_for_bounds2i_from_point(
            p in point2_i32(), delta in 0..2i32,
        ) {
            let b1 = Bounds2::from(p);
            let b2 = b1.expand(delta);
            prop_assert_eq!(b2.p_min.x, p.x - delta);
            prop_assert_eq!(b2.p_min.y, p.y - delta);
            prop_assert_eq!(b2.p_max.x, p.x + delta);
            prop_assert_eq!(b2.p_max.y, p.y + delta);
        }

        #[test]
        fn expand_returns_bounds2i_for_non_empty_bounds2i(
            p1 in point2_i32(), p2 in point2_i32(), delta in 0..2i32,
        ) {
            let b1 = Bounds2::new(p1, p2);
            let b2 = b1.expand(delta);
            prop_assert_eq!(b2.p_min.x, b1.p_min.x - delta);
            prop_assert_eq!(b2.p_min.y, b1.p_min.y - delta);
            prop_assert_eq!(b2.p_max.x, b1.p_max.x + delta);
            prop_assert_eq!(b2.p_max.y, b1.p_max.y + delta);
        }

        #[test]
        fn expand_returns_empty_when_bounds2f_is_empty(delta in 0.0..100.0f32) {
            let b1 = Bounds2::<f32>::empty();
            let b2 = b1.expand(delta);
            prop_assert_eq!(b2.p_min.x, b1.p_min.x - delta);
            prop_assert_eq!(b2.p_min.y, b1.p_min.y - delta);
            prop_assert_eq!(b2.p_max.x, b1.p_max.x + delta);
            prop_assert_eq!(b2.p_max.y, b1.p_max.y + delta);
        }

        #[test]
        fn expand_returns_non_empty_bounds2f_for_bounds2f_from_point(
            p in point2_f32(), delta in 0.0..100.0f32,
        ) {
            let b1 = Bounds2::from(p);
            let b2 = b1.expand(delta);
            prop_assert_eq!(b2.p_min.x, p.x - delta);
            prop_assert_eq!(b2.p_min.y, p.y - delta);
            prop_assert_eq!(b2.p_max.x, p.x + delta);
            prop_assert_eq!(b2.p_max.y, p.y + delta);
        }

        #[test]
        fn expand_returns_non_empty_bounds2f_for_non_empty_bounds2f(
            p1 in point2_f32(), p2 in point2_f32(), delta in 0.0..100.0f32,
        ) {
            let b1 = Bounds2::new(p1, p2);
            let b2 = b1.expand(delta);
            prop_assert_eq!(b2.p_min.x, b1.p_min.x - delta);
            prop_assert_eq!(b2.p_min.y, b1.p_min.y - delta);
            prop_assert_eq!(b2.p_max.x, b1.p_max.x + delta);
            prop_assert_eq!(b2.p_max.y, b1.p_max.y + delta);
        }

        #[test]
        fn union_empty_with_bounds2i_from_point_returns_latter(
            p in point2_i32(),
        ) {
            let empty = Bounds2::<i32>::empty();
            prop_assert_eq!(empty.union(&p), Bounds2::from(p));
        }

        #[test]
        fn union_empty_with_bounds2f_from_point_returns_latter(
            p in point2_f32(),
        ) {
            let empty = Bounds2::<f32>::empty();
            prop_assert_eq!(empty.union(&p), Bounds2::from(p));
        }

        #[test]
        fn union_empty_with_non_empty_bounds2i_returns_non_empty(
            p1 in point2_i32(), p2 in point2_i32(),
        ) {
            let empty = Bounds2::<i32>::empty();
            let non_empty = Bounds2::new(p1, p2);
            prop_assert_eq!(empty.union(&non_empty), non_empty);
            prop_assert_eq!(non_empty.union(&empty), non_empty);
        }

        #[test]
        fn union_empty_with_non_empty_bounds2f_returns_non_empty(
            p1 in point2_f32(), p2 in point2_f32(),
        ) {
            let empty = Bounds2::<f32>::empty();
            let non_empty = Bounds2::new(p1, p2);
            prop_assert_eq!(empty.union(&non_empty), non_empty);
            prop_assert_eq!(non_empty.union(&empty), non_empty);
        }

        #[test]
        fn union_non_empty_bounds2i_with_exterior_point_returns_non_empty_bounds2i(
            p in point2_i32(),
            dx in 1..10i32, dy in 1..10i32,
            s in 1..10i32, t in 0.0..1.0f32,
        ) {
            let v = Vector2::new(dx, dy);
            let b = Bounds2::new(p - v, p + v);

            let y = lerp(t - 1.0, b.p_min.y as f32, b.p_max.y as f32).round() as i32;
            prop_assert_eq!(
                b.union(&Point2::new(b.p_min.x - s, y)),
                Bounds2::new(Point2::new(b.p_min.x - s, y), Point2::new(b.p_max.x, b.p_max.y))
            );
            prop_assert_eq!(
                b.union(&Point2::new(b.p_max.x + s, y)),
                Bounds2::new(Point2::new(b.p_min.x, y), Point2::new(b.p_max.x + s, b.p_max.y))
            );

            let y = lerp(t, b.p_min.y as f32, b.p_max.y as f32).round() as i32;
            prop_assert_eq!(
                b.union(&Point2::new(b.p_min.x - s, y)),
                Bounds2::new(Point2::new(b.p_min.x - s, b.p_min.y), Point2::new(b.p_max.x, b.p_max.y))
            );
            prop_assert_eq!(
                b.union(&Point2::new(b.p_max.x + s, y)),
                Bounds2::new(Point2::new(b.p_min.x, b.p_min.y), Point2::new(b.p_max.x + s, b.p_max.y))
            );

            let y = lerp(t + 1.0, b.p_min.y as f32, b.p_max.y as f32).round() as i32;
            prop_assert_eq!(
                b.union(&Point2::new(b.p_min.x - s, y)),
                Bounds2::new(Point2::new(b.p_min.x - s, b.p_min.y), Point2::new(b.p_max.x, y))
            );
            prop_assert_eq!(
                b.union(&Point2::new(b.p_max.x + s, y)),
                Bounds2::new(Point2::new(b.p_min.x, b.p_min.y), Point2::new(b.p_max.x + s, y))
            );

            let x = lerp(t, b.p_min.x as f32, b.p_max.x as f32).round() as i32;
            prop_assert_eq!(
                b.union(&Point2::new(x, b.p_min.y - s)),
                Bounds2::new(Point2::new(b.p_min.x, b.p_min.y - s), Point2::new(b.p_max.x, b.p_max.y))
            );
            prop_assert_eq!(
                b.union(&Point2::new(x, b.p_max.y + s)),
                Bounds2::new(Point2::new(b.p_min.x, b.p_min.y), Point2::new(b.p_max.x, b.p_max.y + s))
            );
        }

        #[test]
        fn union_non_empty_bounds2f_with_exterior_point_returns_non_empty_bounds2f(
            p in point2_f32(),
            dx in 0.001..10.0f32, dy in 0.001..10.0f32,
            s in 0.0..1.0f32, t in 0.0..1.0f32,
        ) {
            let v = Vector2::new(dx, dy);
            let b = Bounds2::new(p - v, p + v);

            let y = lerp(t - 1.0, b.p_min.y, b.p_max.y);
            prop_assert_eq!(
                b.union(&Point2::new(b.p_min.x - s, y)),
                Bounds2::new(Point2::new(b.p_min.x - s, y), Point2::new(b.p_max.x, b.p_max.y))
            );
            prop_assert_eq!(
                b.union(&Point2::new(b.p_max.x + s, y)),
                Bounds2::new(Point2::new(b.p_min.x, y), Point2::new(b.p_max.x + s, b.p_max.y))
            );

            let y = lerp(t, b.p_min.y, b.p_max.y);
            prop_assert_eq!(
                b.union(&Point2::new(b.p_min.x - s, y)),
                Bounds2::new(Point2::new(b.p_min.x - s, b.p_min.y), Point2::new(b.p_max.x, b.p_max.y))
            );
            prop_assert_eq!(
                b.union(&Point2::new(b.p_max.x + s, y)),
                Bounds2::new(Point2::new(b.p_min.x, b.p_min.y), Point2::new(b.p_max.x + s, b.p_max.y))
            );

            let y = lerp(t + 1.0, b.p_min.y, b.p_max.y);
            prop_assert_eq!(
                b.union(&Point2::new(b.p_min.x - s, y)),
                Bounds2::new(Point2::new(b.p_min.x - s, b.p_min.y), Point2::new(b.p_max.x, y))
            );
            prop_assert_eq!(
                b.union(&Point2::new(b.p_max.x + s, y)),
                Bounds2::new(Point2::new(b.p_min.x, b.p_min.y), Point2::new(b.p_max.x + s, y))
            );

            let x = lerp(t, b.p_min.x, b.p_max.x);
            prop_assert_eq!(
                b.union(&Point2::new(x, b.p_min.y - s)),
                Bounds2::new(Point2::new(b.p_min.x, b.p_min.y - s), Point2::new(b.p_max.x, b.p_max.y))
            );
            prop_assert_eq!(
                b.union(&Point2::new(x, b.p_max.y + s)),
                Bounds2::new(Point2::new(b.p_min.x, b.p_min.y), Point2::new(b.p_max.x, b.p_max.y + s))
            );
        }

        #[test]
        fn union_non_empty_bounds2i_with_interior_point_returns_same_bounds2i(
            p in point2_i32(),
            dx in 1..10i32, dy in 1..10i32,
            tx in 0.0..1.0f32, ty in 0.0..1.0f32,
        ) {
            let v = Vector2::new(dx, dy);
            let b = Bounds2::new(p - v, p + v);

            let x = lerp(tx, b.p_min.x as f32, b.p_max.x as f32).round() as i32;
            let y = lerp(ty, b.p_min.y as f32, b.p_max.y as f32).round() as i32;
            prop_assert_eq!(b.union(&Point2::new(x, y)), b);
        }

        #[test]
        fn union_non_empty_bounds2f_with_interior_point_returns_same_bounds2f(
            p in point2_f32(),
            dx in 0.001..10.0f32, dy in 0.001..10.0f32,
            tx in 0.0..1.0f32, ty in 0.0..1.0f32,
        ) {
            let v = Vector2::new(dx, dy);
            let b = Bounds2::new(p - v, p + v);

            let x = lerp(tx, b.p_min.x, b.p_max.x);
            let y = lerp(ty, b.p_min.y, b.p_max.y);
            prop_assert_eq!(b.union(&Point2::new(x, y)), b);
        }

        #[test]
        fn union_non_empty_non_overlapping_bounds2f_returns_non_empty_bounds2f(
            p in point2_f32(),
            dx in 0.001..10.0f32, dy in 0.001..10.0f32,
            s in 0.002..10.0f32, t1 in -1.0..2.0f32, t2 in -1.0..2.0f32,
        ) {
            let v = Vector2::new(dx, dy);
            let b1 = Bounds2::new(p - v, p + v);

            let x1 = lerp(t1, b1.p_min.x, b1.p_max.x);
            let x2 = lerp(t2, b1.p_min.x, b1.p_max.x);
            let y1 = lerp(t1, b1.p_min.y, b1.p_max.y);
            let y2 = lerp(t2, b1.p_min.y, b1.p_max.y);

            let tests = vec!(
                (
                    Bounds2::new(
                        Point2::new(b1.p_min.x - s, y1),
                        Point2::new(b1.p_min.x - 0.001, y2),
                    ),
                    Bounds2::new(
                        Point2::new(b1.p_min.x - s, min(b1.p_min.y, min(y1, y2))),
                        Point2::new(b1.p_max.x, max(b1.p_max.y, max(y1, y2))),
                    )
                ),
                (
                    Bounds2::new(
                        Point2::new(b1.p_max.x + 0.001, y1),
                        Point2::new(b1.p_max.x + s, y2),
                    ),
                    Bounds2::new(
                        Point2::new(b1.p_min.x, min(b1.p_min.y, min(y1, y2))),
                        Point2::new(b1.p_max.x + s, max(b1.p_max.y, max(y1, y2))),
                    ),
                ),
                (
                    Bounds2::new(
                        Point2::new(x1, b1.p_min.y - s),
                        Point2::new(x2, b1.p_min.y - 0.001),
                    ),
                    Bounds2::new(
                        Point2::new(min(b1.p_min.x, min(x1, x2)), b1.p_min.y - s),
                        Point2::new(max(b1.p_max.x, max(x1, x2)), b1.p_max.y),
                    )
                ),
                (
                    Bounds2::new(
                        Point2::new(x1, b1.p_max.y + 0.001),
                        Point2::new(x2, b1.p_max.y + s),
                    ),
                    Bounds2::new(
                        Point2::new(min(b1.p_min.x, min(x1, x2)), b1.p_min.y),
                        Point2::new(max(b1.p_max.x, max(x1, x2)), b1.p_max.y + s),
                    )
                )
            );
            for (b2, b3) in tests {
                prop_assert_eq!(b1.union(&b2), b3);
                prop_assert_eq!(b2.union(&b1), b3);
            }
        }

        #[test]
        fn union_non_empty_overlapping_bounds2f_returns_non_empty_bounds2f(
            p in point2_f32(),
            dx in 0.001..10.0f32, dy in 0.001..10.0f32,
            t1 in -2.0..1.0f32, t2 in 0.0..2.0f32,
            s1 in -2.0..1.0f32, s2 in 0.0..2.0f32,
        ) {
            let v = Vector2::new(dx, dy);
            let b1 = Bounds2::new(p - v, p + v);

            let x1 = lerp(t1, b1.p_min.x, b1.p_max.x);
            let x2 = lerp(t2, b1.p_min.x, b1.p_max.x);
            let y1 = lerp(s1, b1.p_min.y, b1.p_max.y);
            let y2 = lerp(s2, b1.p_min.y, b1.p_max.y);

            let b2 = Bounds2::new(Point2::new(x1, y1), Point2::new(x2, y2));
            let b3 = Bounds2::new(
                Point2::new(min(b1.p_min.x, min(x1, x2)), min(b1.p_min.y, min(y1, y2))),
                Point2::new(max(b1.p_max.x, max(x1, x2)), max(b1.p_max.y, max(y1, y2))),
            );

            prop_assert_eq!(b1.union(&b2), b3);
            prop_assert_eq!(b2.union(&b1), b3);
        }

        #[test]
        fn intersect_empty_with_non_empty_bounds2i_returns_empty(
            p1 in point2_i32(), p2 in point2_i32(),
        ) {
            let empty = Bounds2::<i32>::empty();
            let non_empty = Bounds2::new(p1, p2);
            prop_assert!(empty.intersect(&non_empty).is_empty());
            prop_assert!(non_empty.intersect(&empty).is_empty());
        }

        #[test]
        fn intersect_empty_with_non_empty_bounds2f_returns_empty(
            p1 in point2_f32(), p2 in point2_f32(),
        ) {
            let empty = Bounds2::<f32>::empty();
            let non_empty = Bounds2::new(p1, p2);
            prop_assert!(empty.intersect(&non_empty).is_empty());
            prop_assert!(non_empty.intersect(&empty).is_empty());
        }

        #[test]
        fn intersect_non_empty_non_overlapping_bounds2f_returns_empty(
            p in point2_f32(),
            dx in 0.001..10.0f32, dy in 0.001..10.0f32,
            s in 0.002..10.0f32, t1 in -1.0..2.0f32, t2 in -1.0..2.0f32,
        ) {
            let v = Vector2::new(dx, dy);
            let b1 = Bounds2::new(p - v, p + v);

            let x1 = lerp(t1, b1.p_min.x, b1.p_max.x);
            let x2 = lerp(t2, b1.p_min.x, b1.p_max.x);
            let y1 = lerp(t1, b1.p_min.y, b1.p_max.y);
            let y2 = lerp(t2, b1.p_min.y, b1.p_max.y);

            let tests = vec!(
                Bounds2::new(
                    Point2::new(b1.p_min.x - s, y1),
                    Point2::new(b1.p_min.x - 0.001, y2),
                ),
                Bounds2::new(
                    Point2::new(b1.p_max.x + 0.001, y1),
                    Point2::new(b1.p_max.x + s, y2),
                ),
                Bounds2::new(
                    Point2::new(x1, b1.p_min.y - s),
                    Point2::new(x2, b1.p_min.y - 0.001),
                ),
                Bounds2::new(
                    Point2::new(x1, b1.p_max.y + 0.001),
                    Point2::new(x2, b1.p_max.y + s),
                ),
            );
            for b2 in tests {
                prop_assert!(b1.intersect(&b2).is_empty());
                prop_assert!(b2.intersect(&b1).is_empty());
            }
        }

        #[test]
        fn intersect_non_empty_overlapping_bounds2f_returns_non_empty(
            p in point2_f32(),
            dx in 0.001..10.0f32, dy in 0.001..10.0f32,
            t1 in -2.0..1.0f32, t2 in 0.0..2.0f32,
            s1 in -2.0..1.0f32, s2 in 0.0..2.0f32,
        ) {
            let v = Vector2::new(dx, dy);
            let b1 = Bounds2::new(p - v, p + v);

            let x1 = lerp(t1, b1.p_min.x, b1.p_max.x);
            let x2 = lerp(t2, b1.p_min.x, b1.p_max.x);
            let y1 = lerp(s1, b1.p_min.y, b1.p_max.y);
            let y2 = lerp(s2, b1.p_min.y, b1.p_max.y);

            let b2 = Bounds2::new(Point2::new(x1, y1), Point2::new(x2, y2));
            let b3 = Bounds2::new(
                Point2::new(max(b1.p_min.x, min(x1, x2)), max(b1.p_min.y, min(y1, y2))),
                Point2::new(min(b1.p_max.x, max(x1, x2)), min(b1.p_max.y, max(y1, y2))),
            );

            prop_assert_eq!(b1.intersect(&b2), b3);
            prop_assert_eq!(b2.intersect(&b1), b3);
        }

        #[test]
        fn iterate_bounds2i_returns_grid_points_left_to_right_x_and_top_to_bottom_y(
            p in point2_i32(), dx in 1..10i32, dy in 1..10i32,
        ) {
            let b = Bounds2::new(p, p + Vector2::new(dx, dy));
            let mut iter = b.into_iter();

            for y in 0..dy {
                for x in 0..dx {
                    prop_assert_eq!(iter.next(), Some(Point2::new(p.x + x, p.y + y)));
                }
            }

            prop_assert!(iter.next().is_none());
        }

        #[test]
        fn iterate_bounds2i_with_0_in_one_dimension_returns_grid_points_along_the_other(
            p in point2_i32(), d in 1..10i32,
        ) {
            let b1 = Bounds2::new(p, p + Vector2::new(0, d));
            let b2 = Bounds2::new(p, p + Vector2::new(d, 0));

            let mut iter1 = b1.into_iter();
            let mut iter2 = b2.into_iter();

            for i in 0..d {
                prop_assert_eq!(iter1.next(), Some(Point2::new(p.x, p.y + i)));
                prop_assert_eq!(iter2.next(), Some(Point2::new(p.x + i, p.y)));
            }

            prop_assert!(iter1.next().is_none());
            prop_assert!(iter2.next().is_none());
        }
    }
}
