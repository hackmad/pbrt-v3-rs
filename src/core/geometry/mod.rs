//! Geometry
use super::pbrt::*;

// Define macros for property based testing.
#[cfg(test)]
#[macro_export]
macro_rules! prop_range {
    ($name: ident, $t: ty, $r: expr) => {
        prop_compose! {
            fn $name()(f in $r) -> $t {
                f
            }
        }
    };
}

#[cfg(test)]
#[macro_export]
macro_rules! prop_non_zero_range {
    ($name: ident, $t: ty, $r: expr) => {
        prop_compose! {
            fn $name()(f in $r.prop_filter("non-zero", |x| !(*x).is_zero())) -> $t {
                f
            }
        }
    };
}

#[cfg(test)]
#[macro_export]
macro_rules! prop_vector2 {
    ($name: ident, $t: ty, $xr: expr, $yr: expr) => {
        prop_compose! {
            fn $name()(x in $xr, y in $yr) -> Vector2<$t> {
                Vector2 { x, y }
            }
        }
    };
}

#[cfg(test)]
#[macro_export]
macro_rules! prop_vector3 {
    ($name: ident, $t: ty, $xr: expr, $yr: expr, $zr: expr) => {
        prop_compose! {
            fn $name()(x in $xr, y in $yr, z in $zr) -> Vector3<$t> {
                Vector3 { x, y, z }
            }
        }
    };
}

#[cfg(test)]
#[macro_export]
macro_rules! prop_normal3 {
    ($name: ident, $t: ty, $xr: expr, $yr: expr, $zr: expr) => {
        prop_compose! {
            fn $name()(x in $xr, y in $yr, z in $zr) -> Normal3<$t> {
                Normal3 { x, y, z }
            }
        }
    };
}

#[cfg(test)]
#[macro_export]
macro_rules! prop_point2 {
    ($name: ident, $t: ty, $xr: expr, $yr: expr) => {
        prop_compose! {
            fn $name()(x in $xr, y in $yr) -> Point2<$t> {
                Point2 { x, y }
            }
        }
    };
}

#[cfg(test)]
#[macro_export]
macro_rules! prop_point3 {
    ($name: ident, $t: ty, $xr: expr, $yr: expr, $zr: expr) => {
        prop_compose! {
            fn $name()(x in $xr, y in $yr, z in $zr) -> Point3<$t> {
                Point3 { x, y, z }
            }
        }
    };
}

mod common;
mod coordinate_system;
mod normal;
mod point2;
mod point3;
mod ray;
mod vector2;
mod vector3;

// Re-export
pub use common::*;
pub use coordinate_system::*;
pub use normal::*;
pub use point2::*;
pub use point3::*;
pub use ray::*;
pub use vector2::*;
pub use vector3::*;
