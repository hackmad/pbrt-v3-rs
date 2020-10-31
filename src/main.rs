#![feature(iter_partition_in_place)] // Can be removed after Rust 1.49 release.

extern crate float_cmp;
extern crate num_traits;
extern crate typed_arena;

mod accelerators;
mod core;
mod shapes;

use self::core::geometry::*;

fn main() {
    let mut v = Vector2::<f32>::zero();
    v += Vector2f { x: 0.2, y: 0.1 };
    let p = Point2f { x: 0.5, y: 0.6 };

    println!("{:?} - {:?} = {:?}!", p, v, p - v);
    println!("{:?} + {:?} = {:?}!", p, v, p + v);
}
