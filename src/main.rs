#![feature(iter_partition_in_place)] // Can be removed after Rust 1.49 release.

#[macro_use]
extern crate lazy_static;

extern crate exr;
extern crate float_cmp;
extern crate hexf;
extern crate image;
extern crate num_traits;
extern crate rand;
extern crate rand_pcg;
extern crate typed_arena;

mod accelerators;
mod cameras;
mod core;
mod filters;
mod lights;
mod samplers;
mod shapes;

use self::core::geometry::*;

fn main() {
    let mut v = Vector2::<f32>::zero();
    v += Vector2f { x: 0.2, y: 0.1 };
    let p = Point2f { x: 0.5, y: 0.6 };

    println!("{:?} - {:?} = {:?}!", p, v, p - v);
    println!("{:?} + {:?} = {:?}!", p, v, p + v);
}
