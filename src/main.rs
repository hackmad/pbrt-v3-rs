extern crate float_cmp;
extern crate num_traits;

mod core;

use self::core::geometry::*;

fn main() {
    let mut v = zero_vector2();
    v += Vector2f { x: 0.2, y: 0.1 };
    let p = Point2f { x: 0.5, y: 0.6 };

    println!("{:?} - {:?} = {:?}!", p, v, p - v);
    println!("{:?} + {:?} = {:?}!", p, v, p + v);
}
