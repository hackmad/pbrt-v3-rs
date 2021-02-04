#![feature(iter_partition_in_place)] // Can be removed once stable

extern crate byteorder;
extern crate clap;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate exr;
extern crate float_cmp;
#[macro_use]
extern crate hexf;
extern crate image;
#[macro_use]
extern crate lazy_static;
extern crate nom;
extern crate num_traits;
extern crate ordered_float;
extern crate rand;
extern crate rand_pcg;
extern crate typed_arena;

mod accelerators;
mod cameras;
mod core;
mod filters;
mod lights;
mod materials;
mod samplers;
mod shapes;
mod textures;

fn main() {
    env_logger::init();
    info!("Hello, World!");
}
