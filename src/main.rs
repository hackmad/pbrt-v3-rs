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
extern crate itertools;
#[macro_use]
extern crate lazy_static;
extern crate num_traits;
extern crate ordered_float;
extern crate pest;
#[macro_use]
extern crate pest_derive;
extern crate rand;
extern crate rand_pcg;

mod accelerators;
mod cameras;
mod core;
mod filters;
mod integrators;
mod lights;
mod materials;
mod samplers;
mod shapes;
mod textures;

use crate::core::api::*;
use crate::core::app::*;
use crate::core::parsers::*;

fn main() {
    // Initialize `env_logger`.
    env_logger::init();

    // Load the program options.
    let options = OPTIONS.clone();

    // Configure number of threads.
    rayon::ThreadPoolBuilder::new()
        .num_threads(options.n_threads)
        .build_global()
        .unwrap();

    // Initialize PBRT API.
    let mut api = Api::new();
    api.pbrt_init();

    // Process scene description.
    for path in options.paths.iter() {
        let parser = PbrtFileParser::new(path);
        match parser.parse(&mut api) {
            Ok(_) => (),
            Err(err) => error!("{}", err),
        }
    }

    api.pbrt_cleanup();
}
