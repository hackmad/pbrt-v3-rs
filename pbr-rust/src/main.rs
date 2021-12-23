#[macro_use]
extern crate log;

use api::parser::*;
use api::*;
use core::app::*;

#[cfg(feature = "dhat-rs")]
use dhat::{Dhat, DhatAlloc};

#[cfg(feature = "dhat-rs")]
#[global_allocator]
static ALLOCATOR: DhatAlloc = DhatAlloc;

fn main() {
    #[cfg(feature = "dhat-rs")]
    let _dhat = Dhat::start_heap_profiling();

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
