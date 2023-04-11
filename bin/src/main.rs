#[macro_use]
extern crate log;

use api::*;
use core::app::*;
use core::fileutil::*;

#[cfg(all(feature = "dhat-rs", feature = "jemalloc"))]
compile_error!("feature 'dhat-rs' and feature 'jemalloc' cannot be enabled at the same time");

#[cfg(feature = "dhat-rs")]
use dhat::{Dhat, DhatAlloc};

#[cfg(feature = "dhat-rs")]
#[global_allocator]
static ALLOCATOR: DhatAlloc = DhatAlloc;

#[cfg(feature = "jemalloc")]
#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(feature = "jemalloc")]
#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static ALLOCATOR: Jemalloc = Jemalloc;

fn main() {
    #[cfg(feature = "dhat-rs")]
    let _dhat = Dhat::start_heap_profiling();

    // Initialize `env_logger`.
    env_logger::init();

    // Initialize PBRT API.
    let mut api = Api::new();
    api.pbrt_init();

    // Process scene description.
    OPTIONS.paths.iter().for_each(|path| {
        match absolute_path(path).and_then(|abs_path| parser::parse(&abs_path, &mut api)) {
            Err(e) => error!("{}", e),
            _ => (),
        }
    });

    api.pbrt_cleanup();
}
