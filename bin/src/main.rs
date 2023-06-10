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
    for path in options().paths.iter() {
        // In case of error report it and continue.
        if let Err(e) = render(path, &mut api) {
            error!("{e}");
        }
    }

    api.pbrt_cleanup();
}

fn render(path: &str, api: &mut Api) -> Result<(), String> {
    // Get absolute path to the scene file.
    let abs_path = absolute_path(path)?;

    // Get path of scene file.
    let scene_path = parent_path(&abs_path).ok_or(format!("Invalid '{path}' for rendering."))?;
    api.set_current_working_dir(&scene_path);

    // Parse and render scene.
    parser::parse(&abs_path, "", api)
}
