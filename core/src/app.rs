//! Application related stuff

use crate::pbrt::Float;
use clap::Parser;

lazy_static! {
    /// The global application options.
    pub static ref OPTIONS: Options = Options::parse();
}

/// System wide options.
#[derive(Parser, Clone)]
#[command(author, version, about, long_about = None)]
pub struct Options {
    /// Number of threads to use for rendering.
    #[arg(
        long = "nthreads",
        short = 't',
        value_name = "NUM",
        default_value_t = 1,
        help = "Use specified number of threads for rendering."
    )]
    n_threads: usize,

    /// Automatically reduce a number of quality settings to render more quickly.
    #[arg(
        long = "quick",
        help = "Automatically reduce a number of quality settings to render more quickly."
    )]
    pub quick_render: bool,

    /// Suppress all text output other than error messages.:
    #[arg(long, help = "Suppress all text output other than error messages.")]
    pub quiet: bool,

    /// Path to the image file.
    #[arg(
        long = "outfile",
        short = 'o',
        value_name = "FILE",
        help = "Write the final image to the given filename."
    )]
    pub image_file: Option<String>,

    /// The crop window x0, x1, y0, y1.
    #[arg(
        long = "cropwindow",
        short = 'c',
        value_name = "FLOAT",
        number_of_values = 4,
        help = "Specify an image crop window (x0 x1 y0 y1)."
    )]
    pub crop_window: Vec<Float>,

    /// Input file paths. Empty vector implies read from stdin.
    #[arg(help = "Input files")]
    pub paths: Vec<String>,

    /// Tile size.
    #[arg(
        long = "tilesize",
        short = 'p',
        value_name = "NUM",
        default_value_t = 16,
        help = "Size in pixels of square tiles rendered per thread."
    )]
    pub tile_size: usize,

    /// Write SPPM radius image
    #[arg(long = "sppm-radius", help = "SPPM radius image file path")]
    pub sppm_radius: Option<String>,

    /// Path prefix to the mipmap file for debugging purposses.
    #[arg(
        long = "mipmap",
        short = 'm',
        value_name = "FILE",
        help = "Write the mipmap images to the given file path prefix."
    )]
    pub mipmap_file_prefix: Option<String>,
}

impl Options {
    /// Returns the number of threads to use.
    pub fn threads(&self) -> usize {
        let max_threads = std::thread::available_parallelism().map_or(1, |n| n.get());
        match self.n_threads {
            0 => {
                warn!("Invalid nthreads");
                1
            }
            n if n > max_threads => {
                warn!("Num threads > max logical CPUs {}", max_threads);
                max_threads
            }
            n => n,
        }
    }
}
