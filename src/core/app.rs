//! Application related stuff

#![allow(dead_code)]
use crate::core::pbrt::Float;
use clap::*;

lazy_static! {
    /// The global application options.
    pub static ref OPTIONS: Options = Options::new();
}

/// System wide options.
#[derive(Clone, Debug)]
pub struct Options {
    /// Number of threads to use for rendering.
    pub n_threads: usize,

    /// Automatically reduce a number of quality settings to render more quickly.
    pub quick_render: bool,

    /// Suppress all text output other than error messages.:
    pub quiet: bool,

    /// Print a reformatted version of the input file(s) to the standard output
    /// instead of rendering an image.
    pub cat: bool,

    /// Print a reformatted version of the input file(s) to standard output and
    /// convert all triangle meshes to PLY files. Does not render an image.
    pub to_ply: bool,

    /// Path to the image file.
    pub image_file: String,

    /// The crop window x0, x1, y0, y1.
    pub crop_window: [[Float; 2]; 2],
}

impl Options {
    /// Loads the command line options.
    pub fn new() -> Self {
        let matches = app_from_crate!()
            .arg(
                Arg::with_name("nthreads")
                    .short("t")
                    .long("nthreads")
                    .value_name("NUM")
                    .default_value("1")
                    .takes_value(true)
                    .help("Use specified number of threads for rendering."),
            )
            .arg(
                Arg::with_name("outfile")
                    .short("o")
                    .long("outfile")
                    .value_name("FILE")
                    .takes_value(true)
                    .help("Write the final image to the given filename."),
            )
            .arg(
                Arg::with_name("cropwindow")
                    .short("cw")
                    .long("cropwindow")
                    .value_name("x0 y0 x1 y1")
                    .number_of_values(4)
                    .takes_value(true)
                    .help("Specify an image crop window."),
            )
            .arg(
                Arg::with_name("quick")
                    .long("quick")
                    .takes_value(false)
                    .default_value("false")
                    .help(
                        "Automatically reduce a number of quality settings to \
                        render more quickly.",
                    ),
            )
            .arg(
                Arg::with_name("quiet")
                    .long("quiet")
                    .takes_value(false)
                    .default_value("false")
                    .help("Suppress all text output other than error messages."),
            )
            .arg(
                Arg::with_name("minloglevel")
                    .long("minloglevel")
                    .takes_value(true)
                    .default_value("0")
                    .help(
                        "Log messages at or above this level (0 -> INFO, 1 -> \
                        WARNING, 2 -> ERROR, 3 -> FATAL).",
                    ),
            )
            .arg(
                Arg::with_name("cat")
                    .long("cat")
                    .takes_value(false)
                    .default_value("false")
                    .help(
                        "Print a reformatted version of the input file(s) to \
                    standard output. Does not render an image.",
                    ),
            )
            .arg(
                Arg::with_name("toply")
                    .long("toply")
                    .takes_value(false)
                    .default_value("false")
                    .help(
                        "Print a reformatted version of the input file(s) to \
                    standard output and convert all triangle meshes to PLY files. \
                    Does not render an image.",
                    ),
            )
            .get_matches();

        let max_threads = num_cpus::get();
        let n_threads = match matches.value_of("nthreads") {
            Some(s) => {
                let n = s.parse::<usize>().expect("Invalid nthreads");

                if n == 0 {
                    panic!("Invalid nthreads");
                } else if n > max_threads {
                    panic!(format!("Num threads > max logical CPUs {}", max_threads));
                }

                n
            }

            _ => 1,
        };

        let image_file = match matches.value_of("outfile") {
            Some(s) => s.to_string(),
            _ => panic!("Missing outfile"),
        };

        let crop_window = match matches.values_of("cropwindow") {
            Some(s) => {
                let v: Vec<&str> = s.collect();
                [
                    [
                        v[0].parse::<Float>().expect("Invalid cropwindow.x0"),
                        v[1].parse::<Float>().expect("Invalid cropwindow.y0"),
                    ],
                    [
                        v[2].parse::<Float>().expect("Invalid cropwindow.x1"),
                        v[3].parse::<Float>().expect("Invalid cropwindow.y1"),
                    ],
                ]
            }
            _ => [[0.0, 1.0], [0.0, 1.0]],
        };

        let quick_render = match matches.value_of("quick") {
            Some(s) => s.parse::<bool>().expect("Invalid quick"),
            _ => false,
        };

        let quiet = match matches.value_of("quiet") {
            Some(s) => s.parse::<bool>().expect("Invalid quiet"),
            _ => false,
        };

        let cat = match matches.value_of("cat") {
            Some(s) => s.parse::<bool>().expect("Invalid cat"),
            _ => false,
        };

        let to_ply = match matches.value_of("toply") {
            Some(s) => s.parse::<bool>().expect("Invalid toply"),
            _ => false,
        };

        Self {
            n_threads,
            quick_render,
            quiet,
            cat,
            to_ply,
            image_file,
            crop_window,
        }
    }
}
