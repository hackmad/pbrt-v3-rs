//! Image I/O

#![allow(dead_code)]
use super::geometry::*;
use super::pbrt::*;
use super::spectrum::*;
use exr::prelude as exrs;
use exr::prelude::*;
use image::*;
use std::ffi::OsStr;
use std::path::Path;

/// Stores RGB image data.
pub struct RGBImage {
    /// The pixels.
    pixels: Vec<RGBSpectrum>,

    /// Image resolution.
    resolution: Point2i,
}

/// Read an image.
///
/// * `path` - Input file path.
pub fn read_image(path: &str) -> RGBImage {
    match get_extension_from_filename(path) {
        Some(".exr") => read_exr(path),
        Some(_extension) => read_8_bit(path),
        None => panic!("Can't determine file type from suffix of filename {}", path),
    }
}

/// Read a single layer OpenEXR file.
///
/// * `path` - Input file path.
fn read_exr(path: &str) -> RGBImage {
    let reader = exrs::read()
        .no_deep_data()
        .largest_resolution_level()
        .rgba_channels(
            |layer_info: &exrs::RgbaChannelsInfo| {
                let width = layer_info.resolution.width();
                let height = layer_info.resolution.height();
                RGBImage {
                    pixels: vec![RGBSpectrum::default(); width * height],
                    resolution: Point2i::new(width as Int, height as Int),
                }
            },
            |img: &mut RGBImage, position: exrs::Vec2<usize>, pixel: exrs::RgbaPixel| {
                let offset = position.y() & (img.resolution.x as usize) + position.x();
                img.pixels[offset] = RGBSpectrum::from(vec![
                    pixel.red.to_f32(),
                    pixel.green.to_f32(),
                    pixel.blue.to_f32(),
                ]);
            },
        )
        .first_valid_layer()
        .all_attributes();

    // Return the `RGBImage`.
    let image: Image<Layer<RgbaChannels<RGBImage>>> = reader
        .from_file(path)
        .expect(format!("Unable to read file {}", path).as_ref());
    image.layer_data.channel_data.storage
}

/// Read an 8-bit image format.
///
/// * `path` - Input file path.
fn read_8_bit(path: &str) -> RGBImage {
    // Read image and convert to RGB.
    let img: RgbImage = open(path)
        .expect(format!("Unable to open {}", path).as_ref())
        .into_rgb8();

    // Read metadata.
    let width = img.width() as usize;
    let height = img.height() as usize;
    let resolution = Point2i::new(width as Int, height as Int);

    // Iterate over the coordinates and pixels of the image
    let pixels: Vec<RGBSpectrum> = img
        .pixels()
        .map(|rgb_u8| {
            RGBSpectrum::from(vec![
                rgb_u8[0] as Float / 255.0,
                rgb_u8[1] as Float / 255.0,
                rgb_u8[2] as Float / 255.0,
            ])
        })
        .collect();

    // Return the `RGBImage`.
    RGBImage { pixels, resolution }
}

/// Write the output image to given path.
///
/// * `path`             - Output file path.
/// * `rgb`              - Floating point RGB pixel data.
/// * `output_bounds`    - The bounds for the image output.
pub fn write_image(path: &str, rgb: &[Float], output_bounds: &Bounds2i) {
    let resolution = output_bounds.diagonal();
    let res_x = resolution.x as u32;
    let res_y = resolution.y as u32;

    match get_extension_from_filename(path) {
        Some(".exr") => write_exr(path, rgb, res_x, res_y),
        Some(".tga") => write_8_bit(path, rgb, res_x, res_y, ImageFormat::Tga),
        Some(".png") => write_8_bit(path, rgb, res_x, res_y, ImageFormat::Png),
        Some(extension) => eprintln!("Extension {} is not supported", extension),
        None => eprintln!("Can't determine file type from suffix of filename {}", path),
    }
}

/// Retrieve the extension from a file path.
///
/// * `path` - The file path.
fn get_extension_from_filename(path: &str) -> Option<&str> {
    Path::new(path).extension().and_then(OsStr::to_str)
}

/// Writes the image in OpenEXR format.
///
/// * `path`        - Output file path.
/// * `rgb`         - Floating point RGB pixel data.
/// * `res_x`       - X resolution.
/// * `res_y`       - Y resolution.
fn write_exr(path: &str, rgb: &[Float], res_x: u32, res_y: u32) {
    println!("Writing image {} with resolution {}x{}", path, res_x, res_y);

    write_rgb_f32_file(
        String::from(path),
        (res_x as usize, res_y as usize),
        |x, y| {
            let offset = y * (res_x as usize) + x;
            (rgb[offset], rgb[offset + 1], rgb[offset + 2])
        },
    )
    .expect("Error saving output file");
}

/// Writes the image in an 8-bit image format.
///
/// * `path`         - Output file path.
/// * `rgb`          - Floating point RGB pixel data.
/// * `res_x`        - X resolution.
/// * `res_y`        - Y resolution.
/// * `image_format` - Image format.
fn write_8_bit(path: &str, rgb: &[Float], res_x: u32, res_y: u32, image_format: ImageFormat) {
    println!("Writing image {} with resolution {}x{}", path, res_x, res_y);

    // Allocate an image buffer.
    let mut imgbuf = ImageBuffer::new(res_x, res_y);
    let mut offset = 0;
    for y in 0..res_y {
        for x in 0..res_x {
            // 8-bit format; apply gamma and clamp.
            let rgb = apply_gamma(&[rgb[offset], rgb[offset + 1], rgb[offset + 2]]);
            imgbuf.put_pixel(x, res_y - 1 - y, Rgb(rgb));
            offset += 3;
        }
    }

    // Write the output file.
    imgbuf
        .save_with_format(String::from(path), image_format)
        .expect("Error saving output file");
}

/// Apply gamma correction to a RGB floating point pixel and return the
/// clamped 8-bit values.
///
/// * `rgb` - RGB floating point pixel value.
#[inline]
fn apply_gamma(rgb: &[Float; 3]) -> [u8; 3] {
    [clamp_byte(rgb[0]), clamp_byte(rgb[1]), clamp_byte(rgb[2])]
}

/// Clamp floating point value to 8-bit range [0, 255].
///
/// * `v` - Value to clamp.
#[inline]
fn clamp_byte(v: Float) -> u8 {
    clamp(255.0 * gamma_correct(v) + 0.5, 0.0, 255.0) as u8
}
