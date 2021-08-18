//! Image I/O

#![allow(dead_code)]
use crate::core::geometry::*;
use crate::core::pbrt::*;
use crate::core::spectrum::*;
use exr::prelude as exrs;
use exr::prelude::*;
use image::*;
use regex::Regex;
use std::result::Result;

/// Stores RGB image data.
pub struct RGBImage {
    /// The pixels.
    pub pixels: Vec<RGBSpectrum>,

    /// Image resolution.
    pub resolution: Point2<usize>,
}

/// Read an image.
///
/// * `path` - Input file path.
pub fn read_image(path: &str) -> Result<RGBImage, String> {
    match get_extension_from_filename(path) {
        Some(".exr") => read_exr(path),
        Some(_extension) => read_8_bit(path),
        None => Err(format!(
            "Can't determine file type from suffix of filename {}.",
            path
        )),
    }
}

/// Read a single layer OpenEXR file.
///
/// * `path` - Input file path.
fn read_exr(path: &str) -> Result<RGBImage, String> {
    let reader = exrs::read()
        .no_deep_data()
        .largest_resolution_level()
        .rgba_channels(
            |resolution, _channels| {
                let width = resolution.width();
                let height = resolution.height();
                RGBImage {
                    pixels: vec![RGBSpectrum::default(); width * height],
                    resolution: Point2::new(width, height),
                }
            },
            |img, position, (r, g, b, _a): (f32, f32, f32, f32)| {
                let offset = position.y() * img.resolution.x + position.x();
                img.pixels[offset] = RGBSpectrum::from(vec![r, g, b]);
            },
        )
        .first_valid_layer()
        .all_attributes();

    // Return the `RGBImage`.
    match reader.from_file(path) {
        Ok(image) => Ok(image.layer_data.channel_data.pixels),
        Err(err) => Err(format!("{:}", err)),
    }
}

/// Read an 8-bit image format.
///
/// * `path` - Input file path.
fn read_8_bit(path: &str) -> Result<RGBImage, String> {
    // Read image and convert to RGB.
    let img: RgbImage = match open(path) {
        Ok(i) => i.into_rgb8(),
        Err(err) => return Err(format!("{:}", err)),
    };

    // Read metadata.
    let width = img.width() as usize;
    let height = img.height() as usize;
    let resolution = Point2::new(width, height);

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
    Ok(RGBImage { pixels, resolution })
}

/// Write the output image to given path.
///
/// * `path`             - Output file path.
/// * `rgb`              - Floating point RGB pixel data.
/// * `output_bounds`    - The bounds for the image output.
pub fn write_image(path: &str, rgb: &[Float], output_bounds: &Bounds2i) -> Result<(), String> {
    let resolution = output_bounds.diagonal();
    let res_x = resolution.x as u32;
    let res_y = resolution.y as u32;

    match get_extension_from_filename(path) {
        Some(".exr") => write_exr(path, rgb, res_x, res_y),
        Some(".tga") => write_8_bit(path, rgb, res_x, res_y, ImageFormat::Tga),
        Some(".png") => write_8_bit(path, rgb, res_x, res_y, ImageFormat::Png),
        Some(extension) => Err(format!("Extension {} is not supported", extension)),
        None => Err(format!(
            "Can't determine file type from suffix of filename {}",
            path
        )),
    }
}

lazy_static! {
    /// Regular expression for extracting the file extension. This will
    /// match the last occurrence of a period followed by no periods or
    /// slashes (could be tightened to exclude other illegal characters
    /// but the code that reads files will bomb anyway).
    static ref REGEX_FILE_EXT: Regex = Regex::new(r"(\.[^./\\]+)$").unwrap();
}

/// Retrieve the extension from a file path.
///
/// * `path` - The file path.
fn get_extension_from_filename(path: &str) -> Option<&str> {
    REGEX_FILE_EXT
        .captures(path)
        .map(|c| c.get(1).map_or("", |m| m.as_str()))
}

/// Writes the image in OpenEXR format.
///
/// * `path`        - Output file path.
/// * `rgb`         - Floating point RGB pixel data.
/// * `res_x`       - X resolution.
/// * `res_y`       - Y resolution.
fn write_exr(path: &str, rgb: &[Float], res_x: u32, res_y: u32) -> Result<(), String> {
    info!("Writing image {} with resolution {}x{}", path, res_x, res_y);

    match write_rgb_file(
        String::from(path),
        res_x as usize,
        res_y as usize,
        |x, y| {
            let offset = 3 * (y * (res_x as usize) + x);
            (rgb[offset], rgb[offset + 1], rgb[offset + 2])
        },
    ) {
        Ok(()) => Ok(()),
        Err(err) => Err(format!("Error saving output image {}. {:}.", path, err)),
    }
}

/// Writes the image in an 8-bit image format.
///
/// * `path`         - Output file path.
/// * `rgb`          - Floating point RGB pixel data.
/// * `res_x`        - X resolution.
/// * `res_y`        - Y resolution.
/// * `image_format` - Image format.
fn write_8_bit(
    path: &str,
    rgb: &[Float],
    res_x: u32,
    res_y: u32,
    image_format: ImageFormat,
) -> std::result::Result<(), String> {
    info!("Writing image {} with resolution {}x{}", path, res_x, res_y);

    // Allocate an image buffer.
    let mut imgbuf = ImageBuffer::new(res_x, res_y);
    let mut offset = 0;
    for y in 0..res_y {
        for x in 0..res_x {
            // 8-bit format; apply gamma and clamp.
            let rgb = apply_gamma(&[rgb[offset], rgb[offset + 1], rgb[offset + 2]]);
            imgbuf.put_pixel(x, y, Rgb(rgb));
            offset += 3;
        }
    }

    // Write the output file.
    match imgbuf.save_with_format(String::from(path), image_format) {
        Ok(()) => Ok(()),
        Err(err) => Err(format!("Error saving output image {}. {:}.", path, err)),
    }
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
