//! Image I/O

use crate::geometry::*;
use crate::pbrt::*;
use crate::spectrum::*;
use byteorder::{BigEndian, LittleEndian, ReadBytesExt, WriteBytesExt};
use exr::prelude as exrs;
use exr::prelude::*;
use image::*;
use regex::Regex;
use std::fs::File;
use std::result::Result;
use std::sync::OnceLock;

/// Stores RGB image data.
pub struct RGBImage {
    /// The pixels.
    pub pixels: Vec<RGBSpectrum>,

    /// Image resolution.
    pub resolution: Point2<usize>,
}

impl RGBImage {
    /// Creates a new `RGBImage` from pixel data.
    ///
    /// * `pixels` - RGB pixel data.
    /// * `width`  - Width of image.
    /// * `height` - Height of image.
    pub fn new(pixels: Vec<RGBSpectrum>, width: usize, height: usize) -> Self {
        assert_eq!(width * height, pixels.len());
        Self {
            pixels,
            resolution: Point2::new(width, height),
        }
    }
}

/// Read an image.
///
/// * `path` - Input file path.
pub fn read_image(path: &str) -> Result<RGBImage, String> {
    match get_extension_from_filename(path) {
        Some(".exr") => read_exr(path),
        Some(".pfm") => read_pfm(path),
        Some(_extension) => read_8_bit(path),
        None => Err(format!("Can't determine file type from suffix of filename {path}.")),
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
        Ok(image) => {
            let pixels = image.layer_data.channel_data.pixels;
            info!(
                "Read EXR image {path} ({} x {})",
                pixels.resolution.x, pixels.resolution.y,
            );
            Ok(pixels)
        }
        Err(err) => Err(err.to_string()),
    }
}

/// Returns true if the character is a space, newline or tab.
///
/// * `c` - The character to check.
#[inline(always)]
fn is_white_space(c: char) -> bool {
    c == ' ' || c == '\n' || c == '\t'
}

// Reads a "word" from the file;  i.e. it keeps reading until whitespace is reached or maximum length is reached.
// Returns the string read *not* including the whitespace if successful; otherwise an error is returned.
//
// * `file` - File to read.
// * `len`  - Maximum number of bytes to read.
fn read_pfm_word(file: &mut File, len: usize) -> Result<String, String> {
    if len == 0 {
        return Err("read_pfm_word(): invalid length 0.".to_string());
    }

    let mut s = String::new();

    let mut c = file.read_u8().map_err(|e| format!("read_pfm_word(): {e}"))?;
    while !is_white_space(c as char) && s.len() < len {
        s.push(c as char);
        c = file.read_u8().map_err(|e| format!("read_pfm_word(): {e}"))?;
    }

    if s.len() <= len {
        Ok(s)
    } else {
        Err("read_pfm_word(): filled buffer before finding whitespace.".to_string())
    }
}

/// Read a PFM (Portable FloatMap) file.
///
/// * `path` - Input file path.
fn read_pfm(path: &str) -> Result<RGBImage, String> {
    let mut file = File::open(path).map_err(|e| format!("read_pfm(): Error reading PFM file '{path}': {e}"))?;

    // Read either "Pf" or "PF".
    let ty = read_pfm_word(&mut file, 2)?;

    let n_channels = match ty.as_str() {
        "Pf" => Ok(1),
        "PF" => Ok(3),
        s => Err(format!("read_pfm(): Invalid PFM type '{s}'")),
    }?;

    // Read the rest of the header.
    // Read width.
    let width = read_pfm_word(&mut file, 80)?
        .parse::<usize>()
        .map_err(|e| format!("Error parsing PFM width: {e}"))?;

    // Read height.
    let height = read_pfm_word(&mut file, 80)?
        .parse::<usize>()
        .map_err(|e| format!("Error parsing PFM height: {e}"))?;

    // Read scale.
    let mut scale = read_pfm_word(&mut file, 80)?
        .parse::<f32>()
        .map_err(|e| format!("Error parsing PFM endian: {e}"))?;

    let file_little_endian = scale < 0.0;

    scale = scale.abs();

    // Read the data.
    let n_floats = n_channels * width * height;
    let mut data = vec![0.0_f32; n_floats];
    for y in (0..=height - 1).rev() {
        let i = y * width * n_channels;

        for j in 0..width * n_channels {
            // Apply endian conversian and endian if appropriate.
            let f = if file_little_endian {
                file.read_f32::<LittleEndian>()
            } else {
                file.read_f32::<BigEndian>()
            }
            .map_err(|e| format!("Error reading PFM pixel data y={y}, i={i}, j={j}: {e}"))?;
            data[i + j] = f * scale;
        }
    }

    // Create RGBs.
    let rgb: Vec<RGBSpectrum> = if n_channels == 1 {
        (0..width * height).map(|i| RGBSpectrum::from(data[i])).collect()
    } else {
        (0..width * height)
            .map(|i| RGBSpectrum::from_rgb(&[data[3 * i], data[3 * i + 1], data[3 * i + 2]], None))
            .collect()
    };

    info!("Read PFM image {path} ({width} x {height} x {n_channels})");

    Ok(RGBImage::new(rgb, width, height))
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

    // Iterate over the coordinates and pixels of the image.
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

    info!("Read 8-bit image {path} ({width} x {height})");

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
        Some(".pfm") => write_pfm(path, rgb, res_x, res_y),
        Some(extension) => Err(format!("Extension {extension} is not supported")),
        None => Err(format!("Can't determine file type from suffix of filename {path}")),
    }
}

/// Returns regular expression for extracting the file extension. This will match the last occurrence of a period
/// followed by no periods or slashes (could be tightened to exclude other illegal characters but the code that reads
/// files will bomb anyway).
fn regex_file_ext() -> &'static Regex {
    static DATA: OnceLock<Regex> = OnceLock::new();
    DATA.get_or_init(|| Regex::new(r"(\.[^./\\]+)$").unwrap())
}

/// Retrieve the extension from a file path.
///
/// * `path` - The file path.
fn get_extension_from_filename(path: &str) -> Option<&str> {
    regex_file_ext()
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

    let size = Vec2(res_x as usize, res_y as usize);

    let layer1 = Layer::new(
        size,
        LayerAttributes::named("render"),
        Encoding::SMALL_LOSSLESS,
        SpecificChannels::rgb(|pos: Vec2<usize>| {
            let offset = 3 * (pos.1 * (res_x as usize) + pos.0);
            (rgb[offset], rgb[offset + 1], rgb[offset + 2])
        }),
    );

    let attributes = ImageAttributes::new(IntegerBounds::from_dimensions(size));
    match Image::empty(attributes).with_layer(layer1).write().to_file(path) {
        Ok(()) => Ok(()),
        Err(err) => Err(format!("Error saving output image {path}: {err}")),
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
    info!("Writing image {path} with resolution {res_x}x{res_y}");

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
        Err(err) => Err(format!("Error saving output image {path}: {err}.")),
    }
}

/// Writes a string characters as bytes to a PFM file.
///
// * `file` - File to write.
// * `s`    - String to write.
fn write_pfm_string(file: &mut File, s: &str) -> Result<(), String> {
    for c in s.chars() {
        file.write_u8(c as u8).map_err(|e| e.to_string())?;
    }
    Ok(())
}

/// Writes the image in PFM (Portable FloatMap) format.
///
/// * `path`  - Output file path.
/// * `rgb`   - Floating point RGB pixel data.
/// * `res_x` - X resolution.
/// * `res_y` - Y resolution.
fn write_pfm(path: &str, rgb: &[Float], res_x: u32, res_y: u32) -> Result<(), String> {
    info!("Writing image {} with resolution {}x{}", path, res_x, res_y);

    let width = res_x as usize;
    let height = res_y as usize;

    let mut file = File::create(path).map_err(|e| format!("write_pfm(): Error writing PFM file '{path}': {e}"))?;

    // Only write 3 channel PFMs here.
    write_pfm_string(&mut file, "PF\n").map_err(|e| format!("write_pfm(): Error writing PFM type '{path}': {e}"))?;

    // Write the width and height, which must be positive.
    write_pfm_string(&mut file, format!("{width} {height}\n").as_str())
        .map_err(|e| format!("write_pfm(): Error writing PFM dimensions '{path}': {e}"))?;

    // Write the scale, which encodes endianness.
    let scale = if *is_big_endian() { 1.0 } else { -1.0 };
    write_pfm_string(&mut file, format!("{scale}\n").as_str())
        .map_err(|e| format!("write_pfm(): Error writing PFM scale '{path}': {e}"))?;

    // Write the data from bottom left to upper right as specified by http://netpbm.sourceforge.net/doc/pfm.html
    // The raster is a sequence of pixels, packed one after another, with no delimiters of any kind. They are grouped by
    // row, with the pixels in each row ordered left to right and the rows ordered bottom to top.
    for y in (0..=height - 1).rev() {
        // In case Float is 'double', copy into a staging buffer that's definitely a 32-bit float.
        for x in 0..3 * width {
            let f = rgb[y * width * 3 + x];
            if *is_big_endian() {
                file.write_f32::<BigEndian>(f)
                    .map_err(|e| format!("write_pfm(): Error writing PFM pixels '{path}': {e}"))?;
            } else {
                file.write_f32::<LittleEndian>(f)
                    .map_err(|e| format!("write_pfm(): Error writing PFM pixels '{path}': {e}"))?;
            }
        }
    }

    Ok(())
}

/// Apply gamma correction to a RGB floating point pixel and return the clamped 8-bit values.
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
