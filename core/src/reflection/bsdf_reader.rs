//! BSDF Reader

use crate::pbrt::is_big_endian;
use byteorder::{BigEndian, LittleEndian, ReadBytesExt};
use std::fs::File;
use std::io::Read;
use std::str;

/// The first 8 byetes of BSDF file are the header `SCATFUN` terminated with char `0x01`.
const EXPECTED_HEADER: [u8; 8] = [b'S', b'C', b'A', b'T', b'F', b'U', b'N', b'\x01'];

/// Opens a file for reading or returns an error if unable to do so.
///
/// * `path` - The file path.
pub fn open_file(path: &str) -> Result<File, String> {
    match File::open(path) {
        Ok(file) => Ok(file),
        Err(err) => Err(format!("Could not open {}. {}", path, err)),
    }
}

/// Interface to add custom helpers for reading BSDF files that supports little and big endian integer format of the
/// system.
///
/// `NOTE`: This is just a convenience way to add helpers to `File`.
pub trait BSDFReader {
    /// Reads the header bytes and compares them to the expected header.
    fn check_header(&mut self) -> Result<(), String>;

    /// Reads one 32-bit unsigned value.
    fn read_i32(&mut self) -> Result<i32, String>;

    /// Reads one 32-bit floating point value.
    fn read_f32(&mut self) -> Result<f32, String>;

    /// Reads given number of 32-bit unsigned values.
    ///
    /// * `count` - Number of values to read.
    fn read_i32_vec(&mut self, count: usize) -> Result<Vec<i32>, String>;

    /// Reads given number of 32-bit floating point values.
    ///
    /// * `count` - Number of values to read.
    fn read_f32_vec(&mut self, count: usize) -> Result<Vec<f32>, String>;
}

impl BSDFReader for File {
    /// Reads the header bytes and compares them to the expected header.
    fn check_header(&mut self) -> Result<(), String> {
        let mut header = [0_u8; 8];
        match self.read_exact(&mut header) {
            Ok(_) => {
                if header == EXPECTED_HEADER {
                    Ok(())
                } else {
                    Err(format!(
                        "Invalid header '{}'. Expected '{}'.",
                        str::from_utf8(&header).unwrap(),
                        str::from_utf8(&EXPECTED_HEADER).unwrap(),
                    ))
                }
            }
            Err(err) => Err(format!("Error reading header {:}", err)),
        }
    }

    /// Reads one 32-bit unsigned value.
    fn read_i32(&mut self) -> Result<i32, String> {
        match if *is_big_endian() {
            ReadBytesExt::read_i32::<BigEndian>(self)
        } else {
            ReadBytesExt::read_i32::<LittleEndian>(self)
        } {
            Ok(v) => Ok(v),
            Err(err) => Err(format!("Error reading one i32. {:}.", err)),
        }
    }

    /// Reads one 32-bit floating point value.
    fn read_f32(&mut self) -> Result<f32, String> {
        match if *is_big_endian() {
            ReadBytesExt::read_f32::<BigEndian>(self)
        } else {
            ReadBytesExt::read_f32::<LittleEndian>(self)
        } {
            Ok(v) => Ok(v),
            Err(err) => Err(format!("Error reading one i32. {:}.", err)),
        }
    }

    /// Reads given number of 32-bit unsigned values.
    ///
    /// * `count` - Number of values to read.
    fn read_i32_vec(&mut self, count: usize) -> Result<Vec<i32>, String> {
        let mut buffer = vec![0_i32; count];
        match if *is_big_endian() {
            ReadBytesExt::read_i32_into::<BigEndian>(self, &mut buffer)
        } else {
            ReadBytesExt::read_i32_into::<LittleEndian>(self, &mut buffer)
        } {
            Ok(_) => Ok(buffer),
            Err(err) => Err(format!("Error reading {} i32. {:}.", count, err)),
        }
    }

    /// Reads given number of 32-bit floating point values.
    ///
    /// * `count` - Number of values to read.
    fn read_f32_vec(&mut self, count: usize) -> Result<Vec<f32>, String> {
        let mut buffer: Vec<f32> = vec![0.0; count];
        match if *is_big_endian() {
            ReadBytesExt::read_f32_into::<BigEndian>(self, &mut buffer)
        } else {
            ReadBytesExt::read_f32_into::<LittleEndian>(self, &mut buffer)
        } {
            Ok(_) => Ok(buffer),
            Err(err) => Err(format!("Error reading {} f32. {:}.", count, err)),
        }
    }
}
