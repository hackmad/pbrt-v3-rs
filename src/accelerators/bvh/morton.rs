//! Morton Codes

#![allow(dead_code)]
use super::{float_to_bits, Vector3f};
use std::cell::RefCell;

/// Stores Morton codes (interleaved bits of coordinate values).
#[derive(Copy, Clone, Default, Debug)]
pub struct MortonPrimitive {
    /// The index of the primitive in a `primitive_info` vector.
    pub primitive_index: usize,

    /// The Morton code.
    pub morton_code: u32,
}

impl MortonPrimitive {
    /// Create a new `MortonPrimitive`.
    ///
    /// * `primitive_index` - The index of the primitive in a `primitive_info` vector.
    /// * `morton_code`     - The Morton code.
    pub fn new(primitive_index: usize, morton_code: u32) -> Self {
        Self {
            primitive_index,
            morton_code,
        }
    }
}

/// Take a 3D coordinate value where each component is a floating-point value
/// between 0 and 2^10 and convert these values to integers and then computes the
/// Morton code by expanding the three 10-bit quantized values so that their i^th
/// bits are at position `3i`, then shifting the bits over one more, the z bits
/// over two more, and ORing together the result.
///
/// * `v` - The coordinate value.
pub fn encode_morton_3(v: &Vector3f) -> u32 {
    debug_assert!(v.x > 0.0);
    debug_assert!(v.y > 0.0);
    debug_assert!(v.z > 0.0);

    (left_shift_3(float_to_bits(v.z)) << 2)
        | (left_shift_3(float_to_bits(v.y)) << 1)
        | left_shift_3(float_to_bits(v.x))
}

pub const N_BITS: usize = 30;
const BITS_PER_PASS: usize = 6;
const N_PASSES: usize = N_BITS / BITS_PER_PASS; // This should divide evenly.
const N_BUCKETS: usize = 1 << BITS_PER_PASS;
const BIT_MASK: usize = (1 << BITS_PER_PASS) - 1;

/// Sorts a list of morton primitives in place using Radix sort.
///
/// * `v` - A `RefCell` containing the morton primitives vector.
pub fn radix_sort(v: &mut RefCell<Vec<MortonPrimitive>>) {
    let n = { v.borrow().len() };
    let mut temp_vector = RefCell::new(vec![MortonPrimitive::default(); n]);

    for pass in 0..N_PASSES {
        // Perform one pass of radix sort, sorting BITS_PER_PASS bits.
        let low_bit = pass * BITS_PER_PASS;

        // Set in and out vector pointers for radix sort pass.
        let (v_in, v_out) = if pass & 1 == 1 {
            (temp_vector.borrow(), v.get_mut())
        } else {
            (v.borrow(), temp_vector.get_mut())
        };

        // Count number of zero bits in array for current radix sort bit.
        let mut bucket_count = [0; N_BUCKETS];
        for mp in v_in.iter() {
            let bucket = ((mp.morton_code >> low_bit) as usize) & BIT_MASK;
            debug_assert!(bucket < N_BUCKETS);
            bucket_count[bucket] += 1;
        }

        // Compute starting index in output array for each bucket.
        let mut out_index = [0; N_BUCKETS];
        for i in 1..N_BUCKETS {
            // want out_index[0] = 0.
            out_index[i] = out_index[i - 1] + bucket_count[i - 1];
        }

        // Store sorted values in output array.
        for mp in v_in.iter() {
            let bucket = ((mp.morton_code >> low_bit) as usize) & BIT_MASK;
            v_out[out_index[bucket]] = *mp;
            out_index[bucket] += 1;
        }
    }

    // Copy final result from temp_vector, if needed.
    if N_PASSES & 1 == 1 {
        *v.get_mut() = temp_vector.into_inner();
    }
}

/// The bit shifts to compute the Morton code for each 3D coordinate are 
/// performed in a series of shifts of power-of-two size. First, bits 8 and 9
/// are shifted 16 places to the left. This places bit 8 in its final position. 
/// Next bits 4 through 7 are shifted 8 places. After shifts of 4 and 2 places 
/// (with appropriate masking so that each bit is shifted the right number of 
/// places in the end), all bits are in the proper position. 
///
/// * `x` - The value.
#[rustfmt::skip]
fn left_shift_3(x: u32) -> u32 {
    debug_assert!(x < 1 << 10);

    let mut x1 = if x == (1 << 10) { x - 1 } else { x };

    x1 = (x1 | (x1 << 16)) & 0b00000011000000000000000011111111;
    // x1 = ---- --98 ---- ---- ---- ---- 7654 3210
    x1 = (x1 | (x1 << 8))  & 0b00000011000000001111000000001111;
    // x1 = ---- --98 ---- ---- 7654 ---- ---- 3210
    x1 = (x1 | (x1 << 4))  & 0b00000011000011000011000011000011;
    // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x1 = (x1 | (x1 << 2))  & 0b00001001001001001001001001001001;
    // x1 = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0

    x1
}
