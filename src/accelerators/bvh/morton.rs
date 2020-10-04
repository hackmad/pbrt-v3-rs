//! Morton Codes

#![allow(dead_code)]
use super::{float_to_bits, Vector3f};

/// Stores Morton codes (interleaved bits of coordinate values).
#[derive(Copy, Clone, Default, Debug)]
pub struct MortonPrimitive {
    /// The index of the primitive in a `primitive_info` vector.
    pub primitive_index: usize,

    /// The Morton code.
    pub morton_code: u32,
}

/// Create a new `MortonPrimitive`.
///
/// * `primitive_index` - The index of the primitive in a `primitive_info` vector.
/// * `morton_code`     - The Morton code.
pub fn morton_primitive(primitive_index: usize, morton_code: u32) -> MortonPrimitive {
    MortonPrimitive {
        primitive_index,
        morton_code,
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
    (left_shift_3(float_to_bits(v.z)) << 2)
        | (left_shift_3(float_to_bits(v.y)) << 1)
        | left_shift_3(float_to_bits(v.x))
}

pub const BITS_PER_PASS: usize = 6;
pub const N_BITS: usize = 30;
pub const N_PASSES: usize = N_BITS / BITS_PER_PASS; // This should divide evenly.
pub const N_BUCKETS: usize = 12;
pub const FIRST_BIT_INDEX: usize = N_BITS - 1 - N_BUCKETS; // The index of the next
                                                           // bit to try splitting

/// Sorts a list of morton primitives in place using Radix sort.
///
/// *NOTE*: This version is copying the entire vector and uses a second vector
/// of same size swapping between their references to finally copy the result.
/// So it'll be a bit more memory intensive.
///
/// * `v` - The morton primitives.
pub fn radix_sort(v: &Vec<MortonPrimitive>) -> Vec<MortonPrimitive> {
    const N_BUCKETS: usize = 1 << BITS_PER_PASS;
    const BIT_MASK: usize = (1 << BITS_PER_PASS) - 1;

    // Copy the input vector and initialize a second one.
    let mut v1 = v.clone();
    let mut v2 = vec![MortonPrimitive::default(); v.len()];

    for pass in 0..N_PASSES {
        // Perform one pass of radix sort, sorting BITS_PER_PASS bits.
        let low_bit = pass * BITS_PER_PASS;

        // Set in and out vector pointers for next radix sort pass.
        let (v_in, v_out) = if pass & 1 == 1 {
            (&mut v2, &mut v1)
        } else {
            (&mut v1, &mut v2)
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
        v2
    } else {
        v1
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
    let mut x1 = if x == (1 << 10) { x - 1 } else { x };

    x1 = (x1 | (x1 << 16)) & 0b00000011000000000000000011111111;
    x1 = (x1 | (x1 << 8))  & 0b00000011000000001111000000001111;
    x1 = (x1 | (x1 << 4))  & 0b00000011000011000011000011000011;
    x1 = (x1 | (x1 << 2))  & 0b00001001001001001001001001001001;

    x1
}
