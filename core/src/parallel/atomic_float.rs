//! AtomicFloat

use crate::pbrt::{bits_to_float, float_to_bits, Float};
use std::sync::atomic::{AtomicU32, Ordering};

/// Implement atomic floating point value using `AtomicU32`.
pub struct AtomicFloat {
    /// Bit representation of floating point value.
    bits: AtomicU32, // Use AtomicU64 when Float = f64
}

impl AtomicFloat {
    /// Create a new `AtomicFloat`.
    ///
    /// * `v` - The value.
    pub fn new(v: Float) -> Self {
        Self {
            bits: AtomicU32::new(float_to_bits(v)),
        }
    }

    /// Add a floating point value.
    ///
    /// * `v` - The value to add.
    pub fn add(&self, v: Float) {
        let mut old_bits: u32 = self.bits.load(Ordering::Relaxed);
        loop {
            let new_bits = float_to_bits(bits_to_float(old_bits) + v);
            let result = self.bits.compare_exchange_weak(
                old_bits,
                new_bits,
                Ordering::SeqCst,
                Ordering::Relaxed,
            );
            match result {
                Ok(_) => break,
                Err(x) => {
                    old_bits = x;
                }
            }
        }
    }

    /// Loads the floating point value.
    ///
    /// * `order` - Memory ordering of this operation
    pub fn load(&self, order: Ordering) -> Float {
        let bits: u32 = self.bits.load(order);
        bits_to_float(bits)
    }

    /// Stores the floating point value.
    ///
    /// * `v`     - The value.
    /// * `order` - Memory ordering of this operation
    pub fn store(&self, v: Float, order: Ordering) {
        self.bits.store(float_to_bits(v), order);
    }
}

impl Default for AtomicFloat {
    /// Returns the "default value" for `AtomicFloat`.
    fn default() -> Self {
        Self {
            bits: AtomicU32::new(0),
        }
    }
}

impl Clone for AtomicFloat {
    fn clone(&self) -> Self {
        let bits: u32 = self.bits.load(Ordering::SeqCst);
        AtomicFloat {
            bits: AtomicU32::new(bits),
        }
    }
}
