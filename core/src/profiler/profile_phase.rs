//! Profile Phase

use super::Prof;
use super::PROFILER_STATE;

pub struct ProfilePhase {
    pub(super) reset: bool,
    pub(super) category_bit: u64,
}

impl ProfilePhase {
    pub fn new(p: Prof) -> Self {
        let category_bit = p.to_bits();

        PROFILER_STATE.with(|s| {
            let reset = *s.borrow() & category_bit == 0;
            *s.borrow_mut() |= category_bit;

            Self { category_bit, reset }
        })
    }
}

impl Drop for ProfilePhase {
    fn drop(&mut self) {
        if self.reset {
            PROFILER_STATE.with(|s| {
                *s.borrow_mut() &= !self.category_bit;
            });
        }
    }
}
