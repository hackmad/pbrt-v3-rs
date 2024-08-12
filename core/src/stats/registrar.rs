//! Statistics Regisration

use super::StatsAccumulator;
use std::sync::LazyLock;
use std::sync::Mutex;

/// The global registrar for statistics.
pub static STATS_REGISTRAR: LazyLock<Mutex<StatsRegistrar>> = LazyLock::new(|| Mutex::new(StatsRegistrar::new()));

/// Registers callback functions for statistics.
pub struct StatsRegistrar {
    /// Callback functions.
    stats_funcs: Vec<fn(&mut StatsAccumulator) -> ()>,
}

impl StatsRegistrar {
    /// Create a new instance of `StatsRegistrar`.
    pub fn new() -> Self {
        Self { stats_funcs: vec![] }
    }

    /// Register a callback function for reporting statistics.
    ///
    /// * `func` - A callback function that takes a `StatsAccumulator` to report statistics.
    pub fn register_stat_func(&mut self, func: fn(&mut StatsAccumulator) -> ()) {
        self.stats_funcs.push(func);
    }

    /// Call all callback functions for reporting statistics.
    ///
    /// * `func` - A callback function that takes a `StatsAccumulator` to report statistics.
    pub fn call_stat_funcs(&self, accum: &mut StatsAccumulator) {
        self.stats_funcs.iter().for_each(|func| {
            func(accum);
        });
    }
}
