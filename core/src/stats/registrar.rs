//! Statistics Regisration

use super::StatsAccumulator;
use std::sync::Mutex;

lazy_static! {
    pub static ref STATS_REGISTRAR: Mutex<StatsRegistrar> = Mutex::new(StatsRegistrar::new());
}

pub struct StatsRegistrar {
    stats_funcs: Vec<fn(&mut StatsAccumulator) -> ()>,
}

impl StatsRegistrar {
    pub fn new() -> Self {
        Self { stats_funcs: vec![] }
    }

    pub fn register_stat_func(&mut self, func: fn(&mut StatsAccumulator) -> ()) {
        self.stats_funcs.push(func);
    }

    pub fn call_stat_funcs(&self, accum: &mut StatsAccumulator) {
        self.stats_funcs.iter().for_each(|func| {
            func(accum);
        });
    }
}
