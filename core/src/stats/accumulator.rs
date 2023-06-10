//! Statistics Accumulator

use crate::pbrt;
use num_traits::{Num, Zero};
use std::collections::HashMap;
use std::ops::AddAssign;
use std::sync::Mutex;
use std::sync::OnceLock;

/// Return the global statistics accumulator.
pub fn stats_accumulator() -> &'static Mutex<StatsAccumulator> {
    static DATA: OnceLock<Mutex<StatsAccumulator>> = OnceLock::new();
    DATA.get_or_init(|| Mutex::new(StatsAccumulator::new()))
}

/// Distribution statistic.
#[derive(Default, Clone)]
pub struct StatsDistribution<T>
where
    T: Num + Default + Copy + Clone,
{
    /// Sum of all values.
    sum: T,

    /// Count of all values.
    count: u64,

    /// Minimum value.
    min: Option<T>,

    /// Maximum value.
    max: Option<T>,
}

impl<T> StatsDistribution<T>
where
    T: Num + Zero + PartialOrd + AddAssign + Default + Copy + Clone,
{
    /// Create a new instance of `StatsDistribution<T>`.
    ///
    /// * `sum`   - Sum of values.
    /// * `count` - Count of values.
    /// * `min`   - Minimum value.
    /// * `max`   - Maximum value.
    pub fn new(sum: T, count: u64, min: T, max: T) -> Self {
        Self {
            sum,
            count,
            min: Some(min),
            max: Some(max),
        }
    }

    /// Accumulate stats.
    ///
    /// * `sum`   - Sum of values.
    /// * `count` - Count of values.
    /// * `min`   - Minimum value.
    /// * `max`   - Maximum value.
    pub fn accumulate(&mut self, distrib: Self) {
        self.sum += distrib.sum;
        self.count += distrib.count;

        if let Some(v) = self.min.as_mut() {
            if let Some(min) = distrib.min {
                *v = pbrt::min(*v, min);
            }
        } else {
            self.min = distrib.min;
        }

        if let Some(v) = self.max.as_mut() {
            if let Some(max) = distrib.max {
                *v = pbrt::max(*v, max);
            }
        } else {
            self.max = distrib.max;
        }
    }

    /// Report a sample value.
    ///
    /// * `val`  - Sample value.
    pub fn report(&mut self, val: T) {
        self.sum += val;
        self.count += 1;

        if let Some(v) = self.min.as_mut() {
            *v = pbrt::min(*v, val);
        } else {
            self.min = Some(val);
        }

        if let Some(v) = self.max.as_mut() {
            *v = pbrt::max(*v, val);
        } else {
            self.max = Some(val);
        }
    }

    /// Clear stats.
    pub fn clear(&mut self) {
        self.sum = T::zero();
        self.count = 0;
        self.min = None;
        self.max = None;
    }
}

/// Aggregate different types of statistics.
pub struct StatsAccumulator {
    /// Counters.
    counters: HashMap<String, i64>,

    /// Memory counters.
    memory_counters: HashMap<String, u64>,

    /// Integer distribution.
    int_distribution: HashMap<String, StatsDistribution<i64>>,

    /// Float distribution.
    float_distribution: HashMap<String, StatsDistribution<f64>>,

    /// Percentages.
    percentages: HashMap<String, (i64, i64)>,

    /// Ratios.
    ratios: HashMap<String, (i64, i64)>,
}

impl StatsAccumulator {
    /// Create a new instance of `StatsAccumulator`.
    pub fn new() -> Self {
        Self {
            counters: HashMap::new(),
            memory_counters: HashMap::new(),
            int_distribution: HashMap::new(),
            float_distribution: HashMap::new(),
            percentages: HashMap::new(),
            ratios: HashMap::new(),
        }
    }

    /// Accumulates a counter value.
    ///
    /// * `name` - Statistic name.
    /// * `val`  - Counter value.
    pub fn report_counter(&mut self, name: &str, val: i64) {
        if let Some(v) = self.counters.get_mut(name) {
            *v += val;
        } else {
            self.counters.insert(name.to_string(), val);
        }
    }

    /// Accumulates a memory counter value.
    ///
    /// * `name` - Statistic name.
    /// * `val`  - Memory counter value.
    pub fn report_memory_counter(&mut self, name: &str, val: u64) {
        if let Some(v) = self.memory_counters.get_mut(name) {
            *v += val;
        } else {
            self.memory_counters.insert(name.to_string(), val);
        }
    }

    /// Accumulates integer point distribution samples.
    ///
    /// * `name`    - Statistic name.
    /// * `distrib` - Distribution.
    pub fn report_int_distribution(&mut self, name: &str, distrib: StatsDistribution<i64>) {
        if let Some(v) = self.int_distribution.get_mut(name) {
            v.accumulate(distrib);
        } else {
            self.int_distribution.insert(name.to_string(), distrib);
        }
    }

    /// Accumulates floating point distribution samples.
    ///
    /// * `name`  - Statistic name.
    /// * `distrib` - Distribution.
    pub fn report_float_distribution(&mut self, name: &str, distrib: StatsDistribution<f64>) {
        if let Some(v) = self.float_distribution.get_mut(name) {
            v.accumulate(distrib);
        } else {
            self.float_distribution.insert(name.to_string(), distrib);
        }
    }

    /// Accumulates a percentage value.
    ///
    /// * `name`  - Statistic name.
    /// * `num`   - Numerator (actual count).
    /// * `denom` - Denominator (total count).
    pub fn report_percentage(&mut self, name: &str, num: i64, denom: i64) {
        if let Some(v) = self.percentages.get_mut(name) {
            v.0 += num;
            v.1 += denom;
        } else {
            self.percentages.insert(name.to_string(), (num, denom));
        }
    }

    /// Accumulates a ratio value.
    ///
    /// * `name`  - Statistic name.
    /// * `num`   - Numerator.
    /// * `denom` - Denominator.
    pub fn report_ratio(&mut self, name: &str, num: i64, denom: i64) {
        if let Some(v) = self.ratios.get_mut(name) {
            v.0 += num;
            v.1 += denom;
        } else {
            self.ratios.insert(name.to_string(), (num, denom));
        }
    }

    /// Prints the report.
    pub fn print(&self) {
        let mut to_print: HashMap<String, Vec<String>> = HashMap::new();

        for (k, v) in self.counters.iter() {
            if *v == 0 {
                continue;
            }

            let (category, title) = get_category_and_title(k);
            let s = format!("{title:-42}               {v:12}");

            if let Some(list) = to_print.get_mut(&category) {
                list.push(s);
            } else {
                to_print.insert(category, vec![s]);
            }
        }

        for (k, v) in self.memory_counters.iter() {
            if *v == 0 {
                continue;
            }
            let (category, title) = get_category_and_title(k);
            let kb = *v as f64 / 1024.0;
            let s = if kb < 1024.0 {
                format!("{title:-42}                  {kb:9.2} kB")
            } else {
                let mib = kb / 1024.0;
                if mib < 1024.0 {
                    format!("{title:-42}                  {mib:9.2} MiB")
                } else {
                    let gib = mib / 1024.0;
                    format!("{title:-42}                  {gib:9.2} GiB")
                }
            };
            if let Some(list) = to_print.get_mut(&category) {
                list.push(s);
            } else {
                to_print.insert(category, vec![s]);
            }
        }

        for (k, v) in self.int_distribution.iter() {
            if v.count == 0 {
                continue;
            }

            let mn = if let Some(m) = v.min { m } else { i64::MAX };
            let mx = if let Some(m) = v.max { m } else { i64::MIN };

            let (category, title) = get_category_and_title(&k);
            let avg = v.sum as f64 / v.count as f64;
            let s = format!("{title:-42}                      {avg:.3} avg [range {mn} - {mx}]");
            if let Some(list) = to_print.get_mut(&category) {
                list.push(s);
            } else {
                to_print.insert(category, vec![s]);
            }
        }

        for (k, v) in self.float_distribution.iter() {
            if v.count == 0 {
                continue;
            }

            let mn = if let Some(m) = v.min { m } else { f64::MAX };
            let mx = if let Some(m) = v.max { m } else { f64::MIN };

            let (category, title) = get_category_and_title(&k);
            let avg = v.sum / v.count as f64;
            let s = format!("{title:-42}                      {avg:.3} avg [range {mn} - {mx}]");
            if let Some(list) = to_print.get_mut(&category) {
                list.push(s);
            } else {
                to_print.insert(category, vec![s]);
            }
        }

        for (k, &(num, denom)) in self.percentages.iter() {
            if denom == 0 {
                continue;
            }
            let (category, title) = get_category_and_title(k);
            let s = format!(
                "{title:-42}{num:12} / {denom:12} ({:.2}%)",
                (100.0 * num as f64) / denom as f64,
            );
            if let Some(list) = to_print.get_mut(&category) {
                list.push(s);
            } else {
                to_print.insert(category, vec![s]);
            }
        }

        for (k, &(num, denom)) in self.ratios.iter() {
            if denom == 0 {
                continue;
            }
            let (category, title) = get_category_and_title(k);
            let s = format!("{title:-42}{num:12} / {denom:12} ({:.2}x)", num as f64 / denom as f64);
            if let Some(list) = to_print.get_mut(&category) {
                list.push(s);
            } else {
                to_print.insert(category, vec![s]);
            }
        }

        println!("Statistics:");
        for (category, items) in to_print {
            println!("  {category}");
            for item in items {
                println!("    {item}");
            }
        }
    }

    /// Clear the accumulated statistics.
    pub fn clear(&mut self) {
        self.counters.clear();
        self.memory_counters.clear();
        self.int_distribution.clear();
        self.float_distribution.clear();
        self.percentages.clear();
        self.ratios.clear();
    }
}

/// Splits a statistic name at the first `/` as the separator and returns category and title. If there is no `/`, then
/// category is the empty string.
///
/// * `s` - The statistic name to split.
fn get_category_and_title(s: &str) -> (String, String) {
    if let Some(slash) = s.find("/") {
        let category = &s[0..slash];
        let title = &s[slash + 1..];
        (category.to_string(), title.to_string())
    } else {
        ("".to_string(), s.to_string())
    }
}
