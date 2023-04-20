//! Statistics Accumulator

use std::collections::HashMap;
use std::sync::Mutex;

lazy_static! {
    pub static ref STATS_ACCUMULATOR: Mutex<StatsAccumulator> = Mutex::new(StatsAccumulator::new());
}

pub struct StatsAccumulator {
    counters: HashMap<String, i64>,
    percentages: HashMap<String, (i64, i64)>,
}

impl StatsAccumulator {
    pub fn new() -> Self {
        Self {
            counters: HashMap::new(),
            percentages: HashMap::new(),
        }
    }

    pub fn report_counter(&mut self, name: &str, val: i64) {
        if let Some(v) = self.counters.get_mut(name) {
            *v += val;
        } else {
            self.counters.insert(name.to_string(), val);
        }
    }

    pub fn report_percentage(&mut self, name: &str, num: i64, denom: i64) {
        if let Some(v) = self.percentages.get_mut(name) {
            v.0 += num;
            v.1 += denom;
        } else {
            self.percentages.insert(name.to_string(), (num, denom));
        }
    }

    pub fn print(&self) {
        let mut to_print: HashMap<String, Vec<String>> = HashMap::new();

        for (k, v) in self.counters.iter() {
            if *v == 0 {
                continue;
            }

            let (category, title) = get_category_and_title(k);
            let s = format!("{:-42}               {:12}", title, v).to_string();

            if let Some(list) = to_print.get_mut(&category) {
                list.push(s);
            } else {
                to_print.insert(category, vec![s]);
            }
        }

        for (k, (num, denom)) in self.percentages.iter() {
            if *denom == 0 {
                continue;
            }
            let (category, title) = get_category_and_title(k);
            let s = format!(
                "{:-42}{:12} / {:12} ({:.2}%)",
                title,
                *num,
                *denom,
                (100.0 * *num as f64) / *denom as f64
            );
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

    pub fn clear(&mut self) {
        self.counters.clear();
    }
}

fn get_category_and_title(s: &str) -> (String, String) {
    if let Some(slash) = s.find("/") {
        let category = &s[0..slash];
        let title = &s[slash + 1..];
        (category.to_string(), title.to_string())
    } else {
        ("".to_string(), s.to_string())
    }
}
