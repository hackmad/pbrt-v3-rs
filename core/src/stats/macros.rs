//! Macros

/// Create a thread local variable to track an `i64` counter across threads.
///
/// * `$title`     - Descriptive title of the statistic that uses `/` as a separator for categories.
///                  For example: "Intersections/Regular ray intersection tests",
/// * `$var`       - An identifier for the thread local variable.
/// * `stats_func` - An identifier for the callback function used by `StatsRegistrar::call_stats_funcs() to report to
///                 `StatsAccumulator`.
#[macro_export]
macro_rules! stat_counter {
    ($title: expr, $var: ident, $stats_func: ident $(,)?) => {
        thread_local! { static $var: std::cell::RefCell<i64> = std::cell::RefCell::new(0); }

        pub fn $stats_func(accum: &mut StatsAccumulator) {
            // Report thread stats.
            let val = $var.with(|v| *v.borrow());
            accum.report_counter($title, val);

            // Reset thread stats.
            $var.with(|v| {
                *v.borrow_mut() = 0;
            });
        }
    };
}

/// Create thread local variables to track an `i64` values for numerator/denominator across threads.
///
/// * `$title`     - Descriptive title of the statistic that uses `/` as a separator for categories.
///                  For example: "Intersections/Regular ray intersection tests",
/// * `$var_num`   - An identifier for the thread local variable for numerator (actual count).
/// * `$var_denom` - An identifier for the thread local variable for denominator (total count).
/// * `stats_func` - An identifier for the callback function used by `StatsRegistrar::call_stats_funcs() to report to
///                 `StatsAccumulator`.
#[macro_export]
macro_rules! stat_percent {
    ($title: expr, $var_num: ident, $var_denom: ident, $stats_func: ident $(,)?) => {
        thread_local! {
            static $var_num: std::cell::RefCell<i64> = std::cell::RefCell::new(0);
            static $var_denom: std::cell::RefCell<i64> = std::cell::RefCell::new(0);
        }

        pub fn $stats_func(accum: &mut StatsAccumulator) {
            // Report thread stats.
            let num = $var_num.with(|v| *v.borrow());
            let denom = $var_denom.with(|v| *v.borrow());
            accum.report_percentage($title, num, denom);

            // Reset thread stats.
            $var_num.with(|v| {
                *v.borrow_mut() = 0;
            });
            $var_denom.with(|v| {
                *v.borrow_mut() = 0;
            });
        }
    };
}

/// Convenience macro to increment a thread local variable for counter/percent statistics.
#[macro_export]
macro_rules! stat_inc {
    ($var: ident) => {
        $var.with(|v| *v.borrow_mut() += 1);
    };
}

/// Convenience macro to decrement a thread local variable for counter/percent statistics.
#[macro_export]
macro_rules! stat_dec {
    ($var: ident) => {
        $var.with(|v| *v.borrow_mut() -= 1);
    };
}

/// Convenience macro to register the callback functions for statistics.
///
/// * `$($func: ident),+` - One or more callback functions created by the `stat_*` macros.
#[macro_export]
macro_rules! register_stats {
    ($($stat_func: ident),+ $(,)?) => {
        lazy_static! {
            /// Used to ensure stats are registered exactly once in the module's private scope.
            static ref IS_STATS_REGISTERED: std::sync::Mutex<bool> = std::sync::Mutex::new(false);
        }

        /// Call this function in a module core/top-level struct to register the statistics. Typically done in `new()`.
        /// Avoid unnecessarily calling it from all structs; for example in `TriangleMesh::new()` vs `Triangle::new()`.
        fn register_stats() {
            let mut is_registered = IS_STATS_REGISTERED.lock().unwrap();
            if !*is_registered  {
                let mut sr = STATS_REGISTRAR.lock().unwrap();
                $(
                    sr.register_stat_func($stat_func);
                )+
                *is_registered = true;
            }
        }
    };
}

/// Convenience function to accumulate thread local statistics in `STATS_ACCUMULATOR`. This will call the registered
/// callbacks created with `stat_*` macros. This should be called at the end of each spawned thread and at the end of
/// rendering a scene from the main thread.
#[macro_export]
macro_rules! report_stats {
    () => {{
        let mut accum = STATS_ACCUMULATOR.lock().unwrap();
        STATS_REGISTRAR.lock().unwrap().call_stat_funcs(&mut accum);
    }};
}

/// Convenience function to print accumulated statistic in `STATS_ACCUMULATOR`.
#[macro_export]
macro_rules! print_stats {
    () => {{
        STATS_ACCUMULATOR.lock().unwrap().print();
    }};
}

/// Convenience function to clear accumulated statistic in `STATS_ACCUMULATOR`.
#[macro_export]
macro_rules! clear_stats {
    () => {{
        STATS_ACCUMULATOR.lock().unwrap().clear();
    }};
}
