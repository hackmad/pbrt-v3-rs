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
        thread_local! { pub(crate) static $var: std::cell::RefCell<i64> = std::cell::RefCell::new(0); }

        pub(crate) fn $stats_func(accum: &mut StatsAccumulator) {
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

/// Create a thread local variable to track memory usage as a counter across threads.
///
/// * `$title`     - Descriptive title of the statistic that uses `/` as a separator for categories.
///                  For example: "Intersections/Regular ray intersection tests",
/// * `$var`       - An identifier for the thread local variable.
/// * `stats_func` - An identifier for the callback function used by `StatsRegistrar::call_stats_funcs() to report to
///                 `StatsAccumulator`.
#[macro_export]
macro_rules! stat_memory_counter {
    ($title: expr, $var: ident, $stats_func: ident $(,)?) => {
        thread_local! { pub(crate) static $var: std::cell::RefCell<u64> = std::cell::RefCell::new(0); }

        pub(crate) fn $stats_func(accum: &mut StatsAccumulator) {
            // Report thread stats.
            let val = $var.with(|v| *v.borrow());
            accum.report_memory_counter($title, val);

            // Reset thread stats.
            $var.with(|v| {
                *v.borrow_mut() = 0;
            });
        }
    };
}

/// Create a thread local variable to track an integer distribution across threads.
///
/// * `$title`     - Descriptive title of the statistic that uses `/` as a separator for categories.
///                  For example: "Intersections/Regular ray intersection tests",
/// * `$var`       - An identifier for the thread local variable.
/// * `stats_func` - An identifier for the callback function used by `StatsRegistrar::call_stats_funcs() to report to
///                 `StatsAccumulator`.
#[macro_export]
macro_rules! stat_int_distribution {
    ($title: expr, $var: ident, $stats_func: ident $(,)?) => {
        thread_local! {
            pub(crate) static $var: std::cell::RefCell<StatsDistribution<i64>> =
                std::cell::RefCell::new(StatsDistribution::default());
        }

        pub(crate) fn $stats_func(accum: &mut StatsAccumulator) {
            // Report thread stats.
            let val = $var.with(|v| v.borrow().clone());
            accum.report_int_distribution($title, val);

            // Reset thread stats.
            $var.with(|v| {
                v.borrow_mut().clear();
            });
        }
    };
}

/// Create a thread local variable to track an float distribution across threads.
///
/// * `$title`     - Descriptive title of the statistic that uses `/` as a separator for categories.
///                  For example: "Intersections/Regular ray intersection tests",
/// * `$var`       - An identifier for the thread local variable.
/// * `stats_func` - An identifier for the callback function used by `StatsRegistrar::call_stats_funcs() to report to
///                 `StatsAccumulator`.
#[macro_export]
macro_rules! stat_float_distribution {
    ($title: expr, $var: ident, $stats_func: ident $(,)?) => {
        thread_local! {
            pub(crate) static $var: std::cell::RefCell<StatsDistribution<f64>> =
                std::cell::RefCell::new(StatsDistribution::default());
        }

        pub(crate) fn $stats_func(accum: &mut StatsAccumulator) {
            // Report thread stats.
            let val = $var.with(|v| v.borrow().clone());
            accum.report_float_distribution($title, val);

            // Reset thread stats.
            $var.with(|v| {
                v.borrow_mut().clear();
            });
        }
    };
}

/// Create thread local variables to track an `i64` values for numerator/denominator as percentage across threads.
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
            pub(crate) static $var_num: std::cell::RefCell<i64> = std::cell::RefCell::new(0);
            pub(crate) static $var_denom: std::cell::RefCell<i64> = std::cell::RefCell::new(0);
        }

        pub(crate) fn $stats_func(accum: &mut StatsAccumulator) {
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

/// Create thread local variables to track an `i64` values for numerator/denominator as ratio across threads.
///
/// * `$title`     - Descriptive title of the statistic that uses `/` as a separator for categories.
///                  For example: "Intersections/Regular ray intersection tests",
/// * `$var_num`   - An identifier for the thread local variable for numerator (actual count).
/// * `$var_denom` - An identifier for the thread local variable for denominator (total count).
/// * `stats_func` - An identifier for the callback function used by `StatsRegistrar::call_stats_funcs() to report to
///                 `StatsAccumulator`.
#[macro_export]
macro_rules! stat_ratio {
    ($title: expr, $var_num: ident, $var_denom: ident, $stats_func: ident $(,)?) => {
        thread_local! {
            pub(crate) static $var_num: std::cell::RefCell<i64> = std::cell::RefCell::new(0);
            pub(crate) static $var_denom: std::cell::RefCell<i64> = std::cell::RefCell::new(0);
        }

        pub(crate) fn $stats_func(accum: &mut StatsAccumulator) {
            // Report thread stats.
            let num = $var_num.with(|v| *v.borrow());
            let denom = $var_denom.with(|v| *v.borrow());
            accum.report_ratio($title, num, denom);

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
    ($var: ident, $e: expr) => {
        $var.with(|v| *v.borrow_mut() += $e);
    };
}

/// Convenience macro to decrement a thread local variable for counter/percent statistics.
#[macro_export]
macro_rules! stat_dec {
    ($var: ident, $e: expr) => {
        $var.with(|v| *v.borrow_mut() -= $e);
    };
}

/// Convenience macro to report a thread local variable for int/float distribution statistics.
#[macro_export]
macro_rules! stat_dist {
    ($var: ident, $e: expr) => {
        $var.with(|v| v.borrow_mut().report($e));
    };
}

/// Convenience macro to register the callback functions for statistics.
///
/// * `$($func: ident),+` - One or more callback functions created by the `stat_*` macros.
#[macro_export]
macro_rules! stat_register_fns {
    ($($stat_func: ident),+ $(,)?) => {
        /// Return the whether stats are registered.
        pub(crate) fn is_stats_registered() -> &'static std::sync::Mutex<bool> {
            static DATA: std::sync::OnceLock<std::sync::Mutex<bool>> = std::sync::OnceLock::new();
            DATA.get_or_init(|| std::sync::Mutex::new(false))
        }

        /// Call this function in a module core/top-level struct to register the statistics. Typically done in:
        /// 1. Type constructor like new(). Avoid unnecessarily calling it from all structs.
        ///    e.g. `TriangleMesh::new()` (but not `Triangle::new()` as this would result in a huge number of calls)
        /// 2. impl From<(&ParamSet, ...)> for ... { ... }.
        ///    e.g. GonioPhotometricLight
        /// 3. pub fn from_props(p: (&ParamSet, ...)) { ... }
        ///    e.g. PLYMesh
        pub(crate) fn register_stats() {
            let mut is_registered = is_stats_registered().lock().unwrap();
            if !*is_registered  {
                let mut sr = stats_registrar().lock().unwrap();
                $(
                    sr.register_stat_func($stat_func);
                )+
                *is_registered = true;
            }
        }
    };
}

/// Convenience function to accumulate thread local statistics in the global `StatsAccumulator`. This will call the
/// registered callbacks created with `stat_*` macros. This should be called at the end of each spawned thread and at
/// the end of rendering a scene from the main thread.
#[macro_export]
macro_rules! report_stats {
    () => {{
        let mut accum = stats_accumulator().lock().unwrap();
        stats_registrar().lock().unwrap().call_stat_funcs(&mut accum);
    }};
}

/// Convenience function to print accumulated statistic in the global `StatsAccumulator`.
#[macro_export]
macro_rules! print_stats {
    () => {{
        stats_accumulator().lock().unwrap().print();
    }};
}

/// Convenience function to clear accumulated statistic in the global `StatsAccumulator`.
#[macro_export]
macro_rules! clear_stats {
    () => {{
        stats_accumulator().lock().unwrap().clear();
    }};
}
