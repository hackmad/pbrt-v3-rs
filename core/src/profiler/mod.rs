//! Profiler

use crate::pbrt::{max, Log2Int};
use itertools::Itertools;
use std::cell::RefCell;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{channel, Sender, TryRecvError};
use std::sync::Mutex;
use std::thread;
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

mod prof;
mod profile_phase;
mod profile_sample;

use profile_sample::*;

pub use prof::*;
pub use profile_phase::*;

thread_local! {
    /// Current profiler state.
    pub static PROFILER_STATE: RefCell<u64> = RefCell::new(0);

    /// Used to send signal to profiler for cleanup.
    static SENDER: RefCell<Option<Sender<()>>> = RefCell::new(None);

    /// Handle for profiler thread.
    static HANDLE: RefCell<Option<JoinHandle<()>>> = RefCell::new(None);

    /// Time when profiler was initialized.
    static PROFILE_START_TIME: RefCell<Instant> = RefCell::new(Instant::now());
}

/// Indicates whether profiler is running or not.
static PROFILER_RUNNING: AtomicBool = AtomicBool::new(false);

/// Hash size for profiler samples.
const PROFILE_HASH_SIZE: usize = 256;

lazy_static! {
    /// Profiler samples.
    static ref PROFILE_SAMPLES: Mutex<Vec< ProfileSample >> =
        Mutex::new(vec![ProfileSample::default(); PROFILE_HASH_SIZE]);
}

/// Initialize the profiler thread.
pub fn init_profiler() {
    if PROFILER_RUNNING.load(Ordering::SeqCst) {
        panic!("Profiler is already running!");
    }

    PROFILE_START_TIME.with(|s| *s.borrow_mut() = Instant::now());

    PROFILER_STATE.with(|s| {
        *s.borrow_mut() = Prof::SceneConstruction.to_bits();
    });

    clear_profiler();

    let (sender, receiver) = channel::<()>();

    SENDER.with(|s| {
        *s.borrow_mut() = Some(sender);
    });

    HANDLE.with(|s| {
        let h = thread::spawn(move || loop {
            match receiver.try_recv() {
                Ok(_) | Err(TryRecvError::Disconnected) => {
                    info!("Terminating profiler thread.");
                    break;
                }
                Err(TryRecvError::Empty) => {
                    report_profile_sample();
                    thread::sleep(Duration::from_millis(10)); // 100 Hz
                }
            }
        });
        *s.borrow_mut() = Some(h);
    });

    PROFILER_RUNNING.store(true, Ordering::SeqCst);
}

/// Report profiler sample.
fn report_profile_sample() {
    let profiler_state = PROFILER_STATE.with(|s| *s.borrow());

    let mut hasher = DefaultHasher::new();
    profiler_state.hash(&mut hasher);
    let mut h = hasher.finish() as usize % (PROFILE_HASH_SIZE - 1);

    let profile_samples = PROFILE_SAMPLES.lock().unwrap();
    let mut count = 0;

    while count < PROFILE_HASH_SIZE
        && profile_samples[h].profiler_state.load(Ordering::SeqCst) != profiler_state
        && profile_samples[h].profiler_state.load(Ordering::SeqCst) != 0
    {
        // Wrap around to the start if we hit the end.
        h = (h + 1) % PROFILE_HASH_SIZE;
        count += 1;
    }

    assert!(count != PROFILE_HASH_SIZE, "Profiler hash table filled up!");

    profile_samples[h]
        .profiler_state
        .store(profiler_state, Ordering::SeqCst);
    profile_samples[h].count.fetch_add(1, Ordering::SeqCst);
}

/// Clear profiler.
pub fn clear_profiler() {
    let mut samples = PROFILE_SAMPLES.lock().unwrap();
    for ps in samples.iter_mut() {
        ps.profiler_state.store(0, Ordering::SeqCst);
        ps.count.store(0, Ordering::SeqCst);
    }
}

/// Cleanup profiler thread.
pub fn cleanup_profiler() {
    if !PROFILER_RUNNING.load(Ordering::SeqCst) {
        panic!("Profiler is not running!");
    }

    SENDER.with(|s| {
        if let Some(sender) = s.borrow().as_ref() {
            sender.send(()).unwrap();
        }
    });

    HANDLE.with(|s| {
        let handle = s.borrow_mut().take();
        match handle {
            Some(h) => match h.join() {
                Ok(_) => info!("Profiler thread finished."),
                Err(e) => error!("Profiler thread didn't join: {:?}.", e),
            },
            _ => {}
        }
    });

    PROFILER_RUNNING.store(false, Ordering::SeqCst);
}

/// Convert a `std::time::Instant` to a human readable `String`.
///
/// * `pct` - Percent value in [0, 100].
/// * `now` - Instant value to convert.
fn time_string(pct: f32, now: Instant) -> String {
    let start = PROFILE_START_TIME.with(|s| *s.borrow());
    let pct = pct / 100.0; // remap passed value to to [0,1]
    let ns = (now - start).as_nanos();
    // milliseconds for this category
    let mut ms = (ns as f32 * pct / 1000000.0) as u64;
    // Peel off hours, minutes, seconds, and remaining milliseconds.
    let h = ms / (3600 * 1000);
    ms -= h * 3600 * 1000;
    let m = ms / (60 * 1000);
    ms -= m * (60 * 1000);
    let s = ms / 1000;
    ms -= s * 1000;
    ms /= 10; // only printing 2 digits of fractional seconds
    format!("{:4}:{:02}:{:02}.{:02}", h, m, s, ms).to_string()
}

/// Print profiler results.
pub fn report_profiler_results() {
    let now = Instant::now();

    let mut overall_count = 0_usize;
    let mut used = 0_usize;

    let profile_samples = PROFILE_SAMPLES.lock().unwrap();
    for ps in profile_samples.iter() {
        let count = ps.count.load(Ordering::SeqCst) as usize;
        if count > 0 {
            overall_count += count;
            used += 1;
        }
    }

    println!("Used: {} / {}  entries in profiler hash table", used, PROFILE_HASH_SIZE);

    let mut flat_results: HashMap<String, u64> = HashMap::new();
    let mut hierarchical_results: HashMap<String, u64> = HashMap::new();

    for ps in profile_samples.iter() {
        let count = ps.count.load(Ordering::SeqCst);
        let state = ps.profiler_state.load(Ordering::SeqCst);
        if count == 0 {
            continue;
        }

        let mut s = String::new();
        for b in 0..NUM_PROF_CATEGORIES {
            if state & (1_u64 << b) > 0 {
                if s.len() > 0 {
                    // contribute to the parents...
                    if let Some(r) = hierarchical_results.get_mut(&s) {
                        *r += count;
                    } else {
                        hierarchical_results.insert(s.clone(), count);
                    }
                    s = format!("{}/", s);
                }
                s = format!("{}{}", s, PROF_NAMES[b]);
            }
        }

        if let Some(r) = hierarchical_results.get_mut(&s) {
            *r += count;
        } else {
            hierarchical_results.insert(s, count);
        }

        let name_index = state.log2int();
        if name_index < NUM_PROF_CATEGORIES as i64 {
            if name_index >= 0 {
                if let Some(r) = flat_results.get_mut(PROF_NAMES[name_index as usize]) {
                    *r += count;
                } else {
                    flat_results.insert(PROF_NAMES[name_index as usize].to_string(), count);
                }
            }
        } else {
            error!(
                "Problem getting profile category. ps.profiler.state = {}, name_index = {}",
                state, name_index
            );
        }
    }

    println!("  Profile");
    for (k, v) in hierarchical_results.into_iter() {
        let pct = (100.0 * v as f32) / overall_count as f32;
        let mut indent = 4_usize;
        let mut slash_index = 0_usize;
        if let Some(idx) = k.rfind("/") {
            let n = k.matches("/").count();
            indent += 2 * n;
            slash_index = idx + 1;
        }

        let to_print = &k[slash_index..];
        let indent_rt = max(0, 67 - to_print.len() - indent);
        println!(
            "{pad:>left$}{s}{pad:>right$} {p:5.2}% ({t})",
            pad = "",
            left = indent,
            s = to_print,
            right = indent_rt,
            p = pct,
            t = time_string(pct, now)
        );
    }

    // Sort the flattened ones by time, longest to shortest.
    let flat_vec: Vec<(String, u64)> = flat_results
        .into_iter()
        .sorted_by(|(_, a), (_, b)| Ord::cmp(a, b))
        .collect();

    println!("  Profile (flattened)");
    for (k, v) in flat_vec.into_iter() {
        let pct = (100.0 * v as f32) / overall_count as f32;
        let indent = 4_usize;
        let to_print = &k;
        let indent_rt = max(0, 67 - to_print.len() - indent);
        println!(
            "{pad:>left$}{s}{pad:>right$} {p:5.2}% ({t})",
            pad = "",
            left = indent,
            s = to_print,
            right = indent_rt,
            p = pct,
            t = time_string(pct, now)
        );
    }
    println!();
}
