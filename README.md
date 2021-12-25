# Physically Based Rendering in Rust

The motivation is to explore the algorithms outlined in the
[book](http://www.pbr-book.org/) while simultaneously learning a new language
like Rust.

## Releases

Completed work from the book will be tagged as a release which can be found
[here](https://github.com/hackmad/pbr_rust/releases).

The images shown below are based on those versions. Not all scenes are
available in each release and may look different due to changes in the
algorithms.

## Building

Build debug profile. The executable will be `target/debug/pbr_rust`.

```bash
cargo build
```

Use `--release` when building/running for faster executable. The executable
will be `target/release/pbr-rust`.

```bash
cargo build --release
```

## Testing

Not everything will be unit tested. The goal was to learn about different
techniques used in Rust for unit testing, property based testing and debugging.

The unit tests can be run as follows:

```
cargo test
```

## Running

This section will be updated as new features get added while progressing
through the book.

The debug version can be run as:

```
cargo run -- <input file>
```

The release version can be run as:

```
cargo run --release -- <input file>
```

## Profiling

Use `--features dhat-rs` to get heap profiling stats. Note that this will be a
lot slower to run.

```
cargo run --release --features dhat-rs -- <input file>
```

This will generate a file `dhat-heap.json` which can be viewed using the DHAT
viewer. There is an [online tool](https://nnethercote.github.io/dh_view/dh_view.html) 
available as well.

## Renders

Coming soon...
