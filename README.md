# Physically Based Rendering in Rust

The motivation is to explore the algorithms outlined in the
[book](http://www.pbr-book.org/) while simultaneously learning a new language
like Rust.

## Releases

Completed work from the book will be tagged as a release which can be found
[here](https://github.com/hackmad/pbr_rust/releases).

## Renders

The images shown below are based on those versions. Not all scenes are
available in each release and may look different due to changes in the
algorithms.

### Shapes

<table>
  <tr>
    <td><a href="renders/shapes/sphere.png"><img src="renders/shapes/sphere.png" width="300px"/></a></td>
    <td><a href="renders/shapes/all-shapes.png"><img src="renders/shapes/all-shapes.png" width="300px"/></a></td>
  </tr>
</table>


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

## Profiling / Performance

### DHAT

Use `--features dhat-rs` to get heap profiling stats. Note that this will be a
lot slower to run.

```
cargo run --release --features dhat-rs -- <input file>
```

This will generate a file `dhat-heap.json` which can be viewed using the DHAT
viewer. There is an [online tool](https://nnethercote.github.io/dh_view/dh_view.html) 
available as well.

### Jemalloc

Use `--features jemalloc` to use jemalloc on Linux/MacOS. On Windows, it will 
use default global allocator. This is mutually exclusive with `dhat-rs` feature.

```
cargo run --release --features jemalloc -- <input file>
```

### Profile guided optimization

See [profile guided optimization](https://doc.rust-lang.org/rustc/profile-guided-optimization.html#a-complete-cargo-workflow).

__NOTE:__ Use of `asdf` is not symlinking toolchain components. Hence the export.

```
# STEP 0: Make sure there is no left-over profiling data from previous runs
rm -rf /tmp/pgo-data

export PATH=~/.asdf/installs/rust/1.57.0/toolchains/1.57.0-x86_64-apple-darwin/lib/rustlib/x86_64-apple-darwin/bin/:$PATH

# STEP 1: Build the instrumented binaries
RUSTFLAGS="-Cprofile-generate=/tmp/pgo-data" cargo build --release --target=x86_64-apple-darwin

# STEP 2: Run the instrumented binaries with some typical data
./target/x86_64-apple-darwin/release/pbr-rust scene1.pbrt
./target/x86_64-apple-darwin/release/pbr-rust scene2.pbrt

# STEP 3: Merge the `.profraw` files into a `.profdata` file
llvm-profdata merge -o /tmp/pgo-data/merged.profdata /tmp/pgo-data

# STEP 4: Use the `.profdata` file for guiding optimizations
RUSTFLAGS="-Cprofile-use=/tmp/pgo-data/merged.profdata" \
    cargo build --release --target=x86_64-apple-darwin
```

### Flamegraph

```
cargo install flamegraph
```

On MacOS DTrace needs root permissions!

```
sudo cargo flamegraph --dev -- target/release/pbr-rust scene.pbrt
```

Ignore this error if you see `cargo-flamegraph.stacks`.

```
[2021-12-25T21:36:44Z ERROR pbr_rust] Error reading file 'target/release/pbr-rust': stream did not contain valid UTF-8
```

The result should be in `flamegraph.svg`.
