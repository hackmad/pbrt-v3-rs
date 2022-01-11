# Physically Based Rendering in Rust

The motivation is to explore the algorithms outlined in the
[book](http://www.pbr-book.org/) while simultaneously learning a new language
like Rust.

## Releases

Completed work from the book will be tagged as a release which can be found
[here](https://github.com/hackmad/pbr_rust/releases).

## Renders

Some scenes include PBRT files from:
- [Scenes for pbrt-v3](https://www.pbrt.org/scenes-v3) under `../pbrt-v3-scenes` 
  relative to this repositories root. 
- [Vripped reconstruction of Dragon Model](http://graphics.stanford.edu/data/3Dscanrep/) 
  under `../dragon_recon` relative to this repositories root. 

The relative paths for `Texture`, `Include` etc are relative to the scene file location.

PNG files can be compressed using `pngquant`:

```
pngquant --ext .png --force renders/shapes/sphere.png
```

OR

```
pngquant --ext .png --force renders/shapes/*.png 
```

### Shapes

<a href="renders/shapes/sphere.png"><img src="renders/shapes/sphere.png" style="height: 100px"/></a>
<a href="renders/shapes/cylinder.png"><img src="renders/shapes/cylinder.png" style="height: 100px;"/></a>
<a href="renders/shapes/all-shapes.png"><img src="renders/shapes/all-shapes.png" style="height: 100px;"/></a>
<a href="renders/shapes/disk.png"><img src="renders/shapes/disk.png" style="height: 100px;"/></a>
<a href="renders/shapes/other-quadrics.png"><img src="renders/shapes/other-quadrics.png" style="height: 100px;"/></a>
<a href="renders/shapes/plymesh.png"><img src="renders/shapes/plymesh.png" style="height: 100px;"/></a>
<a href="renders/shapes/loopsubdiv.png"><img src="renders/shapes/loopsubdiv.png" style="height: 100px;"/></a>

### Textures

<a href="renders/textures/2d-mappings.png"><img src="renders/textures/2d-mappings.png" style="height: 100px"/></a>
<a href="renders/textures/trilinear-filtering.png"><img src="renders/textures/trilinear-filtering.png" style="height: 100px"/></a><a 
   href="renders/textures/ewa-filtering.png"><img src="renders/textures/ewa-filtering.png" style="height: 100px"/></a>
<a href="renders/textures/uv.png"><img src="renders/textures/uv.png" style="height: 100px"/></a>
<a href="renders/textures/2d-checkerboard.png"><img src="renders/textures/2d-checkerboard.png" style="height: 100px"/></a>
<a href="renders/textures/dots.png"><img src="renders/textures/dots.png" style="height: 100px"/></a>
<a href="renders/textures/3d-checkerboard.png"><img src="renders/textures/3d-checkerboard.png" style="height: 100px"/></a>

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
