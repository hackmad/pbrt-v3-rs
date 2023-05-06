# Usage

## Building

Build debug profile. The executable will be `target/debug/pbrt-v3-rs`.

```bash
cargo build
```

Use `--release` when building/running for faster executable. The executable will be `target/release/pbrt-v3-rs`.

```bash
cargo build --release
```

## Testing

The unit tests can be run as follows:

```bash
cargo test
```

## Running

This section will be updated as new features get added while progressing through the book. The debug version can be run
as:

```bash
cargo run -- [OPTIONS] <FILE1> <FILE2> ...
```

The release version can be run as:

```bash
cargo run --release -- [OPTIONS] <FILE1> <FILE2> ...
```

To run the compiled binary directly use:

```bash
./target/release/pbrt-v3-rs [OPTIONS] <FILE1> <FILE2> ...
```

To run all scenes:

```bash
cargo run --release -- [OPTIONS] <FILE1> <FILE2> ...
```

```bash
./target/releases/pbrt-v3-rs [OPTIONS] $(find scenes | grep -v geometry | grep "\.pbrt$")
```

To run examples from `../pbrt-v3-scenes` you can use relative paths to point to the binary. For example, the example in
`../pbrt-v3-scenes/caustic-glass` can be rendered like this:

```bash
cd ../pbrt-v3-scenes/caustic-glass

../../pbrt-v3-rs/target/release/pbrt-v3-rs -t 4 f16-9a.pbrt
```

## Profiling / Performance

### DHAT

Use `--features dhat-rs` to get heap profiling stats. Note that this will be a
lot slower to run.

```bash
cargo run --release --features dhat-rs -- [OPTIONS] <FILE1> <FILE2> ...
```

This will generate a file `dhat-heap.json` which can be viewed using the DHAT
viewer. There is an [online tool](https://nnethercote.github.io/dh_view/dh_view.html)
available as well.

### Jemalloc

Use `--features jemalloc` to use jemalloc on Linux/MacOS. On Windows, it will use default global allocator. This is
mutually exclusive with `dhat-rs` feature.

```bash
cargo run --release --features jemalloc -- [OPTIONS] <FILE1> <FILE2> ...
```

### Profile guided optimization

See [profile guided optimization](https://doc.rust-lang.org/rustc/profile-guided-optimization.html#a-complete-cargo-workflow).

```bash
# STEP 0: Make sure there is no left-over profiling data from previous runs
rm -rf /tmp/pgo-data

# STEP 1: Build the instrumented binaries
RUSTFLAGS="-Cprofile-generate=/tmp/pgo-data" cargo build --release --target=x86_64-apple-darwin

# STEP 2: Run the instrumented binaries with some typical data
./target/x86_64-apple-darwin/release/pbrt-v3-rs scene1.pbrt
./target/x86_64-apple-darwin/release/pbrt-v3-rs scene2.pbrt

# STEP 3: Merge the `.profraw` files into a `.profdata` file
llvm-profdata merge -o /tmp/pgo-data/merged.profdata /tmp/pgo-data

# STEP 4: Use the `.profdata` file for guiding optimizations
RUSTFLAGS="-Cprofile-use=/tmp/pgo-data/merged.profdata" \
    cargo build --release --target=x86_64-apple-darwin
```

### Flamegraph

```bash
cargo install flamegraph
```

On MacOS DTrace needs root permissions!

```bash
sudo cargo flamegraph --dev -- target/release/pbrt-v3-rs scene.pbrt
```

Ignore this error if you see `cargo-flamegraph.stacks`.

```bash
[2021-12-25T21:36:44Z ERROR pbrt-v3-rs] Error reading file 'target/release/pbrt-v3-rs': stream did not contain valid UTF-8
```

The result should be in `flamegraph.svg`.

## Renders

Some scenes include PBRT files from:

- [Scenes for pbrt-v3](https://www.pbrt.org/scenes-v3) under `../pbrt-v3-scenes` relative to this repositories root.
- [Vripped reconstruction of Dragon Model](http://graphics.stanford.edu/data/3Dscanrep/) under `../dragon_recon`
  relative to this repositories root.

**NOTES:**

- The relative paths for `Texture`, `Include` etc are relative to the scene file location.
- The parser was not ported from the original implementation and it uses [pest](https://github.com/pest-parser/pest).
  There is an [online editor](https://pest.rs/#editor) that is useful for testing grammars.
- Additional scenes adapted from `pbrt-v3-scenes` to make it easier to create test renders:
  - `scenes/materials/*.pbrt`

PNG files can be compressed using `pngquant`:

```bash
pngquant --ext .png --force renders/shapes/sphere.png
```

OR

```bash
pngquant --ext .png --force renders/shapes/*.png
```
