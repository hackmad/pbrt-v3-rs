# Physically Based Rendering in Rust

The motivation is to explore the algorithms outlined in the [book](http://www.pbr-book.org/) while simultaneously
learning a new language like Rust.

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

### PBRT v3 Scenes

<a title="barcelona-pavilion" href="renders/pbrt-v3-scenes/barcelona-pavilion/pavilion-day.png"><img src="renders/pbrt-v3-scenes/barcelona-pavilion/pavilion-day.png" style="height: 100px"/></a>
<a title="head" href="renders/pbrt-v3-scenes/head.png"><img src="renders/pbrt-v3-scenes/head.png" style="height: 100px"/></a>
<a title="landscape-view-0" href="renders/pbrt-v3-scenes/landscape/view-0.png"><img src="renders/pbrt-v3-scenes/landscape/view-0.png" style="height: 100px"/></a>
<a title="landscape-view-1" href="renders/pbrt-v3-scenes/landscape/view-1.png"><img src="renders/pbrt-v3-scenes/landscape/view-1.png" style="height: 100px"/></a>
<a title="landscape-view-2" href="renders/pbrt-v3-scenes/landscape/view-2.png"><img src="renders/pbrt-v3-scenes/landscape/view-2.png" style="height: 100px"/></a>
<a title="landscape-view-3" href="renders/pbrt-v3-scenes/landscape/view-3.png"><img src="renders/pbrt-v3-scenes/landscape/view-3.png" style="height: 100px"/></a>
<a title="landscape-view-4" href="renders/pbrt-v3-scenes/landscape/view-4.png"><img src="renders/pbrt-v3-scenes/landscape/view-4.png" style="height: 100px"/></a>
<a title="chopper-titan" href="renders/pbrt-v3-scenes/chopper-titan.png"><img src="renders/pbrt-v3-scenes/chopper-titan.png" style="height: 100px"/></a>
<a title="BMW M6" href="renders/pbrt-v3-scenes/bmw-m6.png"><img src="renders/pbrt-v3-scenes/bmw-m6.png" style="height: 100px"/></a>
<a title="ecosys" href="renders/pbrt-v3-scenes/ecosys.png"><img src="renders/pbrt-v3-scenes/ecosys.png" style="height: 100px"/></a>
<a title="ecosys" href="renders/pbrt-v3-scenes/bathroom.png"><img src="renders/pbrt-v3-scenes/bathroom.png" style="height: 100px"/></a>
<a title="sssdragon" href="renders/pbrt-v3-scenes/sssdragon/f15-7.png"><img src="renders/pbrt-v3-scenes/sssdragon/f15-7.png" style="height: 100px"/></a>
<a title="coffee-splash" href="renders/pbrt-v3-scenes/coffee-splash.png"><img src="renders/pbrt-v3-scenes/coffee-splash.png" style="height: 100px"/></a>
<a title="BDPT breakfast" href="renders/integrators/bdpt/breakfast.png"><img src="renders/integrators/bdpt/breakfast.png" style="height: 100px"/></a>
<a title="bunny-fur" href="renders/pbrt-v3-scenes/bunny-fur/f3-15.png"><img src="renders/pbrt-v3-scenes/bunny-fur/f3-15.png" style="height: 100px"/></a>
<a title="ganesha" href="renders/pbrt-v3-scenes/ganesha.png"><img src="renders/pbrt-v3-scenes/ganesha.png" style="height: 100px"/></a>
<a title="crown" href="renders/pbrt-v3-scenes/crown.png"><img src="renders/pbrt-v3-scenes/crown.png" style="height: 100px"/></a>
<a title="sanmiguel pathtraced" href="renders/pbrt-v3-scenes/sanmiguel.path.png"><img src="renders/pbrt-v3-scenes/sanmiguel.path.png" style="height: 100px"/></a>
<a title="smoke-plume" href="renders/pbrt-v3-scenes/smoke-plume/plume-284.png"><img src="renders/pbrt-v3-scenes/smoke-plume/plume-284.png" style="height: 100px"/></a>
<a title="transparent machines 542" href="renders/pbrt-v3-scenes/transparent-machines/frame542.png"><img src="renders/pbrt-v3-scenes/transparent-machines/frame542.png" style="height: 100px"/></a>
<a title="transparent machines 888" href="renders/pbrt-v3-scenes/transparent-machines/frame888.png"><img src="renders/pbrt-v3-scenes/transparent-machines/frame888.png" style="height: 100px"/></a>
<a title="villa-daylight" href="renders/pbrt-v3-scenes/villa/villa-daylight.png"><img src="renders/pbrt-v3-scenes/villa/villa-daylight.png" style="height: 100px"/></a>

### Shapes

<a title="Sphere" href="renders/shapes/sphere.png"><img src="renders/shapes/sphere.png" style="height: 100px"/></a>
<a title="Cyinder" href="renders/shapes/cylinder.png"><img src="renders/shapes/cylinder.png" style="height: 100px;"/></a>
<a title="Disk" href="renders/shapes/disk.png"><img src="renders/shapes/disk.png" style="height: 100px;"/></a>
<a title="Other Quadrics" href="renders/shapes/other-quadrics.png"><img src="renders/shapes/other-quadrics.png" style="height: 100px;"/></a>
<a title="PLY Mesh" href="renders/shapes/plymesh.png"><img src="renders/shapes/plymesh.png" style="height: 100px;"/></a>
<a title="Loop Subdivision Surface" href="renders/shapes/loopsubdiv.png"><img src="renders/shapes/loopsubdiv.png" style="height: 100px;"/></a>
<a title="Triangles Alpha Mask" href="renders/shapes/triangles-alpha-mask.png"><img src="renders/shapes/triangles-alpha-mask.png" style="height: 100px"/></a>
<a title="All Shapes" href="renders/shapes/all-shapes.png"><img src="renders/shapes/all-shapes.png" style="height: 100px;"/></a>

### Textures

<a title="2D Mappings" href="renders/textures/2d-mappings.png"><img src="renders/textures/2d-mappings.png" style="height: 100px"/></a>
<a title="UV" href="renders/textures/uv.png"><img src="renders/textures/uv.png" style="height: 100px"/></a>
<a title="2D Checkerboard" href="renders/textures/2d-checkerboard.png"><img src="renders/textures/2d-checkerboard.png" style="height: 100px"/></a>
<a title="Dots" href="renders/textures/dots.png"><img src="renders/textures/dots.png" style="height: 100px"/></a>
<a title="Wrinkled" href="renders/textures/wrinkled.png"><img src="renders/textures/wrinkled.png" style="height: 100px"/></a>
<a title="Windy" href="renders/textures/windy.png"><img src="renders/textures/windy.png" style="height: 100px"/></a>
<a title="fBm" href="renders/textures/fbm.png"><img src="renders/textures/fbm.png" style="height: 100px"/></a>
<a title="Marble" href="renders/textures/marble.png"><img src="renders/textures/marble.png" style="height: 100px"/></a>
<a title="Bilerp" href="renders/textures/bilerp.png"><img src="renders/textures/bilerp.png" style="height: 100px"/></a>
<a title="Constant" href="renders/textures/constant.png"><img src="renders/textures/constant.png" style="height: 100px"/></a>
<a title="Mix" href="renders/textures/mix.png"><img src="renders/textures/mix.png" style="height: 100px"/></a>
<a title="Scale" href="renders/textures/scale.png"><img src="renders/textures/scale.png" style="height: 100px"/></a>
<a title="3D Checkerboard" href="renders/textures/3d-checkerboard.png"><img src="renders/textures/3d-checkerboard.png" style="height: 100px"/></a>
<a title="Trilinear Filtering" href="renders/textures/trilinear-filtering.png"><img src="renders/textures/trilinear-filtering.png" style="height: 100px"/></a>
<a title="EWA Filtering" href="renders/textures/ewa-filtering.png"><img src="renders/textures/ewa-filtering.png" style="height: 100px"/></a>

### Materials

<a title="Matte" href="renders/materials/matte.png"><img src="renders/materials/matte.png" style="height: 100px"/></a>
<a title="Glass" href="renders/materials/glass.png"><img src="renders/materials/glass.png" style="height: 100px"/></a>
<a title="Plastic" href="renders/materials/plastic.png"><img src="renders/materials/plastic.png" style="height: 100px"/></a>
<a title="Fourier" href="renders/materials/fourier.png"><img src="renders/materials/fourier.png" style="height: 100px"/></a>
<a title="Mirror" href="renders/materials/mirror.png"><img src="renders/materials/mirror.png" style="height: 100px"/></a>
<a title="Uber" href="renders/materials/uber.png"><img src="renders/materials/uber.png" style="height: 100px"/></a>
<a title="Metal" href="renders/materials/metal.png"><img src="renders/materials/metal.png" style="height: 100px"/></a>
<a title="Translucent" href="renders/materials/translucent.png"><img src="renders/materials/translucent.png" style="height: 100px"/></a>
<a title="Substrate" href="renders/materials/substrate.png"><img src="renders/materials/substrate.png" style="height: 100px"/></a>
<a title="Bump Map" href="renders/materials/bump.png"><img src="renders/materials/bump.png" style="height: 100px"/></a>
<a title="Subsurface" href="renders/materials/subsurface.png"><img src="renders/materials/subsurface.png" style="height: 100px"/></a>

### Media

<a title="Grid Density - smoke" href="renders/media/smoke.png"><img src="renders/media/smoke.png" style="height: 100px"/></a>
<a title="Grid Density - cloud" href="renders/media/cloud.png"><img src="renders/media/cloud.png" style="height: 100px"/></a>
<a title="Homogeneous - spotfog" href="renders/media/spotfog.png"><img src="renders/media/spotfog.png" style="height: 100px"/></a>

### Cameras

<a title="Perspective" href="renders/cameras/perspective.png"><img src="renders/cameras/perspective.png" style="height: 100px"/></a>
<a title="Orthographic" href="renders/cameras/orthographic.png"><img src="renders/cameras/orthographic.png" style="height: 100px"/></a>
<a title="Realistic" href="renders/cameras/realistic.png"><img src="renders/cameras/realistic.png" style="height: 100px"/></a>
<a title="Environment" href="renders/cameras/environment.png"><img src="renders/cameras/environment.png" style="height: 100px"/></a>
<a title="Depth of field" href="renders/cameras/depth-of-field.png"><img src="renders/cameras/depth-of-field.png" style="height: 100px"/></a>

### Lights

<a title="Point" href="renders/lights/point.png"><img src="renders/lights/point.png" style="height: 100px"/></a>
<a title="Diffuse" href="renders/lights/diffuse.png"><img src="renders/lights/diffuse.png" style="height: 100px"/></a>
<a title="Distant" href="renders/lights/distant.png"><img src="renders/lights/distant.png" style="height: 100px"/></a>
<a title="Infinite No Map" href="renders/lights/infinite-no-map.png"><img src="renders/lights/infinite-no-map.png" style="height: 100px"/></a>
<a title="Infinite With Map" href="renders/lights/infinite-with-map.png"><img src="renders/lights/infinite-with-map.png" style="height: 100px"/></a>
<a title="Spot" href="renders/lights/spot.png"><img src="renders/lights/spot.png" style="height: 100px"/></a>
<a title="Projection" href="renders/lights/projection.png"><img src="renders/lights/projection.png" style="height: 100px"/></a>
<a title="Goniophotometric" href="renders/lights/goniometric.png"><img src="renders/lights/goniometric.png" style="height: 100px"/></a>

### Samplers

<a title="(0-2) Sequence" href="renders/samplers/02sequence.png"><img src="renders/samplers/02sequence.png" style="height: 100px"/></a>
<a title="Halton" href="renders/samplers/halton.png"><img src="renders/samplers/halton.png" style="height: 100px"/></a>
<a title="Maximized Minimal Distance" href="renders/samplers/maxmindist.png"><img src="renders/samplers/maxmindist.png" style="height: 100px"/></a>
<a title="Random" href="renders/samplers/random.png"><img src="renders/samplers/random.png" style="height: 100px"/></a>
<a title="Stratified" href="renders/samplers/stratified.png"><img src="renders/samplers/stratified.png" style="height: 100px"/></a>

### Transforms &amp; Object Instancing

<a title="Animation" href="renders/transforms/anim-bluespheres.png"><img src="renders/transforms/anim-bluespheres.png" style="height: 100px"/></a>
<a title="Instances" href="renders/objects/instances.png"><img src="renders/objects/instances.png" style="height: 100px"/></a>

### Integrators

<a title="SPPM 10 iterations caustic-glass" href="renders/integrators/sppm/f16-9a.png"><img src="renders/integrators/sppm/f16-9a.png" style="height: 100px"/></a>
<a title="SPPM 100 iterations caustic-glass" href="renders/integrators/sppm/f16-9b.png"><img src="renders/integrators/sppm/f16-9b.png" style="height: 100px"/></a>
<a title="BDPT caustic-glass" href="renders/integrators/bdpt/f16-22a.png"><img src="renders/integrators/bdpt/f16-22a.png" style="height: 100px"/></a>
<a title="MLT caustic-glass" href="renders/integrators/mlt/f16-22b.png"><img src="renders/integrators/mlt/f16-22b.png" style="height: 100px"/></a>

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

**NOTE:** Use of `asdf` is not symlinking toolchain components. Hence the export.

```bash
# STEP 0: Make sure there is no left-over profiling data from previous runs
rm -rf /tmp/pgo-data

export PATH=~/.asdf/installs/rust/1.57.0/toolchains/1.57.0-x86_64-apple-darwin/lib/rustlib/x86_64-apple-darwin/bin/:$PATH

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
