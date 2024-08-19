# Physically Based Rendering in Rust

The motivation is to explore the algorithms outlined in the [book](http://www.pbr-book.org/) while simultaneously
learning a new language like Rust.

For details on how to use this project please refer to [USAGE.md](USAGE.md).

![GUI](renders/GUI.gif)


## [Scenes for pbrt-v3](https://www.pbrt.org/scenes-v3)

### Barcelona Pavillion

![barcelona-pavilion](renders/pbrt-v3-scenes/barcelona-pavilion/pavilion-day.png)

### Human head model with realistic BSSRDF

![head](renders/pbrt-v3-scenes/head.png)

### Landscape

![landscape-view-0](renders/pbrt-v3-scenes/landscape/view-0.png)
![landscape-view-1](renders/pbrt-v3-scenes/landscape/view-1.png)
![landscape-view-2](renders/pbrt-v3-scenes/landscape/view-2.png)
![landscape-view-3](renders/pbrt-v3-scenes/landscape/view-3.png)
![landscape-view-4](renders/pbrt-v3-scenes/landscape/view-4.png)

### Chopper Titan

![chopper-titan](renders/pbrt-v3-scenes/chopper-titan.png)

### BMW M6

![bmw-m6](renders/pbrt-v3-scenes/bmw-m6.png)

### Ecosys

![ecosys](renders/pbrt-v3-scenes/ecosys.png)

### Bathroom

![bathroom](renders/pbrt-v3-scenes/bathroom.png)

### Dragon model rendered with subsurface scattering

![sssdragon](renders/pbrt-v3-scenes/sssdragon/f15-7.png)

### A splash of coffee in a cup with a spoon

![coffee-splash](renders/pbrt-v3-scenes/coffee-splash.png)

### Indoor scene with chairs around a table

![breakfast-bdpt](renders/integrators/bdpt/breakfast.png)

### Stanford Bunny with fur growing out of it

![bunny-fur](renders/pbrt-v3-scenes/bunny-fur/f3-15.png)

### Ganesha

![ganesha](renders/pbrt-v3-scenes/ganesha.png)

### Detailed model of the Austrian Imperial Crown

![crown](renders/pbrt-v3-scenes/crown.png)

### A complex model inspired by a hotel in San Miguel de Allende, Mexico

![sanmiguel pathtraced](renders/pbrt-v3-scenes/sanmiguel.path.png)

### Smoke simulation

![smoke-plume](renders/pbrt-v3-scenes/smoke-plume/plume-284.png)

### Transparent Machines

![transparent machines 542](renders/pbrt-v3-scenes/transparent-machines/frame542.png)
![transparent machines 888](renders/pbrt-v3-scenes/transparent-machines/frame888.png)

### Modern indoor environment

![villa-daylight](renders/pbrt-v3-scenes/villa/villa-daylight.png)

### Caustic glass

![SPPM 10 iterations](renders/integrators/sppm/f16-9a.png)
![SPPM 100 iterations](renders/integrators/sppm/f16-9b.png)

### A glass sphere in participating media

![BDPT](renders/integrators/bdpt/f16-22a.png)
![MLT](renders/integrators/mlt/f16-22b.png)

## Other scenes and figures

### Shapes

![All Shapes](renders/shapes/all-shapes.png)
![Sphere](renders/shapes/sphere.png)
![Cyinder](renders/shapes/cylinder.png)
![Disk](renders/shapes/disk.png)
![Other Quadrics](renders/shapes/other-quadrics.png)
![PLY Mesh](renders/shapes/plymesh.png)
![Loop Subdivision Surface](renders/shapes/loopsubdiv.png)
![Triangles Alpha Mask](renders/shapes/triangles-alpha-mask.png)

### Textures

![2D Mappings](renders/textures/2d-mappings.png)
![UV](renders/textures/uv.png)
![2D Checkerboard](renders/textures/2d-checkerboard.png)
![Dots](renders/textures/dots.png)
![Wrinkled](renders/textures/wrinkled.png)
![Windy](renders/textures/windy.png)
![fBm](renders/textures/fbm.png)
![Marble](renders/textures/marble.png)
![Bilerp](renders/textures/bilerp.png)
![Constant](renders/textures/constant.png)
![Mix](renders/textures/mix.png)
![Scale](renders/textures/scale.png)
![3D Checkerboard](renders/textures/3d-checkerboard.png)
![Trilinear Filtering](renders/textures/trilinear-filtering.png)
![EWA Filtering](renders/textures/ewa-filtering.png)

### Materials

![Matte](renders/materials/matte.png)
![Glass](renders/materials/glass.png)
![Plastic](renders/materials/plastic.png)
![Fourier](renders/materials/fourier.png)
![Mirror](renders/materials/mirror.png)
![Uber](renders/materials/uber.png)
![Metal](renders/materials/metal.png)
![Translucent](renders/materials/translucent.png)
![Substrate](renders/materials/substrate.png)
![Bump Map](renders/materials/bump.png)
![Subsurface](renders/materials/subsurface.png)

### Media

![Grid Density - smoke](renders/media/smoke.png)
![Grid Density - cloud](renders/media/cloud.png)
![Homogeneous - spotfog](renders/media/spotfog.png)

### Cameras

![Perspective](renders/cameras/perspective.png)
![Orthographic](renders/cameras/orthographic.png)
![Realistic](renders/cameras/realistic.png)
![Environment](renders/cameras/environment.png)
![Depth of field](renders/cameras/depth-of-field.png)

### Lights

![Point](renders/lights/point.png)
![Diffuse](renders/lights/diffuse.png)
![Distant](renders/lights/distant.png)
![Infinite No Map](renders/lights/infinite-no-map.png)
![Infinite With Map](renders/lights/infinite-with-map.png)
![Spot](renders/lights/spot.png)
![Projection](renders/lights/projection.png)
![Goniophotometric](renders/lights/goniometric.png)

### Transforms &amp; Object Instancing

![Animation](renders/transforms/anim-bluespheres.png)
![Instances](renders/objects/instances.png)

### Samplers

![(0-2) Sequence](renders/samplers/02sequence.png)
![Halton](renders/samplers/halton.png)
![Maximized Minimal Distance](renders/samplers/maxmindist.png)
![Random](renders/samplers/random.png)
![Stratified](renders/samplers/stratified.png)
