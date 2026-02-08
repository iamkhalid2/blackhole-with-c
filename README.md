# Black Hole Simulation in Pure C

A 3D black hole simulation with gravitational lensing, rendered entirely in C using only standard libraries.

## Features

- **Schwarzschild black hole physics** - Light bending calculated using general relativity approximation
- **Ray marching** - Adaptive step sizes for accurate light path tracing
- **Procedural starfield** - Hash-based star generation (no random library needed)
- **HDR rendering** - ACES tonemapping and gamma correction
- **PPM output** - Pure text image format, no external libraries

## Building

### Using Make (Linux/macOS/MinGW)

```bash
make
```

### Using GCC directly (Windows/any platform)

```bash
gcc -O2 src/main.c -o blackhole -lm
```

## Running

```bash
./blackhole
```

Output will be saved to `output/blackhole.ppm`

## Converting to PNG (optional)

Using ImageMagick:
```bash
magick output/blackhole.ppm output/blackhole.png
```

Or using ffmpeg:
```bash
ffmpeg -i output/blackhole.ppm output/blackhole.png
```

## Project Structure

```
├── src/
│   ├── main.c           # Entry point, render loop
│   ├── vec3.h           # Vector math
│   ├── ray.h            # Ray structure
│   ├── camera.h         # Camera and viewport
│   ├── framebuffer.h    # Image buffer and PPM output
│   ├── blackhole.h      # Black hole physics
│   └── scene.h          # Starfield and accretion disk
├── output/              # Rendered images
├── Makefile
└── README.md
```

## Physics

The simulation uses a simplified Schwarzschild metric approximation:

- **Event horizon**: Rays within Rs = 2GM/c² are absorbed
- **Light bending**: Based on weak-field deflection angle θ ≈ 4GM/(c²b)
- **Ray marching**: Adaptive step sizes near the black hole for accuracy

## Future Additions

- [ ] Accretion disk with temperature gradient
- [ ] Animation / orbit camera
- [ ] Doppler shift effects
- [ ] Kerr (rotating) black hole
