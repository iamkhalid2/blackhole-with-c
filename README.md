# Supermassive Kerr Black Hole Simulation (Pure C)

![Simulation Output](blackhole_animation.gif)
*(Placeholder for output GIF - render using instructions below)*

A scientifically accurate, high-performance 3D black hole simulation written entirely in **pure C** (C99 standard) with no external library dependencies. It implements general relativity ray-tracing to visualize a rotating (Kerr) black hole with an accretion disk, achieving visual fidelity comparable to *Interstellar*'s Gargantua.

## 🚀 Features

*   **Kerr Metric Physics:** Simulates a rotating black hole ($a=0.95$) with frame dragging (Lense-Thirring effect) and an oblate event horizon.
*   **Accurate Accretion Disk:**
    *   **ISCO Inner Edge:** The disk naturally terminates at the Innermost Stable Circular Orbit (~1.94M), calculated via the Bardeen-Press-Teukolsky formula.
    *   **Thermodynamics:** Uses the Shakura-Sunyaev temperature profile ($T \propto r^{-3/4}$) with Planck blackbody color mapping.
*   **Relativistic Optics:**
    *   **Doppler Beaming:** The approaching side of the disk is brighter ($I \propto \delta^3$) due to relativistic beaming.
    *   **Gravitational Redshift:** Light emitted deep in the gravity well is redshifted and dimmed.
    *   **Gravitational Lensing:** Adaptive ray-tracing simulating curved light paths in strong gravity.
*   **Visual Fidelity:**
    *   **Logarithmic Exposure:** Mimics astronomical imaging to resolve both the blinding core and faint outer disk.
    *   **Procedural Starfield:** High-resolution star background.
    *   **Optimized:** Uses OpenMP (optional) and architectural optimizations (`-O3 -march=native`) for rendering.

## 🔬 The Physics

Great care was taken to ensure the simulation adheres to the laws of General Relativity and Thermodynamics.

### 1. The Kerr Metric & Horizon
Unlike a static Schwarzschild black hole, a rotating Kerr black hole drags spacetime with it. The event horizon $r_+$ is oblate:
$$ r_+ = M (1 + \sqrt{1 - a^2}) $$
Where $a$ is the spin parameter (0.95 in this sim).

### 2. Accretion Disk Limits (ISCO)
The accretion disk cannot exist inside the Innermost Stable Circular Orbit. For a prograde orbit in a Kerr metric, the ISCO radius $r_{isco}$ is determined by the root of:
$$ r^2 - 6Mr + 8a\sqrt{Mr} - 3a^2 = 0 $$
For our spin $a=0.95$, this puts the disk edge at $\approx 1.94M$, much closer (and hotter) than the $6.0M$ limit of a non-rotating hole.

### 3. Relativistic Beaming & Redshift
Photons emitted from the disk are shifted in frequency ($\nu$) and intensity ($I$) by the Doppler factor $\delta$:
$$ \delta = \frac{1}{\gamma (1 - \beta \cos \theta)} $$
The observed intensity $I_{obs}$ scales with the third power of the Doppler factor:
$$ I_{obs} = I_{emit} \cdot \delta^3 $$
This is why the left side of the disk (moving towards the camera) is significantly brighter than the right.

### 4. Temperature & Color
Colors are not artistically chosen but calculated. The disk temperature follows the standard thin-disk model:
$$ T(r) \propto r^{-3/4} $$
This temperature is converted to RGB values using the **Planck Law** for blackbody radiation, producing a physically correct gradient from Blue-White (~15,000K) at the center to Red-Orange (~3,000K) at the edge.

## 🛠️ Building & Running

### Requirements
*   **GCC** (or Clang)
*   **Make** (optional)
*   No external libraries required (uses `math.h`, `stdio.h`)

### Compilation
Use the provided `Makefile` or compile manually with optimizations:
```bash
gcc -O3 -march=native src/main.c -o blackhole -lm
```
*(Note: Add `-fopenmp` if your compiler supports it for multithreaded rendering)*

### Running
```bash
./blackhole
```
The simulation will generate frames in the `output/` directory.

### Configuration
Edit `src/main.c` to change settings:
```c
#define WIDTH 1600         // Resolution
#define HEIGHT 900
#define ANIMATE 1          // 1 = Animation, 0 = Single Frame
#define NUM_FRAMES 30      // Number of frames
#define BH_SPIN 0.95       // Black hole spin (0.0 to 0.99)
```

### Creating a GIF/Video
After rendering frames, use **ImageMagick**:
```bash
magick -delay 10 -loop 0 output/frame_*.ppm blackhole_render.gif
```
Or **FFmpeg**:
```bash
ffmpeg -framerate 30 -i output/frame_%04d.ppm -c:v libx264 -pix_fmt yuv420p output.mp4
```

## 📜 License
MIT License. Free to use for educational and scientific visualization purposes.
