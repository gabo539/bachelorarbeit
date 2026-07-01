# Simplicial Approximations of 3D Hypersurfaces in 4D


This repository is a **fork of [OpenVolumeMesh](https://openvolumemesh.org/)**; all
thesis code lives in `examples/bachelor_thesis/`.

## Dependencies

- A C++17 compiler (GCC ≥ 7, Clang ≥ 5, or MSVC 2017+)
- CMake ≥ 3.14
- [Eigen 3](https://eigen.tuxfamily.org) (system package)
- [Polyscope](https://polyscope.run) — fetched automatically by CMake, so **an
  internet connection is required at configure time**
- OpenGL, used by Polyscope (bundled with macOS and Windows; on Linux install the
  X11/OpenGL dev packages)

Install the system dependencies:

| OS | Command |
|----|---------|
| **macOS** | `brew install cmake eigen` |
| **Linux** (Debian/Ubuntu) | `sudo apt install cmake g++ libeigen3-dev libgl1-mesa-dev xorg-dev` |
| **Windows** | Install CMake + Visual Studio (MSVC); `vcpkg install eigen3` |

## Build

From the repository root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This builds OpenVolumeMesh and, as part of its examples, all thesis executables.
The binaries land in `build/Build/bin/` (with MSVC's multi-config generator:
`build/Build/bin/Release/`).

## Run

Each generator prompts for its parameters on the console.

- **macOS / Linux:** `./build/Build/bin/generate_ditorus`
- **Windows:** `build\Build\bin\generate_ditorus.exe`

**Normal vs. debug mode** — by default a generator just builds the shape and opens
the viewer. Pass `--debug` (or `-d`) to additionally run the topology/Betti/quality
checks and write the quality-metric CSVs to `validation/results/`:

```bash
./build/Build/bin/generate_ditorus --debug
```

**Slice widget** — a slice (clipping) plane is created on startup; enable it and
drag its on-screen gizmo from the **Slice Planes** section of Polyscope's control
window (top-left) to cut into the solid and inspect its interior.

## Files

Interactive generators (open a Polyscope window):

- `executables/generate_spheritorus.cc` — builds the spheritorus S²×S¹ in R⁴ and shows it (`--debug` adds all checks).
- `executables/generate_ditorus.cc` — builds the ditorus S¹×S¹×S¹ in R⁴ and shows it (`--debug` adds all checks).
- `executables/generate_cylinder.cc` — the open cylinder S¹×I in R³.
- `executables/generate_prism.cc` —  one triangle×I prism split into 3 tets, both orientations shown.
- `executables/generate_adjacent_prisms.cc` —  confirms shared diagonals stay conforming.

Analysis programs (write CSVs, no window):

- `validation/mesh_quality_assessment.cc` — the three mesh-quality experiments (subdivision, base resolution, radius ratio).
- `validation/schwarz_lantern.cc` — Schwarz-lantern volume-convergence experiment (idea discarded)
- `validation/timing.cc` —  runtime benchmark.

Pipeline headers (`header/`):

- `ovm_nd_types.hh` — the Mesh3D/4D/5D/6D and vector type aliases.
- `base_shape_generator.hh` — base factors: circle (S¹ n-gon), interval I, octahedron (S²).
- `cartesian_product.hh` —  Cartesian product of the factor complexes.
- `nd_subdivision.hh` — prism→tet split, quad→triangle split, barycentric and uniform 1→8 subdivision, factor-block normalization.
- `orientation.hh` — consistent orientation via BFS over the dual graph, then rebuild.
- `projection.hh` — maps the R⁵/R⁶ product down to the R⁴ spheritorus/ditorus embedding.
- `homology.hh` — Betti numbers from boundary-matrix ranks 

Validation headers (`validation/`):

- `validation.hh` — topology tests (Euler, closedness, links) and geometry tests.
- `quality_metrics.hh` — mean ratio, dihedral angles, volume and vertex degree → CSV.
- `visualize_4d.hh` — Polyscope view of a 4D tetrahedral mesh (one axis dropped to color).
- `results_dir.hh` — compile-time base directory for the CSV output.

Build:

- `CMakeLists.txt` — declares all of the targets above.
