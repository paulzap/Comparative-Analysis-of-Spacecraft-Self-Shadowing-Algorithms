# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Overview

This repository implements comparative analysis of 5 self-shadowing algorithms for spacecraft Solar Radiation Pressure (SRP) modeling. The project consists of:

- **SatSelf-Shadow**: C++ implementation of self-shadowing algorithms with benchmarking
- **SatForm**: Spacecraft geometry generator creating synthetic datasets
- **Python visualization tools**: For 3D rendering and benchmark analysis

## Key Architecture Concepts

### Algorithm Interface Pattern
All self-shadowing algorithms follow a uniform interface:
```cpp
std::vector<int> calculate_labels_<algorithm_name>(
    const std::vector<Triangle>& triangles,
    const std::vector<double>& sun_vector
)
```
Returns: Vector of labels (1 = sunlit, 0 = shadowed) for each triangle.

### Five Core Algorithms
1. **Pairwise Comparison** - O(n²) reference implementation, 100% accurate ground truth
2. **Sweep Line** - O(n log n) spatial decomposition, scales best for large models
3. **Shadow Mapping** - O(n) GPU-inspired, fastest for real-time onboard applications
4. **Uniform Grid** - O(n + k) spatial hashing, optimized for dense geometries
5. **Ray Casting with BVH** - O(n log n) most robust, perfect accuracy with linear scaling

**Critical**: Always benchmark new algorithms against Pairwise Comparison as ground truth.

### Data Flow Architecture
1. SatForm generates parametric spacecraft models → CSV files (vertices + normals)
2. SatelliteDataset loader parses CSV → Triangle structures
3. SatSelf-Shadow runs benchmarks → timing + accuracy metrics → CSV results
4. Python scripts visualize models and performance comparisons

### Coordinate System
- Sun vector defines light direction in 3D space
- Triangles have vertices and outward-facing normals
- 2D projection plane perpendicular to sun vector for most algorithms
- Shadow mapping uses discrete grid representation

## Build and Development Commands

### Windows (Primary Platform)

**Building SatSelf-Shadow:**
```powershell
cd SatSelf-Shadow
# Compile C++ files (main.cpp, ShadowAlgorithms.cpp) with C++17 support
# MSVC: cl /EHsc /std:c++17 /O2 main.cpp ShadowAlgorithms.cpp
# Or use your preferred build system
```

**Building SatForm:**
```powershell
cd SatForm
# Compile C++ files (main.cpp and related headers) with C++17 support
# MSVC: cl /EHsc /std:c++17 /O2 main.cpp
```

### Linux/macOS (CMake - if implemented)

**Building:**
```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

**Note**: Always use Release builds for benchmarking. Debug builds are 10-50x slower and produce misleading performance metrics.

### Running Components

**Generate spacecraft datasets:**
```powershell
cd SatForm
.\main.exe  # or whatever you name the compiled executable
# Outputs to data3d_<scale>/ directories
# Edit main.cpp to configure scale factor, model count, component probabilities
```

**Run self-shadowing benchmarks:**
```powershell
cd SatSelf-Shadow
.\main.exe  # or your executable name
# Processes all datasets in data3d_*/ directories
# Outputs benchmark_results_*.csv with timing and accuracy
```

**Visualize spacecraft models:**
```bash
python SatSelf-Shadow/visualize3d.py --file data3d_0.5/model_0001.csv --show-normals
```

### Python Environment Setup

**Required dependencies:**
```bash
pip install numpy matplotlib pandas
```

## Development Workflow

### Adding New Algorithms

1. **Implement** in `SatSelf-Shadow/ShadowAlgorithms.cpp` following the standard interface
2. **Add header** declaration in `SatSelf-Shadow/ShadowAlgorithms.h`
3. **Register** in benchmark suite in `SatSelf-Shadow/main.cpp`
4. **Document** complexity (time/space) and use cases in code comments
5. **Benchmark** against all datasets (data3d_0.05, data3d_0.5, data3d_1.5)
6. **Validate** accuracy against Pairwise Comparison (target >95% for fast algorithms)
7. **Update** README.md with algorithm description and benchmark results

### Testing Algorithms

**No dedicated test framework** - validation via benchmarking:
- Reference algorithm (Pairwise) provides ground truth labels
- Metrics: Precision, Recall, F1-score computed automatically
- Test across multiple scales: small satellites (0.05m), CubeSats (0.5m), large spacecraft (1.5m)

### Performance Considerations

**Critical performance factors:**
- Spatial acceleration structures (grids, sorting) are essential for O(n log n) performance
- 2D projection reduces 3D shadowing to 2D polygon overlap problem
- Triangle back-face culling: if normal·sun_vector < 0, triangle is already shadowed
- Shadow mapping resolution trade-off: higher resolution = better accuracy, slower computation
- Memory locality matters: process triangles in depth order when possible

**Benchmarking guidelines:**
- Warm up runs before timing measurements
- Average over multiple runs to reduce variance
- Report both mean and standard deviation
- Test on ~1000 triangle models for realistic comparison

## Dataset Specifications

### CSV Format (SatForm Output)
```
v1_x, v1_y, v1_z, v2_x, v2_y, v2_z, v3_x, v3_y, v3_z, n_x, n_y, n_z
```
- Rows: One per triangle
- Columns: 3 vertices (9 floats) + 1 normal (3 floats)
- Precision: Double precision (16 significant digits)
- Normals: Unit vectors pointing outward from spacecraft

### Dataset Scales
- **data3d_0.05**: ~850 triangles, small satellites
- **data3d_0.5**: ~1200 triangles, CubeSats (most common for benchmarking)
- **data3d_1.5**: ~2100 triangles, large spacecraft

## Code Style Standards

### C++ Conventions (from CONTRIBUTING.md)
- Modern C++17 features encouraged
- Const correctness mandatory for algorithm parameters
- Meaningful variable names (avoid single letters except loop counters)
- Comments required for complex geometric operations
- Avoid premature optimization - profile first

### Python Conventions
- PEP 8 style guide
- Type hints for function signatures
- Docstrings with Args/Returns sections
- Modular functions for reusability

## Common Pitfalls

1. **Floating-point precision**: Use epsilon comparisons for geometric tests (e.g., `abs(dot - 0.0) < 1e-10`)
2. **Edge cases**: Coplanar triangles, degenerate triangles (zero area), grazing angles
3. **Debug vs Release**: Never benchmark Debug builds - results are meaningless
4. **Sun vector normalization**: Ensure sun vector is unit length before algorithm input
5. **Triangle winding order**: Normals must point outward; incorrect winding breaks back-face culling

## File Structure Expectations

```
SatSelf-Shadow/
├── ShadowAlgorithms.h          # Algorithm declarations
├── ShadowAlgorithms.cpp        # Algorithm implementations
├── main.cpp                    # Benchmark driver (main entry point)
├── SatelliteDataset.h          # CSV dataset loader
└── visualize3d.py              # 3D model visualization

SatForm/
├── CLASS_Spacecraft.h          # Top-level spacecraft model
├── CLASS_SpacecraftPart.h      # Components (panels, antenna, body)
├── CLASS_Shape.h               # Geometric primitives (boxes, cylinders)
├── main.cpp                    # Dataset generator main
└── visualize3d.py              # Model visualization
```

## Research Context

This project addresses a critical problem in spacecraft dynamics: self-shadowing affects:
- **Orbit determination**: SRP accelerations depend on illuminated surface area
- **Attitude dynamics**: Asymmetric shadowing creates torques
- **Thermal analysis**: Shadow patterns determine temperature distribution
- **Power systems**: Solar panel output depends on illuminated area

Target use case: **Real-time onboard computation** with constrained computational resources, requiring algorithms that balance accuracy (>95%) with speed (<50ms for 1000 triangles).

## Key Algorithm Selection Guide

- **Ground truth validation**: Use Pairwise Comparison
- **Precise orbit determination & attitude control**: Use Ray Casting with BVH (most robust, perfect accuracy)
- **Real-time or resource-constrained systems**: Use Shadow Mapping (fastest, 95-97% accuracy)
- **Large complex models (>5000 triangles)**: Use Sweep Line (best scaling)
- **Dense geometries**: Use Uniform Grid (optimized for dense configurations)
