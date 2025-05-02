# METADA (Modern Data Assimilation System for Earth Sciences)

### Documentation Status:
[![Documentation Status](https://readthedocs.org/projects/modern-data-assimilation-system-for-earth-sciences/badge/?version=latest)](https://modern-data-assimilation-system-for-earth-sciences.readthedocs.io/en/latest/?badge=latest)

### Continuous Integration:
| Platform      | Status |
| ------------- | ------ |
| Linux-GNU     | ![Linux-GNU Build](https://github.com/open-x-da/metada/workflows/Linux-GNU/badge.svg) |
| Linux-Intel   | [![Debug+Release](https://img.shields.io/github/actions/workflow/status/open-x-da/metada/intel.yaml?label=Debug%2BRelease)](https://github.com/open-x-da/metada/actions/workflows/intel.yaml) |
| macOS-Clang   | ![macOS-Clang Build](https://github.com/open-x-da/metada/workflows/macOS-Clang/badge.svg) |
| Windows-GNU   | [![Debug+Release](https://img.shields.io/github/actions/workflow/status/open-x-da/metada/windows-gnu.yaml?label=Debug%2BRelease)](https://github.com/open-x-da/metada/actions/workflows/windows-gnu.yaml) |
| Windows-Clang | ![Windows-Clang Build](https://github.com/open-x-da/metada/workflows/Windows-Clang/badge.svg) |
| Code Coverage | [![codecov](https://codecov.io/gh/open-x-da/metada/graph/badge.svg?token=QVL2X0P6UO)](https://codecov.io/gh/open-x-da/metada) |

### License:
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

(C) Copyright 2024 METADA Contributors

### Description:

METADA is a high-performance, cross-platform data assimilation framework designed for Earth science applications. Built with modern C++ and optimized for performance, METADA provides a robust infrastructure for implementing, testing, and deploying data assimilation algorithms across diverse computing environments.

The system features:
- Cross-platform compatibility (Linux, macOS, Windows)
- Performance-optimized implementations of state-of-the-art data assimilation algorithms
- Flexible interfaces for integrating various Earth system models
- Optional GPU acceleration for compute-intensive operations
- Comprehensive testing infrastructure with optional code coverage analysis
- Modern CMake-based build system with cross-compiler support

### Key Features:

- **Advanced DA Algorithms**: Implementation of Ensemble Kalman Filter (EnKF), 4D-Var, Particle Filters, and hybrid methods
- **High Performance Computing**: 
  - Multi-threaded parallelism for shared-memory systems
  - Optional GPU acceleration via CUDA
  - Optimized linear algebra operations using modern C++ libraries
- **Cross-Platform**: Fully tested on Linux, macOS, and Windows (with MinGW 13.2)
- **Flexible Model Integration**: 
  - Generic interfaces for Earth system models
  - Support for structured and unstructured grids
  - Configurable observation operators
- **Modern Software Design**: 
  - C++17 with modular architecture
  - Template metaprogramming for compile-time optimizations
  - Policy-based design for algorithm customization
- **Python Integration**: Comprehensive Python bindings using pybind11
- **Containerization**: Ready-to-use Docker and Singularity containers
- **Quality Assurance**:
  - Comprehensive unit and integration testing
  - Automated CI/CD pipeline
  - Optional code coverage analysis (when lcov/genhtml are available)
- **Dependency Management**: Integrated CMake-based dependency handling

### Installation:

For detailed installation instructions, please visit our [documentation](https://modern-data-assimilation-system-for-earth-sciences.readthedocs.io/en/latest/index.html).

### Development:

- **Compiler Support**: MinGW GCC 14.2, Clang, GCC on Linux/macOS
- **Coverage Testing**: Automatically enabled when lcov/genhtml are available
- **CI/CD**: Automated testing on Linux, macOS, and Windows
