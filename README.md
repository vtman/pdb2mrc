# pdb2mrc - Convert PDB to cryo-EM Density Maps

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux-blue)](https://github.com/yourusername/pdb2mrc)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

A high-performance C++ library and command-line tool for generating cryo-EM density maps from atomic models (PDB files) at specified resolutions. Supports multiple map generation algorithms commonly used in structural biology.

## Features

- **Multiple Map Generation Methods**:
  - **Peng1996**: Classic 5-Gaussian scattering factors from Peng et al. (1996)
  - **ChimeraX**: UCSF ChimeraX `molmap` algorithm with Gaussian blur
  - **Situs**: Multi-kernel density projection (Gaussian, triangular, Epanechnikov)
  - **EMmer**: International Tables Vol. C coefficients with Refmac blur

- **Resolution Criteria**:
  - Rayleigh criterion (default): $\sigma = R/1.665$
  - ChimeraX: $\sigma = R/(\pi\sqrt{2})$
  - EMAN2: $\sigma = R/(\pi\sqrt{8})$
  - FSC-based: $R_{FSC=0.143}$ and $R_{FSC=0.5}$

- **Amplitude Scaling Modes**:
  - Peng1996 $f_e(0)$ values
  - Atomic number (Z) - EMAN2 style
  - ChimeraX molmap scaling
  - EMmer International Tables coefficients

- **Performance Optimizations**:
  - OpenMP parallelization
  - Intel MKL FFT for convolution
  - Memory-efficient per-element processing
  - 64-bit support for large grids

- **Output Format**:
  - MRC/CCP4 format (32-bit float)
  - Proper header with origin and voxel size
  - Machine-independent byte ordering

## Installation

### Prerequisites

- **Compiler**: C++11 compatible (GCC 4.8+, Clang 3.3+, MSVC 2015+)
- **Intel MKL** (2020 or later) or **OpenBLAS** + FFTW (optional)
- **CMake** 3.12+ (for build system)
- **Intel IPP** (optional, for optimized memory operations)

### Building from Source

```bash
git clone https://github.com/yourusername/pdb2mrc.git
cd pdb2mrc
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
