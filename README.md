# pdb2mrc - Convert PDB to cryo-EM Density Maps

A high-performance C++ library and command-line tool for generating cryo-EM density maps from atomic models (PDB files) at specified resolutions. The program implements multiple map generation algorithms commonly used in structural biology, including Peng1996 (International Tables), ChimeraX molmap, Situs, and EMmer methods.

## Table of Contents
- [Theory](#theory)
  - [Resolution Definition](#resolution-definition)
  - [Peng1996 / International Tables Method](#peng1996--international-tables-method)
  - [ChimeraX molmap Method](#chimerax-molmap-method)
  - [Situs Method](#situs-method)
  - [EMmer / GEMMI Method](#emmer--gemmi-method)
- [Installation](#installation)
- [Command Line Usage](#command-line-usage)
- [Examples](#examples)
- [Performance](#performance)
- [References](#references)
- [Citation](#citation)

## Theory

### Resolution Definition

In pdb2mrc, we define resolution based on the **two-atom criterion**: two atoms placed at a distance equal to the target resolution, after applying the blurring function, should produce a density map where the two peaks are just barely distinguishable as separate blobs. This follows the Rayleigh criterion in optics:

<img src="https://latex.codecogs.com/svg.latex?R_{\text{Rayleigh}}&space;=&space;\frac{0.61\lambda}{NA}" title="R_{\text{Rayleigh}} = \frac{0.61\lambda}{NA}" />

For Gaussian blurring, this translates to:

<img src="https://latex.codecogs.com/svg.latex?\sigma&space;=&space;\frac{R}{1.665}" title="\sigma = \frac{R}{1.665}" />

where <i>R</i> is the target resolution and <i>&sigma;</i> is the standard deviation of the Gaussian kernel.

Different software packages use slightly different criteria:

| Criterion | Formula | Reference |
|-----------|---------|-----------|
| Rayleigh | <img src="https://latex.codecogs.com/svg.latex?\sigma&space;=&space;R/1.665"> | Standard optics |
| ChimeraX | <img src="https://latex.codecogs.com/svg.latex?\sigma&space;=&space;R/(\pi\sqrt{2})"> | Goddard et al. (2018) |
| EMAN2 | <img src="https://latex.codecogs.com/svg.latex?\sigma&space;=&space;R/(\pi\sqrt{8})"> | Tang et al. (2007) |
| FSC=0.143 | <img src="https://latex.codecogs.com/svg.latex?\sigma&space;=&space;R/(1.1\times1.665)"> | Rosenthal & Henderson (2003) |
| FSC=0.5 | <img src="https://latex.codecogs.com/svg.latex?\sigma&space;=&space;R/(1.3\times1.665)"> | Conventional |

### Peng1996 / International Tables Method

This method uses the 5-Gaussian parameterization of electron scattering factors from the **International Tables for Crystallography** [1,2]. The scattering factor for an element as a function of resolution is:

<img src="https://latex.codecogs.com/svg.latex?f_e(s)&space;=&space;\sum_{i=1}^{5}&space;a_i&space;\exp(-b_i&space;s^2)" title="f_e(s) = \sum_{i=1}^{5} a_i \exp(-b_i s^2)" />

where <img src="https://latex.codecogs.com/svg.latex?s&space;=&space;\sin\theta/\lambda" title="s = \sin\theta/\lambda" /> is the scattering vector. The real-space density is obtained by inverse Fourier transform:

<img src="https://latex.codecogs.com/svg.latex?\rho(r)&space;=&space;\mathcal{F}^{-1}[f_e(s)]&space;=&space;\sum_{i=1}^{5}&space;\frac{a_i}{(2\pi\sigma_i^2)^{3/2}}&space;\exp\left(-\frac{r^2}{2\sigma_i^2}\right)" title="\rho(r) = \mathcal{F}^{-1}[f_e(s)] = \sum_{i=1}^{5} \frac{a_i}{(2\pi\sigma_i^2)^{3/2}} \exp\left(-\frac{r^2}{2\sigma_i^2}\right)" />

The total width includes both the intrinsic atomic scattering and the resolution broadening:

<img src="https://latex.codecogs.com/svg.latex?\sigma_i^2&space;=&space;\frac{b_i}{4\pi^2}&space;+&space;\sigma_{\text{res}}^2" title="\sigma_i^2 = \frac{b_i}{4\pi^2} + \sigma_{\text{res}}^2" />

where <img src="https://latex.codecogs.com/svg.latex?\sigma_{\text{res}}" title="\sigma_{\text{res}}" /> is determined by the target resolution using the chosen criterion.

**Amplitude Modes**:
- **Peng1996**: Uses the sum of coefficients <img src="https://latex.codecogs.com/svg.latex?\sum&space;a_i" title="\sum a_i" /> as the atomic scattering power at zero angle
- **Atomic Number**: Scales by Z (EMAN2-style), useful for mass-weighted maps

### ChimeraX molmap Method

UCSF ChimeraX implements a simplified approach using a single Gaussian per atom [3]:

<img src="https://latex.codecogs.com/svg.latex?\rho(r)&space;=&space;\frac{Z}{(2\pi\sigma^2)^{3/2}}&space;\exp\left(-\frac{r^2}{2\sigma^2}\right)" title="\rho(r) = \frac{Z}{(2\pi\sigma^2)^{3/2}} \exp\left(-\frac{r^2}{2\sigma^2}\right)" />

with the resolution-dependent width:

<img src="https://latex.codecogs.com/svg.latex?\sigma&space;=&space;\frac{R}{\pi\sqrt{2}}" title="\sigma = \frac{R}{\pi\sqrt{2}}" />

The algorithm uses a cutoff at <i>n&sigma;</i> (default 5.0) for efficiency:

<img src="https://latex.codecogs.com/svg.latex?\rho(r)&space;=&space;0&space;\quad\text{for}\quad&space;r&space;>&space;n\sigma" title="\rho(r) = 0 \quad\text{for}\quad r > n\sigma" />

### Situs Method

The Situs package [4] offers multiple kernel types with flexible resolution definitions. The user can specify resolution either as:
- **Half-max radius** <i>r<sub>h</sub></i> (positive value)
- **2&sigma;** (negative value)

The kernel functions are:

| Kernel Type | Function <i>K(r)</i> | Half-max relation |
|-------------|----------------------|-------------------|
| Gaussian | <img src="https://latex.codecogs.com/svg.latex?\exp(-1.5&space;r^2/\sigma^2)" title="\exp(-1.5 r^2/\sigma^2)" /> | <img src="https://latex.codecogs.com/svg.latex?r_h&space;=&space;\sigma\sqrt{\ln&space;2&space;/&space;1.5}" title="r_h = \sigma\sqrt{\ln 2 / 1.5}" /> |
| Triangular | <img src="https://latex.codecogs.com/svg.latex?1&space;-&space;r/(2r_h)" title="1 - r/(2r_h)" /> | <i>r<sub>h</sub></i> = half-max radius |
| Semi-Epanechnikov | <img src="https://latex.codecogs.com/svg.latex?1&space;-&space;(r/r_h)^{1.5}" title="1 - (r/r_h)^{1.5}" /> | <i>r<sub>h</sub></i> = half-max radius |
| Epanechnikov | <img src="https://latex.codecogs.com/svg.latex?1&space;-&space;(r/r_h)^2" title="1 - (r/r_h)^2" /> | <i>r<sub>h</sub></i> = half-max radius |
| Hard Sphere | <img src="https://latex.codecogs.com/svg.latex?1&space;-&space;(r/r_h)^{60}" title="1 - (r/r_h)^{60}" /> | <i>r<sub>h</sub></i> = half-max radius |

The map generation follows a two-step process:
1. **Projection**: Atoms are projected onto a lattice using trilinear interpolation
2. **Convolution**: The lattice is convolved with the chosen kernel

### EMmer / GEMMI Method

This method, inspired by the GEMMI library [5] and EMmer, uses the complete International Tables Vol. C coefficients (c4322.lib) with Refmac-compatible blur [6]:

The effective B-factor including resolution-dependent blur:

<img src="https://latex.codecogs.com/svg.latex?B_{\text{eff}}&space;=&space;\frac{8\pi^2}{1.1}&space;\left(\frac{d_{\min}}{2R}\right)^2&space;-&space;B_{\min}" title="B_{\text{eff}} = \frac{8\pi^2}{1.1} \left(\frac{d_{\min}}{2R}\right)^2 - B_{\min}" />

where:
- <i>d<sub>min</sub></i> is the high-resolution limit
- <i>R</i> = 1.5 is the Shannon rate (default)
- <i>B<sub>min</sub></i> is the minimum B-factor in the structure

The atomic scattering factors are then:

<img src="https://latex.codecogs.com/svg.latex?f_e(s)&space;=&space;\sum_{i=1}^{5}&space;a_i&space;\exp\left(-(b_i&space;+&space;B_{\text{eff}})&space;s^2\right)" title="f_e(s) = \sum_{i=1}^{5} a_i \exp\left(-(b_i + B_{\text{eff}}) s^2\right)" />

This produces maps compatible with Refmac sharpening/blurring conventions.

## Installation

### Prerequisites

- **Compiler**: C++11 compatible (GCC 4.8+, Clang 3.3+, MSVC 2015+)
- **Intel MKL** (2020 or later) - for FFT operations
- **Intel IPP** (optional, for optimized memory operations)
- **CMake** 3.12+ (for build system)
- **OpenMP** (for parallel processing)

### Building from Source

```bash
git clone https://github.com/yourusername/pdb2mrc.git
cd pdb2mrc
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
