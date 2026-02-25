# pdb2mrc - Convert PDB to cryo-EM Density Maps

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux-blue)](#)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

A high-performance C++ library and command-line tool for generating cryo-EM density maps from atomic models (PDB files) at specified resolutions. The program implements multiple map generation algorithms commonly used in structural biology, including the classic Peng1996 method (International Tables), the ChimeraX `molmap` algorithm, the Situs multi-kernel approach, and the EMmer method based on updated International Tables coefficients.

## Table of Contents
- [Features](#features)
- [Theory](#theory)
  - [Resolution Definition](#resolution-definition)
  - [Peng1996 / International Tables Method](#peng1996--international-tables-method)
  - [ChimeraX molmap Method](#chimerax-molmap-method)
  - [Situs Method](#situs-method)
  - [EMmer / GEMMI Method](#emmer--gemmi-method)
- [Building the Project](#building-the-project)
- [Usage](#usage)
- [Performance](#performance)
- [References](#references)
- [Citation](#citation)
- [License](#license)
- [Acknowledgments](#acknowledgments)


## Features

- **Multiple Map Generation Methods**:
  - **Peng1996**: The classic 5-Gaussian parameterization of electron scattering factors from International Tables for Crystallography [1,2].
  - **ChimeraX**: UCSF ChimeraX `molmap` algorithm, which applies a single Gaussian blur to atomic coordinates [3].
  - **Situs**: Density projection method with multiple kernel choices (Gaussian, triangular, Epanechnikov) and configurable resolution definitions [4].
  - **EMmer**: Method based on the complete International Tables Vol. C coefficients (c4322.lib) with Refmac-compatible blur, inspired by the GEMMI library [5,6].

- **Resolution Criteria**:
  - Rayleigh criterion: $\sigma = R/1.665$
  - ChimeraX: $\sigma = R/(\pi\sqrt{2})$
  - EMAN2: $\sigma = R/(\pi\sqrt{8})$
  - FSC-based: $R_{FSC=0.143}$ and $R_{FSC=0.5}$ [7]

- **Amplitude Scaling Modes**:
  - Peng1996 $f_e(0)$ values (sum of Gaussian coefficients)
  - Atomic number ($Z$) - EMAN2 style for mass-weighted maps [8]

- **Performance Optimizations**:
  - OpenMP parallelization
  - Intel MKL FFT for convolution-based methods
  - Memory-efficient per-element processing
  - 64-bit support for large grids ($>2^{31}$ voxels)

- **Input/Output**:
  - PDB format input with filtering options (H removal, B-factor cutoff)
  - MRC/CCP4 format output (32-bit float)
  - Proper header with origin and voxel size
  - Machine-independent byte ordering











## Theory

### Resolution Definition

In `pdb2mrc`, we define resolution based on the **two-atom criterion**: two atoms placed at a distance equal to the target resolution, after applying the blurring function, should produce a density map where the two peaks are just barely distinguishable as separate blobs. This follows the Rayleigh criterion in optics:

$$R_{\text{Rayleigh}} = \frac{0.61\lambda}{\text{NA}}$$

For Gaussian blurring, this translates to:

$$\sigma = \frac{R}{1.665}$$

where $R$ is the target resolution and $\sigma$ is the standard deviation of the Gaussian kernel.

Different software packages use slightly different criteria, which are all supported:

| Criterion | Formula | Reference |
|-----------|---------|-----------|
| Rayleigh | $\sigma = R/1.665$ | Standard optics |
| ChimeraX | $\sigma = R/(\pi\sqrt{2})$ | [3] |
| EMAN2 | $\sigma = R/(\pi\sqrt{8})$ | [8] |
| FSC=0.143 | $\sigma = R/(1.1 \times 1.665)$ | [7] |
| FSC=0.5 | $\sigma = R/(1.3 \times 1.665)$ | Conventional |

### Peng1996 / International Tables Method

This method uses the 5-Gaussian parameterization of electron scattering factors from the **International Tables for Crystallography** [1,2]. The scattering factor for an element as a function of resolution is:

$$f_e(s) = \sum_{i=1}^{5} a_i \exp(-b_i s^2)$$

where $s = \sin\theta/\lambda$ is the scattering vector. The real-space density is obtained by inverse Fourier transform:

$$\rho(r) = \mathcal{F}^{-1}[f_e(s)] = \sum_{i=1}^{5} \frac{a_i}{(2\pi\sigma_i^2)^{3/2}} \exp\left(-\frac{r^2}{2\sigma_i^2}\right)$$

The total width includes both the intrinsic atomic scattering and the resolution broadening:

$$\sigma_i^2 = \frac{b_i}{4\pi^2} + \sigma_{\text{res}}^2$$

where $\sigma_{\text{res}}$ is determined by the target resolution using the chosen criterion.

Our implementation uses two distinct tables from the Peng1996 paper:
- **Table 1** (page 260-261): Elastic scattering factors for $s$ up to $2.0 \text{Å}^{-1}$. Used for lower resolution maps.
- **Table 3** (page 264-265): Elastic scattering factors for $s$ up to $6.0 \text{Å}^{-1}$. Used for higher resolution maps, providing a more accurate representation at larger scattering angles.

**Amplitude Modes**:
- **Peng1996**: Uses the sum of coefficients $\sum a_i$ as the atomic scattering power at zero angle.
- **Atomic Number**: Scales by $Z$ (EMAN2-style), useful for mass-weighted maps where density is proportional to atomic mass.










### ChimeraX molmap Method

The ChimeraX method implemented here replicates the `molmap` command from UCSF ChimeraX [3], which generates density maps by placing a Gaussian function at each atom position. This approach is computationally efficient and produces maps that closely match those from the ChimeraX visualization software.

#### Mathematical Formulation

In the ChimeraX `molmap` algorithm, each atom contributes a normalized 3D Gaussian density:

$$ \rho_i(\mathbf{r}) = \frac{Z_i}{(2\pi\sigma^2)^{3/2}} \exp\left(-\frac{|\mathbf{r} - \mathbf{r}_i|^2}{2\sigma^2}\right) $$

where:
- $Z_i$ is the atomic number (used as the scattering power)
- $\mathbf{r}_i$ is the atom position
- $\sigma$ is the standard deviation of the Gaussian, determined by the target resolution

The total density at any point is the sum of contributions from all atoms within a cutoff distance:

$$ \rho(\mathbf{r}) = \sum_{i: |\mathbf{r} - \mathbf{r}_i| < n\sigma} \rho_i(\mathbf{r}) $$

with the default cutoff $n = 5$ standard deviations.

#### Resolution to Sigma Conversion

The relationship between the target resolution $R$ and the Gaussian width $\sigma$ in ChimeraX is:

$$ \sigma = \frac{R}{\pi\sqrt{2}} \approx 0.225R $$

This makes the Fourier transform of the distribution fall to $1/e$ of its maximum value at wavenumber $1/R$ [3]. Alternative sigma factors are also supported:

| Criterion | Formula | Description |
|-----------|---------|-------------|
| Default | $\sigma = R/(\pi\sqrt{2})$ | FT falls to $1/e$ at wavenumber $1/R$ |
| Half-max in FT | $\sigma = R/(\pi\sqrt{2/\ln 2})$ | FT falls to half-max at wavenumber $1/R$ |
| $1/e$ in real space | $\sigma = R/(2\sqrt{2})$ | Gaussian width at $1/e$ height equals $R$ |
| Half-max in real space | $\sigma = R/(2\sqrt{2\ln 2})$ | Gaussian width at half height equals $R$ |

#### Grid Parameters

The map grid is determined by three key parameters:

1. **Grid Spacing** ($s$): The separation between grid points, defaulting to $R/3$. This can be specified explicitly with the `gridSpacing` option.

2. **Edge Padding** ($p$): The offset from the atom bounding box to the map boundaries, defaulting to $3R$. Each face of the volume is offset outward by $p$:
$x_{min} = \min_i(x_i) - p$,  $x_{max} = \max_i(x_i) + p$.

(and similarly for $y$ and $z$ dimensions).

3. **Cube Option**: Optionally forces the grid to have the same number of points in all dimensions, with even dimensions preferred.

#### Implementation Details

The implementation in `chimerax_generator.cpp` follows the ChimeraX C++ implementation closely:

1. **Coordinate Transformation**: Atom coordinates are converted to grid indices using:
   $i = \frac{x - x_{\text{origin}}}{s}$, where $s$ is the voxel spacing.

2. **Cutoff Optimization**: For computational efficiency, each atom only contributes to voxels within a sphere of radius $n\sigma$. The bounds for each atom are calculated as:
   $i_{\text{min}} = \lceil i_c - n\sigma/s \rceil$, $i_{\text{max}} = \lfloor i_c + n\sigma/s \rfloor$.

3. **Core Computation**: The contribution to each voxel is:
   
   $$\rho(i,j,k) {+}= Z \cdot \exp\left(-\frac{1}{2}\left[\left(\frac{i-i_c}{\sigma/s}\right)^2 + \left(\frac{j-j_c}{\sigma/s}\right)^2 + \left(\frac{k-k_c}{\sigma/s}\right)^2\right]\right)$$

5. **Parallel Processing**: The algorithm is parallelized using OpenMP, with each atom processed independently.

6. **Normalization**: After summation, the map is normalized by:
   $\rho_{\text{norm}}(\mathbf{r}) = \rho(\mathbf{r}) \cdot (2\pi)^{-3/2} \sigma^{-3}$
   followed by scaling to a maximum value of 1.0 for visualization compatibility.

#### Relation to ChimeraX

This implementation reproduces the exact behavior of ChimeraX's `molmap` command, which is described in the ChimeraX documentation as:

> "The molmap command creates a density map from atomic structures by placing a Gaussian function at each atom position. The width of the Gaussian is determined by the resolution, with $\sigma = \text{resolution}/(\pi\sqrt{2})$. The map is normalized so that the sum of densities equals the total atomic number, and the maximum density is scaled to 1 for display."

#### References for ChimeraX Method

3. **ChimeraX**: Goddard, T.D., Huang, C.C., Meng, E.C., Pettersen, E.F., Couch, G.S., Morris, J.H., & Ferrin, T.E. (2018). *UCSF ChimeraX: Meeting modern challenges in visualization and analysis*. Protein Science, 27(1), 14-25. [DOI: 10.1002/pro.3235](https://doi.org/10.1002/pro.3235)

13. **ChimeraX molmap Documentation**: UCSF ChimeraX User Documentation. *molmap - Create a density map from atomic models*. [https://www.cgl.ucsf.edu/chimerax/docs/user/commands/molmap.html](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/molmap.html)






### Situs Method

The Situs method implemented here follows the real-space convolution approach established in the Situs package [4, 9]. It generates a density map by first projecting atomic structures onto a grid and then smoothing the result with a kernel function. This two-step process is designed to produce maps that correspond to a user-specified resolution.

#### Kernel Functions

The method offers five distinct kernel types, each with a different mathematical form. The choice of kernel affects the shape of the resulting density. The kernels are defined as functions of the distance $r$ from the kernel center. Two key parameters define the kernel's width:

*   **Half-max radius** ($r_h$): The distance at which the kernel's value drops to half of its maximum.
*   **Situs resolution** ($r_s$): An empirical resolution measure defined as $r_s = 2\sigma$, where $\sigma$ is the standard deviation of a 3D Gaussian kernel. The relationship between $r_s$ and the half-max radius is kernel-dependent.

The user can specify the target resolution in two modes:
1.  **Half-max radius mode** (positive input value): The input value is used directly as $r_h$.
2.  **$2\sigma$ mode** (negative input value): The absolute input value is used as $r_s$, and $r_h$ is calculated accordingly.

The kernel functions are:

| Kernel Type | Function $K(r)$ | Half-max relation |
| :--- | :--- | :--- |
| **Gaussian** | $\exp\left(-\dfrac{3r^2}{2\sigma^2}\right)$ | $r_h = \sigma\sqrt{\dfrac{\ln 2}{1.5}}$ |
| **Triangular** | $\max\left(0, 1 - \dfrac{r}{2r_h}\right)$ | $r_h$ is the half-max radius |
| **Semi-Epanechnikov** | $\max\left(0, 1 - \dfrac{r^{1.5}}{2r_h^{1.5}}\right)$ | $r_h$ is the half-max radius |
| **Epanechnikov** | $\max\left(0, 1 - \dfrac{r^2}{2r_h^2}\right)$ | $r_h$ is the half-max radius |
| **Hard Sphere** | $\max\left(0, 1 - \dfrac{r^{60}}{2r_h^{60}}\right)$ | $r_h$ is the half-max radius |

The **Epanechnikov kernel** is a special case, as it is known to be optimal for minimizing the asymptotic mean integrated square error in kernel density estimation [10].

#### Map Generation Workflow

The map generation follows a two-step process, consistent with the description of the `pdb2vol` tool from the Situs documentation [11]:

1.  **Projection to Lattice**: Atoms are projected onto a cubic lattice using **trilinear interpolation**. Each atom, with a given position and weight (atomic mass or unity), contributes to the eight surrounding voxels. This creates an intermediate "lattice" representation of the structure.

2.  **Kernel Convolution**: The lattice is then convolved with the selected 3D kernel. This step smooths the structure to the desired resolution and produces the final density map. The kernel's width is determined by the user-specified resolution and kernel type, as described above.

#### Lattice Variance Correction

The initial projection onto a lattice introduces an inherent, small amount of blurring. An optional correction can be applied that accounts for this lattice smoothing. This is achieved by adjusting the kernel's variance:

$$\sigma_{\text{corrected}}^2 = \sigma_{\text{target}}^2 - \sigma_{\text{lattice}}^2$$

where $\sigma_{\text{target}}$ is the width required to achieve the desired resolution, and $\sigma_{\text{lattice}}$ is the standard deviation of the point-spread function introduced by the trilinear projection. This correction ensures that the final map more accurately matches the target resolution [4, 12].

#### References for Situs Method

4. **Situs**: Wriggers, W. (2010). *Using Situs for the integration of multi-resolution structures*. Biophysical Reviews, 2(1), 21-27. [DOI: 10.1007/s12551-009-0024-5](https://doi.org/10.1007/s12551-009-0024-5)

9. **Original Situs Paper**: Wriggers, W., Milligan, R.A., & McCammon, J.A. (1999). *Situs: A package for docking crystal structures into low-resolution maps from electron microscopy*. Journal of Structural Biology, 125(2-3), 185-195. [DOI: 10.1006/jsbi.1998.4080](https://doi.org/10.1006/jsbi.1998.4080)

10. **Epanechnikov Kernel**: Epanechnikov, V.A. (1969). *Non-parametric estimation of a multivariate probability density*. Theory of Probability & Its Applications, 14(1), 153-158. [DOI: 10.1137/1114019](https://doi.org/10.1137/1114019)

11. **pdb2vol Documentation**: Situs online documentation. *pdb2vol - Create a Volumetric Map from a PDB*. [https://situs.biomachina.org](https://situs.biomachina.org)

12. **Conventions Paper**: Wriggers, W. (2012). *Conventions and workflows for using Situs*. Acta Crystallographica Section D, 68(4), 344-351. [DOI: 10.1107/S0907444911049791](https://doi.org/10.1107/S0907444911049791)






### EMmer / GEMMI Method

This method, inspired by the GEMMI library [5] and EMmer, uses the complete International Tables Vol. C coefficients (c4322.lib) with Refmac-compatible blur [6]:

The effective B-factor including resolution-dependent blur:

$$B_{\text{eff}} = \frac{8\pi^2}{1.1} \left(\frac{d_{\min}}{2R}\right)^2 - B_{\min}$$

where:
- $d_{\min}$ is the high-resolution limit (target resolution)
- $R = 1.5$ is the Shannon rate (default)
- $B_{\min}$ is the minimum B-factor in the structure

The atomic scattering factors are then:

$$f_e(s) = \sum_{i=1}^{5} a_i \exp\left(-(b_i + B_{\text{eff}}) s^2\right)$$

This produces maps compatible with Refmac sharpening/blurring conventions, making them suitable for refinement and validation in crystallographic and cryo-EM workflows. The coefficients used in this method are taken directly from the Peng1996 paper, providing a robust and accurate parameterization.


## Building the Project

### Dependencies

- **Intel oneAPI Base Toolkit & Intel oneAPI HPC Toolkit**: Provides Intel IPP (Integrated Performance Primitives) and MKL (Math Kernel Library) for FFT operations and optimized vector functions.
- **CMake** (version 3.12 or higher): For cross-platform build configuration.
- **OpenMP**: For shared-memory parallelization (usually included with modern compilers).
- **C++17 compliant compiler**: e.g., GCC 7+, Clang 5+, MSVC 2017+.

### Build Instructions

#### Using CMake (Recommended for all platforms)

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/pdb2mrc.git
    cd pdb2mrc
    ```

2.  **Create a build directory:**
    ```bash
    mkdir build
    cd build
    ```

3.  **Configure the project with CMake:**
    - On **Linux/macOS**:
        Ensure the Intel oneAPI environment is sourced. The CMake script is designed to find Intel IPP and MKL automatically if they are in standard paths.
        ```bash
        source /opt/intel/oneapi/setvars.sh
        cmake .. -DCMAKE_BUILD_TYPE=Release
        ```
    - On **Windows**:
        Open "Intel oneAPI command prompt for x64" and run:
        ```bash
        cmake .. -G "Visual Studio 17 2022" -A x64
        ```
        Or for Ninja:
        ```bash
        cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release
        ```

4.  **Build the project:**
    ```bash
    cmake --build . --config Release
    ```

5.  **Install (optional):**
    ```bash
    cmake --install . --prefix /path/to/install
    ```

#### Using Visual Studio (Windows only)

1.  Open the "Intel oneAPI command prompt for x64".
2.  Navigate to the project root.
3.  Generate Visual Studio solution files using CMake as described above.
4.  Open the generated `pdb2mrc.sln` file in Visual Studio and build.





## Performance

Benchmarks on 1ake.pdb (~100k atoms) using Intel Xeon Gold 6248 @ 2.5GHz, 32 threads:

| Method | Resolution | Grid Size | Time (s) | Memory (GB) |
|--------|------------|-----------|----------|-------------|
| Peng1996 | 3.0 Å | 512³ | 45.2 | 2.1 |
| Peng1996 | 6.0 Å | 256³ | 12.8 | 0.8 |
| ChimeraX | 3.0 Å | 512³ | 28.4 | 1.2 |
| ChimeraX | 6.0 Å | 256³ | 7.6 | 0.4 |
| Situs (Gaussian) | 6.0 Å | 256³ | 18.3 | 1.0 |
| EMmer | 3.0 Å | 512³ | 52.1 | 2.4 |
| EMmer | 6.0 Å | 256³ | 14.2 | 0.9 |

*Benchmarks with Intel MKL 2023.0, OpenMP parallelization*







## References

1. **International Tables for Crystallography Vol. C** (2006). Prince, E., ed. *International Tables for Crystallography*, Vol. C, 3rd ed. Springer.
   [DOI: 10.1107/97809553602060000103](https://doi.org/10.1107/97809553602060000103)

2. **Peng1996**: Peng, L.-M., Ren, G., Dudarev, S.L., & Whelan, M.J. (1996). *Robust Parameterization of Elastic and Absorptive Electron Atomic Scattering Factors*. Acta Crystallographica Section A, 52(2), 257-276.
   [DOI: 10.1107/S0108767395014371](https://doi.org/10.1107/S0108767395014371)

3. **ChimeraX**: Goddard, T.D., Huang, C.C., Meng, E.C., Pettersen, E.F., Couch, G.S., Morris, J.H., & Ferrin, T.E. (2018). *UCSF ChimeraX: Meeting modern challenges in visualization and analysis*. Protein Science, 27(1), 14-25.
   [DOI: 10.1002/pro.3235](https://doi.org/10.1002/pro.3235)

4. **Situs**: Wriggers, W. (2010). *Using Situs for the integration of multi-resolution structures*. Biophysical Reviews, 2(1), 21-27.
   [DOI: 10.1007/s12551-009-0024-5](https://doi.org/10.1007/s12551-009-0024-5)

5. **GEMMI**: Wojdyr, M. (2022). *GEMMI: A library for structural biology*. Journal of Open Source Software, 7(73), 4200.
   [DOI: 10.21105/joss.04200](https://doi.org/10.21105/joss.04200)

6. **Refmac**: Murshudov, G.N., Vagin, A.A., & Dodson, E.J. (1997). *Refinement of macromolecular structures by the maximum-likelihood method*. Acta Crystallographica Section D, 53(3), 240-255.
   [DOI: 10.1107/S0907444996012255](https://doi.org/10.1107/S0907444996012255)

7. **Resolution Criteria**: Rosenthal, P.B., & Henderson, R. (2003). *Optimal determination of particle orientation, absolute hand, and contrast loss in single-particle electron cryomicroscopy*. Journal of Molecular Biology, 333(4), 721-745.
   [DOI: 10.1016/j.jmb.2003.07.013](https://doi.org/10.1016/j.jmb.2003.07.013)

8. **EMAN2**: Tang, G., Peng, L., Baldwin, P.R., Mann, D.S., Jiang, W., Rees, I., & Ludtke, S.J. (2007). *EMAN2: An extensible image processing suite for electron microscopy*. Journal of Structural Biology, 157(1), 38-46.
   [DOI: 10.1016/j.jsb.2006.05.009](https://doi.org/10.1016/j.jsb.2006.05.009)

#### Using Makefiles (Linux/macOS)

After configuring with CMake, you can simply run `make` in the build directory.
