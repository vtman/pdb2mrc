# pdb2mrc - Convert PDB to cryo-EM Density Maps

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux-blue)](#)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

A high-performance C++ library and command-line tool for generating cryo-EM density maps from atomic models (PDB files) at specified resolutions. The program implements multiple **real-space** map generation algorithms commonly used in structural biology, including the classic Peng1996 method (International Tables), the ChimeraX `molmap` algorithm, the Situs multi-kernel approach, and the EMmer method based on updated International Tables coefficients.

## Table of Contents
- [Features](#features)
- [Theory](#theory)
  - [Resolution Definition](#resolution-definition)
  - [EMmer / GEMMI Method](#emmer--gemmi-method)
  - [Peng1996 / International Tables Method](#peng1996--international-tables-method)
  - [ChimeraX molmap Method](#chimerax-molmap-method)
  - [Situs Method](#situs-method)
- [Input Parameters](#input-parameters)
- [Examples](#examples)
- [References](#references)
- [Citation](#citation)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Features

- **Multiple Real-Space Map Generation Methods**:
  - **EMmer/GEMMI**: Method based on the complete International Tables Vol. C coefficients (c4322.lib) with Refmac-compatible blur, inspired by the GEMMI library [Wojdyr2022]
  - **Peng1996**: The classic 5-Gaussian parameterization of electron scattering factors from International Tables for Crystallography [Peng1996]
  - **ChimeraX**: UCSF ChimeraX `molmap` algorithm, which applies a single Gaussian blur to atomic coordinates [Goddard2018, Pettersen2020]
  - **Situs**: Density projection method with multiple kernel choices (Gaussian, triangular, Epanechnikov) and configurable resolution definitions [Wriggers2010, Wriggers1999]

- **Resolution Criteria** (all real-space definitions):
  - Rayleigh criterion: $\sigma = R/1.665$ [Rayleigh1879]
  - ChimeraX: $\sigma = R/(\pi\sqrt{2})$ [Goddard2018]
  - EMAN2: $\sigma = R/(\pi\sqrt{8})$ [Tang2007]

- **Amplitude Scaling Modes**:
  - Peng1996 $f_e(0)$ values (sum of Gaussian coefficients) [Peng1996]
  - Atomic number ($Z$) - EMAN2 style for mass-weighted maps [Tang2007]

- **Performance Optimizations**:
  - OpenMP parallelization
  - Intel MKL FFT for convolution-based methods
  - Memory-efficient per-element processing
  - 64-bit support for large grids ($>2^{31}$ voxels)

- **Input/Output**:
  - PDB format input with filtering options (H removal, B-factor cutoff)
  - MRC/CCP4 format output (32-bit float) [Cheng2015]
  - Proper header with origin and voxel size
  - Machine-independent byte ordering

## Theory

### Resolution Definition

In pdb2mrc, we define resolution based on the **two-atom criterion**: two atoms placed at a distance equal to the target resolution, after applying the blurring function, should produce a density map where the two peaks are just barely distinguishable as separate blobs. This follows the Rayleigh criterion in optics [Rayleigh1879]:

$$R_{\text{Rayleigh}} = \frac{0.61\lambda}{\text{NA}}$$

For Gaussian blurring, this translates to:

$$\sigma = \frac{R}{1.665}$$

where $R$ is the target resolution and $\sigma$ is the standard deviation of the Gaussian kernel.

Different software packages use slightly different criteria, all of which are implemented:

| Criterion | Formula | Reference |
|-----------|---------|-----------|
| Rayleigh | $\sigma = R/1.665$ | Rayleigh1879 |
| ChimeraX | $\sigma = R/(\pi\sqrt{2})$ | Goddard2018 |
| EMAN2 | $\sigma = R/(\pi\sqrt{8})$ | Tang2007 |




### Default Method: Peng1996 with Per-Element FFT Convolution

The default map generation method in pdb2mrc combines the accurate Peng1996 scattering factors with an efficient FFT-based convolution approach that processes each chemical element separately. This design choice enables both physical accuracy and computational efficiency.

#### Core Strategy: Decoupling B-factors from Scattering Factors

A key insight in our implementation is the separation of atomic scattering factors from temperature factors (B-factors). The electron scattering factor for an atom can be written as:

$$f_e(s) = \underbrace{\sum_{i=1}^{5} a_i \exp(-b_i s^2)}_{\text{element-specific}} \times \underbrace{\exp(-B s^2 / 4)}_{\text{atom-specific B-factor}}$$

This factorization is crucial because:
- The element-specific part depends only on the chemical element and target resolution
- The B-factor part varies per atom and includes thermal motion

In real space, this factorization becomes even more powerful. The density contribution from an atom is:

$$\rho(\mathbf{r}) = \underbrace{\mathcal{F}^{-1}\left[\sum_{i=1}^{5} a_i \exp(-b_i s^2)\right]}_{\text{element kernel}} * \underbrace{\mathcal{F}^{-1}\left[\exp(-B s^2 / 4)\right]}_{\text{B-factor Gaussian}}$$

where $*$ denotes convolution. This means the total map can be built by:
1. Precomputing an "element kernel" - the real-space density of a single atom of that element at the target resolution
2. Convolving a point-spread function that depends on each atom's B-factor

#### Per-Element Processing Workflow

Our implementation takes advantage of this factorization through a multi-step process:

1. **Element Classification**: All atoms are grouped by their chemical element (C, N, O, etc.)

2. **Profile Precomputation**: For each unique element present, a radial density profile is precomputed using the Peng1996 coefficients at the target resolution:

   $$\rho_{\text{element}}(r) = \sum_{i=1}^{5} \frac{a_i}{(2\pi(\sigma_i^2 + \sigma_{\text{res}}^2))^{3/2}} \exp\left(-\frac{r^2}{2(\sigma_i^2 + \sigma_{\text{res}}^2)}\right)$$

   where $\sigma_i^2 = b_i/(4\pi^2)$ and $\sigma_{\text{res}}$ is determined by the target resolution.

3. **Grid Projection (Real Space)**: Atoms are projected onto the 3D grid using trilinear interpolation, but **without any B-factor blurring**. This creates a "delta-map" where each atom contributes only to its eight neighboring voxels:

   $$\rho_{\text{delta}}(\mathbf{r}) = \sum_j \delta(\mathbf{r} - \mathbf{r}_j) \cdot \text{occupancy}_j$$

   This step is computationally cheap and preserves the exact atomic positions.

4. **Element-Specific Convolution (Fourier Space)**: For each element type, the delta-map is convolved with the precomputed element kernel using FFT:

   $$\rho_{\text{element}}(\mathbf{r}) = \rho_{\text{delta,element}}(\mathbf{r}) * \rho_{\text{element,kernel}}(\mathbf{r})$$

   The convolution is performed efficiently in Fourier space:
   - Forward FFT of the delta-map
   - Forward FFT of the element kernel
   - Complex multiplication in Fourier space
   - Inverse FFT back to real space

5. **Summation**: The per-element maps are summed to produce the final density map.

#### Why This Approach is Efficient

The per-element FFT convolution strategy offers several advantages:

| Aspect | Traditional Approach | pdb2mrc Approach |
|--------|---------------------|------------------|
| **B-factor handling** | Different kernel per atom | Handled in real space, separated from element kernel |
| **FFT operations** | One per atom (intractable) | One per element type (typically 5-10) |
| **Memory usage** | Store all atom contributions | One grid per element type |
| **Parallelization** | Difficult to load balance | Natural per-element parallelism |

For a typical protein with 10,000 atoms and 8 unique elements:
- Traditional approach: 10,000 FFTs (impossible)
- pdb2mrc approach: 8 FFTs (efficient)

#### Handling B-factors

B-factors are incorporated through a two-stage process:

1. **During projection**: Each atom's contribution is scaled by its occupancy and stored in the per-element delta-map

2. **During convolution**: The element kernel already includes the resolution-dependent broadening. B-factors represent additional Gaussian broadening, which in Fourier space corresponds to multiplication by $\exp(-B s^2 / 4)$

However, since B-factors vary per atom, they cannot be incorporated into the element kernel. Instead, they would require per-atom convolution - which is why in the default Peng1996 method, B-factors are **ignored** (or a default value is used) to maintain the efficiency of per-element convolution.

This design choice reflects a fundamental trade-off:
- **Without B-factors**: Fast, element-based FFT convolution (default)
- **With B-factors**: Slower but more accurate for high-resolution work where thermal motion matters

#### Memory and Performance Optimization

The implementation includes several optimizations:

- **Profile truncation**: Element kernels are truncated at a radius where density falls below $10^{-3}$ of the maximum, reducing kernel size
- **Workspace reuse**: FFT buffers and descriptors are allocated once and reused for all elements
- **OpenMP parallelization**: Per-element processing and FFT operations are parallelized where possible
- **64-bit support**: Handles grids larger than $2^{31}$ voxels for high-resolution maps of large complexes

This approach makes it practical to generate high-quality density maps from PDB files on standard workstations, with typical runtimes of seconds to minutes depending on grid size and atom count.






### Kernel Discretization and Normalization

For an element with Peng1996 coefficients, the continuous kernel in real space is:

$$K_{\text{element}}(\mathbf{r}) = \sum_{i=1}^{5} \frac{a_i}{(2\pi\sigma_i^2)^{3/2}} \exp\left(-\frac{|\mathbf{r}|^2}{2\sigma_i^2}\right)$$

where $\sigma_i^2 = \frac{b_i}{4\pi^2} + \sigma_{\text{res}}^2$.

The integral of this kernel over all space equals the element's scattering power at zero angle:

$$\iiint K_{\text{element}}(\mathbf{r}) \, d^3\mathbf{r} = \sum_{i=1}^{5} a_i = f_e(0)$$

This value varies by element (e.g., ~2.808 for Carbon, ~3.85 for Nitrogen, ~4.76 for Oxygen) and represents the total scattering power of the atom.

#### Discretization and Normalization

Our implementation follows a simple two-step process to create a properly normalized discrete kernel for convolution:

**Step 1: Evaluate at grid points**

For each grid point $(i,j,k)$ within the kernel extent (typically 3-5$\sigma$), compute the kernel value directly using the continuous function:

$$K_{ijk} = \sum_{m=1}^{5} \frac{a_m}{(2\pi\sigma_m^2)^{3/2}} \exp\left(-\frac{(i h)^2 + (j h)^2 + (k h)^2}{2\sigma_m^2}\right)$$

where $h$ is the grid spacing in Ångströms.

**Step 2: Calculate sum and normalize**

Compute the total sum of all kernel values: $S = \sum_{i,j,k} K_{ijk}$

Then scale all kernel values so their sum equals the theoretical $f_e(0)$:

$$K_{ijk}^{\text{norm}} = K_{ijk} \times \frac{f_e(0)}{S}$$

This normalization ensures that when the kernel is used in convolution with a delta-map containing $N$ atoms of this element, the total density in the output map will be $N \times f_e(0)$, correctly representing the total scattering power.









### EMmer / GEMMI Method

The EMmer method, inspired by the GEMMI library [Wojdyr2022], uses the complete International Tables Vol. C coefficients (c4322.lib) to generate real-space density maps through direct summation of Gaussian functions. This approach provides highly accurate electron scattering factors by incorporating the full 5-Gaussian parameterization with temperature factor optimization.

#### Gaussian Coefficient Database

The method employs the comprehensive coefficient table from International Tables Vol. C [InternationalTables2006], containing parameters for elements up to Californium (Z=98). For each element, five Gaussian terms capture the scattering behavior:

$$f_e(s) = \sum_{i=1}^{5} a_i \exp(-b_i s^2)$$

where $s = \sin\theta/\lambda$ is the scattering vector. The coefficients $a_i$ and $b_i$ are empirically determined to fit experimental scattering factors across the full resolution range.

#### Real-Space Density Formulation

In real space, each Gaussian term in the scattering factor corresponds to a Gaussian density distribution. The total density contribution from an atom at position $\mathbf{r}_j$ is:

$$\rho_j(\mathbf{r}) = \sum_{i=1}^{5} \frac{a_i}{(2\pi\sigma_i^2)^{3/2}} \exp\left(-\frac{|\mathbf{r} - \mathbf{r}_j|^2}{2\sigma_i^2}\right)$$

The width of each Gaussian component combines the intrinsic atomic scattering width with thermal motion (B-factor) and optional resolution-dependent blur:

$$\sigma_i^2 = \frac{b_i}{4\pi^2} + \frac{B_j}{8\pi^2} + \sigma_{\text{blur}}^2$$

where:
- $b_i$ are the coefficients from International Tables
- $B_j$ is the atomic B-factor from the PDB file
- $\sigma_{\text{blur}}$ is an optional additional blurring term

#### Refmac-Compatible Blur

A key feature of the EMmer method is the optional application of Refmac-compatible blur [Murshudov1997, Murshudov2011]. When enabled, an effective B-factor is calculated to match the conventions used in crystallographic refinement:

$$B_{\text{eff}} = \frac{8\pi^2}{1.1} \left(\frac{d_{\min}}{2R}\right)^2 - B_{\min}$$

where:
- $d_{\min}$ is the high-resolution limit (target resolution)
- $R = 1.5$ is the Shannon rate (configurable via `--emmer-rate`)
- $B_{\min}$ is the minimum B-factor in the structure

This effective B-factor is then incorporated into the Gaussian widths, producing maps that are directly compatible with Refmac-sharpened maps used in crystallography and cryo-EM refinement workflows.

#### Direct Real-Space Summation

Unlike Fourier-based methods, EMmer generates maps through direct summation in real space, which offers several advantages:

1. **Natural handling of B-factors**: Individual atomic B-factors are directly incorporated into the Gaussian widths
2. **Accurate near-field behavior**: The multi-Gaussian representation correctly models the complex density near atomic nuclei
3. **No FFT artifacts**: Direct summation avoids Gibbs phenomena and aliasing issues
4. **Parallelizable**: The atom-centered contributions can be computed independently and summed

For each atom, the contribution is calculated only within a sphere where the density exceeds a cutoff threshold ($10^{-5}$ by default), determined by analyzing the decay of each Gaussian component:

$$r_{\text{cut}} = \max_i \sqrt{\frac{\ln(a_i / \text{cutoff})}{-w_i}}$$

where $w_i$ is the width parameter for the i-th Gaussian.

#### Symmetry Expansion

When working with crystallographic data, the method can optionally apply space group symmetry to expand the asymmetric unit atoms throughout the unit cell. This ensures that the final map includes all symmetry-related density contributions, essential for crystal structure visualization and validation.

#### Output Alignment

For compatibility with standard visualization tools, the EMmer method includes an optional output alignment step that applies the MRC/CCP4 convention for map orientation. This ensures that the generated maps display correctly in programs like ChimeraX, Coot, and PyMOL without requiring manual reorientation.




### Peng1996 / International Tables Method

This method uses the 5-Gaussian parameterization of electron scattering factors from the **International Tables for Crystallography** [InternationalTables2006] as parameterized by Peng and colleagues [Peng1996]. The scattering factor for an element as a function of resolution is:

$$f_e(s) = \sum_{i=1}^{5} a_i \exp(-b_i s^2)$$

where $s = \sin\theta/\lambda$ is the scattering vector. The real-space density is obtained by inverse Fourier transform:

$$\rho(r) = \mathcal{F}^{-1}[f_e(s)] = \sum_{i=1}^{5} \frac{a_i}{(2\pi\sigma_i^2)^{3/2}} \exp\left(-\frac{r^2}{2\sigma_i^2}\right)$$

The total width includes both the intrinsic atomic scattering and the resolution broadening:

$$\sigma_i^2 = \frac{b_i}{4\pi^2} + \sigma_{\text{res}}^2$$

where $\sigma_{\text{res}}$ is determined by the target resolution using the chosen criterion.

Our implementation uses two distinct tables from the Peng1996 paper:
- **Table 1** (pages 260-261): Elastic scattering factors for $s$ up to $2.0 \text{Å}^{-1}$. Used for lower resolution maps.
- **Table 3** (pages 264-265): Elastic scattering factors for $s$ up to $6.0 \text{Å}^{-1}$. Used for higher resolution maps, providing a more accurate representation at larger scattering angles.

**Amplitude Modes**:
- **Peng1996**: Uses the sum of coefficients $\sum a_i$ as the atomic scattering power at zero angle [Peng1996].
- **Atomic Number**: Scales by $Z$ (EMAN2-style), useful for mass-weighted maps where density is proportional to atomic mass [Tang2007].

### ChimeraX molmap Method

The ChimeraX method implemented here replicates the `molmap` command from UCSF ChimeraX [Goddard2018, Pettersen2020], which generates density maps by placing Gaussian functions at each atom position. The implementation is based on the actual ChimeraX C++ and Python code [ChimeraXSource], ensuring compatibility with maps produced by ChimeraX.

#### Mathematical Formulation

In the ChimeraX `molmap` algorithm, each atom contributes a normalized 3D Gaussian density:

$$\rho_i(\mathbf{r}) = \frac{Z_i}{(2\pi\sigma^2)^{3/2}} \exp\left(-\frac{|\mathbf{r} - \mathbf{r}_i|^2}{2\sigma^2}\right)$$

where:
- $Z_i$ is the atomic number (element number) used as the scattering power
- $\mathbf{r}_i$ is the atom position in Ångströms
- $\sigma$ is the standard deviation of the Gaussian, determined by the target resolution

The total density at any grid point is the sum of contributions from all atoms within a cutoff distance:

$$\rho(\mathbf{r}) = \sum_{i: |\mathbf{r} - \mathbf{r}_i| < n\sigma} \rho_i(\mathbf{r})$$

with the default cutoff $n = 5$ standard deviations, as implemented in the C++ function `sum_of_gaussians()` [ChimeraXSource].

#### Resolution to Sigma Conversion

The relationship between the target resolution $R$ and the Gaussian width $\sigma$ in ChimeraX is:

$$\sigma = \frac{R}{\pi\sqrt{2}} \approx 0.225R$$

This relationship is defined by the `sigma_factor` parameter, with the default value $1/(\pi\sqrt{2})$. The Fourier transform of the Gaussian falls to $1/e$ of its maximum at wavenumber $1/R$ [Goddard2018].

#### Grid Generation

The map grid is generated using the following steps:

1. **Bounding Box**: Calculate the minimum bounding box of all atom coordinates: $x_{\min} = \min_i(x_i)$, $x_{\max} = \max_i(x_i)$ (and similarly for $y$ and $z$ dimensions).

2. **Padding**: Add padding on all sides: $x_{\min} \leftarrow x_{\min} - p$, $x_{\max} \leftarrow x_{\max} + p$, where the default padding is $p = 3R$ (controlled by `edge_padding`).

3. **Grid Dimensions**: Calculate the number of grid points: $n_x = \lceil (x_{\max} - x_{\min})/s \rceil + 1$, where the grid spacing $s$ defaults to $R/3$ (controlled by `grid_spacing`).

4. **Cubic Option**: If cube=True, all dimensions are set to the maximum and made even: $n = \max(n_x, n_y, n_z)$, $n \leftarrow n + (n \bmod 2)$.

5. **Origin**: Set the grid origin to the minimum coordinates: $\mathbf{r}_{\text{origin}} = (x_{\min}, y_{\min}, z_{\min})$.

#### Core Algorithm Implementation

The core computation in `gaussian.cpp` [ChimeraXSource] implements an optimized summation. For each atom, the algorithm:

1. Calculates bounds in grid coordinates: $i_{\min} = \lceil i_c - n \cdot \sigma/s \rceil$, $i_{\max} = \lfloor i_c + n \cdot \sigma/s \rfloor$.

2. Computes contributions within the bounding box using nested loops ordered for cache efficiency $(k, j, i)$:

   $$\rho(i,j,k) \mathrel{+}= Z \cdot \exp\left(-\frac{1}{2}\left[\left(\frac{i-i_c}{\sigma/s}\right)^2 + \left(\frac{j-j_c}{\sigma/s}\right)^2 + \left(\frac{k-k_c}{\sigma/s}\right)^2\right]\right)$$

#### Normalization

After summation, the map is normalized by:

$$\rho_{\text{norm}}(\mathbf{r}) = \rho(\mathbf{r}) \cdot (2\pi)^{-3/2} \sigma^{-3}$$

This normalization ensures that the integral of each Gaussian equals its atomic number $Z$. The final map is then scaled to a maximum value of 1.0 for visualization compatibility [Pettersen2020].

#### Balls Mode

An alternative representation using "balls" with Gaussian falloff is available when balls=True. In this mode, each atom contributes a constant value of 1 within its van der Waals radius $r_{\text{vdW}}$, with a Gaussian falloff outside:

$$\rho_i(\mathbf{r}) = \begin{cases} 
1 & |\mathbf{r} - \mathbf{r}_i| \leq r_{\text{vdW}} \\
\exp\left(-\frac{1}{2}\left(\frac{|\mathbf{r} - \mathbf{r}_i| - r_{\text{vdW}}}{\sigma}\right)^2\right) & |\mathbf{r} - \mathbf{r}_i| > r_{\text{vdW}}
\end{cases}$$

This mode produces maps where isosurfaces approximate the van der Waals envelope when contoured at low levels [Goddard2018].

### Situs Method

The Situs method implemented here follows the real-space convolution approach established in the Situs package [Wriggers2010, Wriggers1999]. It generates a density map by first projecting atomic structures onto a grid and then smoothing the result with a kernel function. This two-step process is designed to produce maps that correspond to a user-specified resolution.

#### Kernel Functions

The method offers five distinct kernel types, each with a different mathematical form. The choice of kernel affects the shape of the resulting density. The kernels are defined as functions of the distance $r$ from the kernel center. Two key parameters define the kernel's width:

- **Half-max radius** ($r_h$): The distance at which the kernel's value drops to half of its maximum.
- **Situs resolution** ($r_s$): An empirical resolution measure defined as $r_s = 2\sigma$, where $\sigma$ is the standard deviation of a 3D Gaussian kernel. The relationship between $r_s$ and the half-max radius is kernel-dependent.

The user can specify the target resolution in two modes:
1. **Half-max radius mode** (positive input value): The input value is used directly as $r_h$.
2. **$2\sigma$ mode** (negative input value): The absolute input value is used as $r_s$, and $r_h$ is calculated accordingly.

The kernel functions are:

| Kernel Type | Function $K(r)$ | Half-max relation |
|-------------|---------------|-------------------|
| **Gaussian** | $\exp\left(-\dfrac{3r^2}{2\sigma^2}\right)$ | $r_h = \sigma\sqrt{\dfrac{\ln 2}{1.5}}$ |
| **Triangular** | $\max\left(0, 1 - \dfrac{r}{2r_h}\right)$ | $r_h$ is the half-max radius |
| **Semi-Epanechnikov** | $\max\left(0, 1 - \dfrac{r^{1.5}}{2r_h^{1.5}}\right)$ | $r_h$ is the half-max radius |
| **Epanechnikov** | $\max\left(0, 1 - \dfrac{r^2}{2r_h^2}\right)$ | $r_h$ is the half-max radius |
| **Hard Sphere** | $\max\left(0, 1 - \dfrac{r^{60}}{2r_h^{60}}\right)$ | $r_h$ is the half-max radius |

The **Epanechnikov kernel** is a special case, as it is known to be optimal for minimizing the asymptotic mean integrated square error in kernel density estimation [Epanechnikov1969].

#### Map Generation Workflow

The map generation follows a two-step process, consistent with the description of the `pdb2vol` tool from the Situs documentation [SitusDoc]:

1. **Projection to Lattice**: Atoms are projected onto a cubic lattice using **trilinear interpolation**. Each atom, with a given position and weight (atomic mass or unity), contributes to the eight surrounding voxels. This creates an intermediate "lattice" representation of the structure.

2. **Kernel Convolution**: The lattice is then convolved with the selected 3D kernel. This step smooths the structure to the desired resolution and produces the final density map. The kernel's width is determined by the user-specified resolution and kernel type, as described above.

#### Lattice Variance Correction

The initial projection onto a lattice introduces an inherent, small amount of blurring. An optional correction can be applied that accounts for this lattice smoothing [Wriggers2012]. This is achieved by adjusting the kernel's variance:

$$\sigma_{\text{corrected}}^2 = \sigma_{\text{target}}^2 - \sigma_{\text{lattice}}^2$$

where $\sigma_{\text{target}}$ is the width required to achieve the desired resolution, and $\sigma_{\text{lattice}}$ is the standard deviation of the point-spread function introduced by the trilinear projection. This correction ensures that the final map more accurately matches the target resolution [Wriggers2010, Wriggers1999].

### EMmer / GEMMI Method

This method, inspired by the GEMMI library [Wojdyr2022] and EMmer, uses the complete International Tables Vol. C coefficients (c4322.lib) with Refmac-compatible blur [Murshudov1997]:

The effective B-factor including resolution-dependent blur:

$$B_{\text{eff}} = \frac{8\pi^2}{1.1} \left(\frac{d_{\min}}{2R}\right)^2 - B_{\min}$$

where:
- $d_{\min}$ is the high-resolution limit (target resolution)
- $R = 1.5$ is the Shannon rate (default)
- $B_{\min}$ is the minimum B-factor in the structure

The atomic scattering factors are then:

$$f_e(s) = \sum_{i=1}^{5} a_i \exp\left(-(b_i + B_{\text{eff}}) s^2\right)$$

This produces maps compatible with Refmac sharpening/blurring conventions [Murshudov1997, Murshudov2011], making them suitable for refinement and validation in crystallographic and cryo-EM workflows. The coefficients used in this method are taken directly from the International Tables Vol. C [InternationalTables2006] as compiled in the GEMMI library's c4322.lib [Wojdyr2022].

## Input Parameters

### Required Parameters

| Parameter | Format | Description |
|-----------|--------|-------------|
| `-i FILE` | string | Input PDB file path |
| `-o FILE` | string | Output MRC file path |

### Method Selection

| Parameter | Format | Description |
|-----------|--------|-------------|
| `--method STR` | string | Generation method: `peng1996` (default), `chimerax`, `situs`, `emmer` |

### Common Parameters

| Parameter | Format | Default | Description |
|-----------|--------|---------|-------------|
| `-r FLOAT` | float | 6.0 | Target resolution in Ångströms |
| `-c STR` | string | `rayleigh` | Resolution criterion: `rayleigh`, `chimerax`, `eman2`, `fsc0143`, `fsc05` |
| `-s FLOAT` | float | auto | Voxel size in Å (default = resolution/3) |
| `-p FLOAT` | float | 3.0 | Padding around atoms in Å |
| `-t INT` | int | 0 | Number of OpenMP threads (0 = auto) |
| `-v` | flag | - | Verbose output |
| `-q` | flag | - | Quiet mode |
| `-h` | flag | - | Show help message |

### PDB Filtering Parameters

| Parameter | Format | Default | Description |
|-----------|--------|---------|-------------|
| `--filter-h` | flag | on | Filter out hydrogen atoms |
| `--no-filter-h` | flag | - | Keep hydrogen atoms |
| `--filter-w` | flag | off | Filter out water molecules |
| `-b FLOAT` | float | 0.0 | B-factor cutoff (exclude atoms with B > value, 0 = no cutoff) |

### Peng1996/AtomicNumber Parameters

| Parameter | Format | Default | Description |
|-----------|--------|---------|-------------|
| `-a STR` | string | `peng1996` | Amplitude mode: `peng1996` or `atomic-number` |
| `--no-bfac` | flag | - | Ignore B-factors from PDB file |

### ChimeraX-Specific Parameters

| Parameter | Format | Default | Description |
|-----------|--------|---------|-------------|
| `--cutoff FLOAT` | float | 5.0 | Cutoff range in sigma ($n\sigma$) |
| `--no-norm` | flag | - | Skip final normalization to maximum 1.0 |

The ChimeraX method uses the following relationships:

| Parameter | Formula | Description |
|-----------|---------|-------------|
| Sigma | $\sigma = R/(\pi\sqrt{2})$ | Gaussian width |
| Cutoff radius | $r_{\text{cut}} = n \cdot \sigma$ | Maximum distance for atom contributions |
| Grid spacing | $s = R/3$ | Default voxel size |

### Situs-Specific Parameters

| Parameter | Format | Default | Description |
|-----------|--------|---------|-------------|
| `--situs-kernel NUM` | int (1-5) | 1 | Kernel type by number |
| `--situs-kernel-type TYPE` | string | `gaussian` | Kernel type by name |
| `--situs-halfmax FLOAT` | float | - | Set resolution as half-max radius (positive value) |
| `--situs-2sigma FLOAT` | float | - | Set resolution as $2\sigma$ (negative value) |
| `--situs-margin INT` | int | 2 | Margin voxels around structure |
| `--situs-mass` | flag | on | Use atomic mass weighting |
| `--situs-no-mass` | flag | - | Use unit weights (each atom = 1) |
| `--situs-correction` | flag | on | Apply lattice smoothing correction |
| `--situs-no-correction` | flag | - | Skip lattice correction |
| `--situs-amplitude FLOAT` | float | 1.0 | Kernel amplitude scaling factor |

#### Situs Kernel Types

| Number | Name | Function |
|--------|------|----------|
| 1 | `gaussian` | $\exp(-1.5 r^2/\sigma^2)$ |
| 2 | `triangular` | $\max(0, 1 - r/(2r_h))$ |
| 3 | `semi-epanechnikov` | $\max(0, 1 - r^{1.5}/(2r_h^{1.5}))$ |
| 4 | `epanechnikov` | $\max(0, 1 - r^2/(2r_h^2))$ |
| 5 | `hard-sphere` | $\max(0, 1 - r^{60}/(2r_h^{60}))$ |

#### Situs Resolution Modes

| Mode | Input | Interpretation |
|------|-------|----------------|
| Half-max radius | positive value | $r_h =$ input value |
| $2\sigma$ mode | negative value | $r_s = \|\text{input}\|$, $\sigma = r_s/2$ |

### EMmer-Specific Parameters

| Parameter | Format | Default | Description |
|-----------|--------|---------|-------------|
| `--emmer-align` | flag | on | Apply MRC output alignment (flip+rotate) |
| `--emmer-no-align` | flag | - | Skip output alignment |
| `--emmer-refmac-blur` | flag | on | Apply Refmac-compatible blur |
| `--emmer-no-blur` | flag | - | Skip Refmac blur |
| `--emmer-blur FLOAT` | float | 0.0 | Manual blur value in Å² (0 = auto) |
| `--emmer-symmetry` | flag | on | Apply space group symmetry |
| `--emmer-no-symmetry` | flag | - | Skip symmetry expansion |
| `--emmer-cutoff FLOAT` | float | 1e-5 | Density cutoff for radius determination |
| `--emmer-rate FLOAT` | float | 1.5 | Shannon rate for grid spacing |

#### EMmer Blur Calculation

When `--emmer-refmac-blur` is enabled and no manual blur is provided, the blur is calculated as:

$$B_{\text{eff}} = \frac{8\pi^2}{1.1} \left(\frac{d_{\min}}{2R}\right)^2 - B_{\min}$$

where:
- $d_{\min}$ is the target resolution
- $R = 1.5$ is the Shannon rate
- $B_{\min}$ is the minimum B-factor in the structure

## Examples

Default Peng1996 mode:

<tt>pdb2mrc -i 1ake.pdb -o 1ake.mrc -r 8.0</tt>

ChimeraX mode with custom cutoff:

<tt>pdb2mrc -i 1ake.pdb -o 1ake_chx.mrc -r 8.0 --method chimerax --cutoff 3.0</tt>

Situs with Epanechnikov kernel:

<tt>pdb2mrc -i 1ake.pdb -o 1ake_situs.mrc --method situs --situs-kernel-type epanechnikov --situs-halfmax 8.0</tt>

EMmer with manual blur:

<tt>pdb2mrc -i 1ake.pdb -o 1ake_emmer.mrc --method emmer -r 8.0 --emmer-blur 25.0 --emmer-no-align</tt>

Full filtering options:

<tt>pdb2mrc -i 1ake.pdb -o 1ake_filtered.mrc -r 6.0 --filter-h --filter-w -b 50.0</tt>

## Building the Project

### Dependencies

- **Intel oneAPI Base Toolkit & Intel oneAPI HPC Toolkit**: Provides Intel IPP (Integrated Performance Primitives) and MKL (Math Kernel Library) for FFT operations and optimized vector functions.
- **CMake** (version 3.12 or higher): For cross-platform build configuration.
- **OpenMP**: For shared-memory parallelization (usually included with modern compilers).
- **C++17 compliant compiler**: e.g., GCC 7+, Clang 5+, MSVC 2017+.

### Build Instructions

<tt>mkdir build</tt>
<tt>cd build</tt>
<tt>cmake ..</tt>
<tt>make</tt>

I'll update the References section with DOIs and links:

## References

[ChimeraXSource] UCSF ChimeraX. (2023). *Gaussian summation implementation*. GitHub repository. https://github.com/ucsf-chimerax/chimerax

[Epanechnikov1969] Epanechnikov, V. A. (1969). Non-parametric estimation of a multivariate probability density. *Theory of Probability & Its Applications*, 14(1), 153-158. https://doi.org/10.1137/1114019

[Goddard2018] Goddard, T. D., Huang, C. C., Meng, E. C., Pettersen, E. F., Couch, G. S., Morris, J. H., & Ferrin, T. E. (2018). UCSF ChimeraX: Meeting modern challenges in visualization and analysis. *Protein Science*, 27(1), 14-25. https://doi.org/10.1002/pro.3235

[InternationalTables2006] Prince, E. (Ed.). (2006). *International Tables for Crystallography, Vol. C: Mathematical, physical and chemical tables* (3rd ed.). Springer. https://doi.org/10.1107/97809553602060000103

[Murshudov1997] Murshudov, G. N., Vagin, A. A., & Dodson, E. J. (1997). Refinement of macromolecular structures by the maximum-likelihood method. *Acta Crystallographica Section D*, 53(3), 240-255. https://doi.org/10.1107/S0907444996012255

[Murshudov2011] Murshudov, G. N., Skubák, P., Lebedev, A. A., Pannu, N. S., Steiner, R. A., Nicholls, R. A., Winn, M. D., Long, F., & Vagin, A. A. (2011). REFMAC5 for the refinement of macromolecular crystal structures. *Acta Crystallographica Section D*, 67(4), 355-367. https://doi.org/10.1107/S0907444911001314

[Peng1996] Peng, L. M., Ren, G., Dudarev, S. L., & Whelan, M. J. (1996). Robust parameterization of elastic and absorptive electron atomic scattering factors. *Acta Crystallographica Section A*, 52(2), 257-276. https://doi.org/10.1107/S0108767395014371

[Pettersen2020] Pettersen, E. F., Goddard, T. D., Huang, C. C., Meng, E. C., Couch, G. S., Croll, T. I., Morris, J. H., & Ferrin, T. E. (2020). UCSF ChimeraX: Structure visualization for researchers, educators, and developers. *Protein Science*, 30(1), 70-82. https://doi.org/10.1002/pro.3943

[Rayleigh1879] Rayleigh, L. (1879). Investigations in optics, with special reference to the spectroscope. *Philosophical Magazine*, 8(49), 261-274. https://doi.org/10.1080/14786447908639684

[SitusDoc] Situs. (2023). *pdb2vol - Create a volumetric map from a PDB*. Online documentation. https://situs.biomachina.org

[Tang2007] Tang, G., Peng, L., Baldwin, P. R., Mann, D. S., Jiang, W., Rees, I., & Ludtke, S. J. (2007). EMAN2: An extensible image processing suite for electron microscopy. *Journal of Structural Biology*, 157(1), 38-46. https://doi.org/10.1016/j.jsb.2006.05.009

[Wojdyr2022] Wojdyr, M. (2022). GEMMI: A library for structural biology. *Journal of Open Source Software*, 7(73), 4200. https://doi.org/10.21105/joss.04200

[Wriggers1999] Wriggers, W., Milligan, R. A., & McCammon, J. A. (1999). Situs: A package for docking crystal structures into low-resolution maps from electron microscopy. *Journal of Structural Biology*, 125(2-3), 185-195. https://doi.org/10.1006/jsbi.1998.4080

[Wriggers2010] Wriggers, W. (2010). Using Situs for the integration of multi-resolution structures. *Biophysical Reviews*, 2(1), 21-27. https://doi.org/10.1007/s12551-009-0024-5

[Wriggers2012] Wriggers, W. (2012). Conventions and workflows for using Situs. *Acta Crystallographica Section D*, 68(4), 344-351. https://doi.org/10.1107/S0907444911049791

[Cheng2015] Cheng, A., Henderson, R., Mastronarde, D., Ludtke, S. J., Schoenmakers, R. H., Short, J., ... & Agard, D. A. (2015). MRC2014: Extensions to the MRC format header for electron cryo-microscopy and tomography. *Journal of Structural Biology*, 192(2), 146-150. https://doi.org/10.1016/j.jsb.2015.04.002


## Citation

If you use pdb2mrc in your research, please cite:

[Your citation information here]

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

We thank the developers of UCSF ChimeraX, Situs, GEMMI, and other structural biology software for their foundational work that inspired this implementation.
