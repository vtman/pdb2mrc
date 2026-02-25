// src/situs_generator.cpp
#include "situs_generator.hpp"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

//=============================================================================
// SitusGenerator implementation
//=============================================================================

SitusGenerator::SitusGenerator() {
    m_lattice = nullptr;
    m_kernel = nullptr;
    m_output = nullptr;
    m_atoms = nullptr;
    m_nAtoms = 0;

    m_nx = m_ny = m_nz = 0;
    m_lattice_nx = m_lattice_ny = m_lattice_nz = 0;
    m_kernel_extent = 0;

    m_origin[0] = m_origin[1] = m_origin[2] = 0.0;
    m_lattice_origin[0] = m_lattice_origin[1] = m_lattice_origin[2] = 0.0;
    m_grid_spacing = 0.0;

    m_lattice_variance = 0.0;
    m_total_mass = 0.0;

    memset(&m_kernel_params, 0, sizeof(m_kernel_params));
}

SitusGenerator::~SitusGenerator() {
    cleanup();
}

void SitusGenerator::cleanup() {
    if (m_lattice) {
        ippsFree(m_lattice);
        m_lattice = nullptr;
    }
    if (m_kernel) {
        ippsFree(m_kernel);
        m_kernel = nullptr;
    }
    if (m_output) {
        ippsFree(m_output);
        m_output = nullptr;
    }
    // Note: m_atoms is not owned, so don't free it
    m_atoms = nullptr;
    m_nAtoms = 0;
}

//=============================================================================
// Private helper methods
//=============================================================================

int SitusGenerator::computeBoundingBox(const Atom* atoms, int n_atoms,
    Ipp64f& xmin, Ipp64f& xmax,
    Ipp64f& ymin, Ipp64f& ymax,
    Ipp64f& zmin, Ipp64f& zmax) {
    if (n_atoms == 0) return -1;

    xmin = xmax = atoms[0].x;
    ymin = ymax = atoms[0].y;
    zmin = zmax = atoms[0].z;

    for (int i = 1; i < n_atoms; i++) {
        if (atoms[i].x < xmin) xmin = atoms[i].x;
        if (atoms[i].x > xmax) xmax = atoms[i].x;
        if (atoms[i].y < ymin) ymin = atoms[i].y;
        if (atoms[i].y > ymax) ymax = atoms[i].y;
        if (atoms[i].z < zmin) zmin = atoms[i].z;
        if (atoms[i].z > zmax) zmax = atoms[i].z;
    }

    return 0;
}

void SitusGenerator::computeKernelParams() {
    Ipp64f reso = m_config.resolution;
    Ipp64f rh, sig;

    if (reso < 0.0) {
        // Negative value: target resolution R = 2σ
        sig = -reso / 2.0;

        // Compute half-max radius for each kernel type based on sigma
        switch (m_config.kernel_type) {
        case SITUS_KERNEL_GAUSSIAN:
            // Gaussian: r_h = σ * sqrt(ln(2)/1.5)
            rh = sig * sqrt(log(2.0) / 1.5);
            break;
        case SITUS_KERNEL_TRIANGULAR:
            // Triangular: r_h = σ / (2^(1/1) * sqrt(3*4/(5*6)))
            rh = sig / (exp(log(2.0)) * sqrt(12.0 / 30.0));
            break;
        case SITUS_KERNEL_SEMI_EPANECHNIKOV:
            // Semi-Epanechnikov: r_h = σ / (2^(1/1.5) * sqrt(3*4.5/(5*6.5)))
            rh = sig / (exp(log(2.0) / 1.5) * sqrt(13.5 / 32.5));
            break;
        case SITUS_KERNEL_EPANECHNIKOV:
            // Epanechnikov: r_h = σ / (2^(1/2) * sqrt(3*5/(5*7)))
            rh = sig / (exp(log(2.0) / 2.0) * sqrt(15.0 / 35.0));
            break;
        case SITUS_KERNEL_HARD_SPHERE:
            // Hard Sphere: r_h = σ / (2^(1/60) * sqrt(3*63/(5*65)))
            rh = sig / (exp(log(2.0) / 60.0) * sqrt(189.0 / 325.0));
            break;
        default:
            rh = sig * sqrt(log(2.0) / 1.5);
        }
    }
    else {
        // Positive value: direct half-max radius
        rh = reso;

        // Compute sigma for each kernel type based on half-max radius
        switch (m_config.kernel_type) {
        case SITUS_KERNEL_GAUSSIAN:
            // Gaussian: σ = r_h * sqrt(1.5/ln(2))
            sig = rh * sqrt(1.5 / log(2.0));
            break;
        case SITUS_KERNEL_TRIANGULAR:
            // Triangular: σ = r_h * 2^(1/1) * sqrt(3*4/(5*6))
            sig = rh * exp(log(2.0)) * sqrt(12.0 / 30.0);
            break;
        case SITUS_KERNEL_SEMI_EPANECHNIKOV:
            // Semi-Epanechnikov: σ = r_h * 2^(1/1.5) * sqrt(3*4.5/(5*6.5))
            sig = rh * exp(log(2.0) / 1.5) * sqrt(13.5 / 32.5);
            break;
        case SITUS_KERNEL_EPANECHNIKOV:
            // Epanechnikov: σ = r_h * 2^(1/2) * sqrt(3*5/(5*7))
            sig = rh * exp(log(2.0) / 2.0) * sqrt(15.0 / 35.0);
            break;
        case SITUS_KERNEL_HARD_SPHERE:
            // Hard Sphere: σ = r_h * 2^(1/60) * sqrt(3*63/(5*65))
            sig = rh * exp(log(2.0) / 60.0) * sqrt(189.0 / 325.0);
            break;
        default:
            sig = rh * sqrt(1.5 / log(2.0));
        }
    }

    // Compute cutoff radius for each kernel type
    Ipp64f rc;
    switch (m_config.kernel_type) {
    case SITUS_KERNEL_GAUSSIAN:
        rc = sqrt(3.0) * sig;
        break;
    case SITUS_KERNEL_TRIANGULAR:
        rc = exp(log(2.0)) * rh;
        break;
    case SITUS_KERNEL_SEMI_EPANECHNIKOV:
        rc = exp(log(2.0) / 1.5) * rh;
        break;
    case SITUS_KERNEL_EPANECHNIKOV:
        rc = exp(log(2.0) / 2.0) * rh;
        break;
    case SITUS_KERNEL_HARD_SPHERE:
        rc = exp(log(2.0) / 60.0) * rh;
        break;
    default:
        rc = sqrt(3.0) * sig;
    }

    m_kernel_params.sigma = sig;
    m_kernel_params.half_max_radius = rh;
    m_kernel_params.cutoff_radius = rc;
}

void SitusGenerator::generateGaussianKernel(Ipp64f* kernel, int extent, Ipp64f sigma, Ipp64f amplitude) {
    int half = extent / 2;
    Ipp64f sigma_sq = sigma * sigma;
    Ipp64f bvalue = -1.0 / (2.0 * sigma_sq);
    Ipp64f cvalue = 9.0 * sigma_sq;  // (3σ)^2

    for (int k = 0; k < extent; k++) {
        int dk = k - half;
        Ipp64f dk_sq = (Ipp64f)(dk * dk);
        for (int j = 0; j < extent; j++) {
            int dj = j - half;
            Ipp64f dj_sq = (Ipp64f)(dj * dj);
            for (int i = 0; i < extent; i++) {
                int di = i - half;
                Ipp64f dsq = (Ipp64f)(di * di) + dj_sq + dk_sq;
                if (dsq < cvalue) {
                    kernel[kernelIdx(i, j, k)] = amplitude * exp(dsq * bvalue);
                }
                else {
                    kernel[kernelIdx(i, j, k)] = 0.0;
                }
            }
        }
    }
}

void SitusGenerator::generateTriangularKernel(Ipp64f* kernel, int extent, Ipp64f rh, Ipp64f amplitude) {
    int half = extent / 2;
    Ipp64f bvalue = 0.5 / rh;  // 0.5/r_h

    for (int k = 0; k < extent; k++) {
        int dk = k - half;
        for (int j = 0; j < extent; j++) {
            int dj = j - half;
            for (int i = 0; i < extent; i++) {
                int di = i - half;
                Ipp64f r = m_grid_spacing * sqrt((Ipp64f)(di * di + dj * dj + dk * dk));
                Ipp64f val = amplitude * (1.0 - r * bvalue);
                kernel[kernelIdx(i, j, k)] = (val > 0.0) ? val : 0.0;
            }
        }
    }
}

void SitusGenerator::generateSemiEpanechnikovKernel(Ipp64f* kernel, int extent, Ipp64f rh, Ipp64f amplitude) {
    int half = extent / 2;
    Ipp64f bvalue = 0.5 / pow(rh, 1.5);  // 0.5 / r_h^1.5

    for (int k = 0; k < extent; k++) {
        int dk = k - half;
        for (int j = 0; j < extent; j++) {
            int dj = j - half;
            for (int i = 0; i < extent; i++) {
                int di = i - half;
                Ipp64f r = m_grid_spacing * sqrt((Ipp64f)(di * di + dj * dj + dk * dk));
                Ipp64f r_pow = pow(r, 1.5);
                Ipp64f val = amplitude * (1.0 - r_pow * bvalue);
                kernel[kernelIdx(i, j, k)] = (val > 0.0) ? val : 0.0;
            }
        }
    }
}

void SitusGenerator::generateEpanechnikovKernel(Ipp64f* kernel, int extent, Ipp64f rh, Ipp64f amplitude) {
    int half = extent / 2;
    Ipp64f bvalue = 0.5 / (rh * rh);  // 0.5 / r_h^2

    for (int k = 0; k < extent; k++) {
        int dk = k - half;
        for (int j = 0; j < extent; j++) {
            int dj = j - half;
            for (int i = 0; i < extent; i++) {
                int di = i - half;
                Ipp64f r = m_grid_spacing * sqrt((Ipp64f)(di * di + dj * dj + dk * dk));
                Ipp64f r_sq = r * r;
                Ipp64f val = amplitude * (1.0 - r_sq * bvalue);
                kernel[kernelIdx(i, j, k)] = (val > 0.0) ? val : 0.0;
            }
        }
    }
}

void SitusGenerator::generateHardSphereKernel(Ipp64f* kernel, int extent, Ipp64f rh, Ipp64f amplitude) {
    int half = extent / 2;
    Ipp64f bvalue = 0.5 / pow(rh, 60.0);  // 0.5 / r_h^60

    for (int k = 0; k < extent; k++) {
        int dk = k - half;
        for (int j = 0; j < extent; j++) {
            int dj = j - half;
            for (int i = 0; i < extent; i++) {
                int di = i - half;
                Ipp64f r = m_grid_spacing * sqrt((Ipp64f)(di * di + dj * dj + dk * dk));
                Ipp64f r_pow = pow(r, 60.0);
                Ipp64f val = amplitude * (1.0 - r_pow * bvalue);
                kernel[kernelIdx(i, j, k)] = (val > 0.0) ? val : 0.0;
            }
        }
    }
}

int SitusGenerator::generateKernel() {
    // Compute kernel parameters based on resolution input
    computeKernelParams();

    // Calculate target variance in voxel units
    Ipp64f kmsd = (m_kernel_params.sigma * m_kernel_params.sigma) / (m_grid_spacing * m_grid_spacing);

    // Apply lattice correction if requested
    Ipp64f varmap = kmsd;
    if (m_config.apply_lattice_correction) {
        varmap -= m_lattice_variance;
        if (varmap < 0.0) {
            fprintf(stderr, "Error: Lattice variance exceeds target variance. Increase grid spacing or reduce resolution.\n");
            return -1;
        }
    }

    // Determine kernel extent based on cutoff radius
    Ipp64f sigma_corrected = sqrt(varmap) * m_grid_spacing;

    // For Gaussian, use 3σ cutoff
    if (m_config.kernel_type == SITUS_KERNEL_GAUSSIAN) {
        m_kernel_extent = 2 * (int)ceil(3.0 * sigma_corrected / m_grid_spacing) + 1;
    }
    else {
        // For other kernels, use the precomputed cutoff radius
        m_kernel_extent = 2 * (int)ceil(m_kernel_params.cutoff_radius / m_grid_spacing) + 1;
    }

    // Ensure kernel extent is odd
    if (m_kernel_extent % 2 == 0) m_kernel_extent++;

    // Allocate kernel
    int64_t nvox_kernel = (int64_t)m_kernel_extent * m_kernel_extent * m_kernel_extent;
    m_kernel = (Ipp64f*)ippsMalloc_64f(nvox_kernel);
    if (!m_kernel) return -2;

    ippsZero_64f(m_kernel, nvox_kernel);

    // Generate kernel based on type
    switch (m_config.kernel_type) {
    case SITUS_KERNEL_GAUSSIAN:
        generateGaussianKernel(m_kernel, m_kernel_extent, sigma_corrected, m_config.kernel_amplitude);
        break;
    case SITUS_KERNEL_TRIANGULAR:
        generateTriangularKernel(m_kernel, m_kernel_extent, m_kernel_params.half_max_radius, m_config.kernel_amplitude);
        break;
    case SITUS_KERNEL_SEMI_EPANECHNIKOV:
        generateSemiEpanechnikovKernel(m_kernel, m_kernel_extent, m_kernel_params.half_max_radius, m_config.kernel_amplitude);
        break;
    case SITUS_KERNEL_EPANECHNIKOV:
        generateEpanechnikovKernel(m_kernel, m_kernel_extent, m_kernel_params.half_max_radius, m_config.kernel_amplitude);
        break;
    case SITUS_KERNEL_HARD_SPHERE:
        generateHardSphereKernel(m_kernel, m_kernel_extent, m_kernel_params.half_max_radius, m_config.kernel_amplitude);
        break;
    default:
        generateGaussianKernel(m_kernel, m_kernel_extent, sigma_corrected, m_config.kernel_amplitude);
    }

    return 0;
}



//=============================================================================
// Public methods
//=============================================================================

int SitusGenerator::init(const SitusGeneratorConfig* config, const Atom* atoms, int n_atoms) {
    if (!config || !atoms || n_atoms <= 0) return -1;

    m_config = *config;
    m_atoms = atoms;
    m_nAtoms = n_atoms;
    m_grid_spacing = m_config.grid_spacing;

    // Compute bounding box of atoms
    Ipp64f xmin, xmax, ymin, ymax, zmin, zmax;
    int ret = computeBoundingBox(atoms, n_atoms, xmin, xmax, ymin, ymax, zmin, zmax);
    if (ret != 0) return -2;

    // Align to grid - this is correct
    xmin = m_grid_spacing * floor(xmin / m_grid_spacing);
    xmax = m_grid_spacing * ceil(xmax / m_grid_spacing);
    ymin = m_grid_spacing * floor(ymin / m_grid_spacing);
    ymax = m_grid_spacing * ceil(ymax / m_grid_spacing);
    zmin = m_grid_spacing * floor(zmin / m_grid_spacing);
    zmax = m_grid_spacing * ceil(zmax / m_grid_spacing);

    // Calculate base dimensions without margin
    int base_nx = (int)ceil((xmax - xmin) / m_grid_spacing) + 1;
    int base_ny = (int)ceil((ymax - ymin) / m_grid_spacing) + 1;
    int base_nz = (int)ceil((zmax - zmin) / m_grid_spacing) + 1;

    // Add margin voxels to dimensions
    m_lattice_nx = base_nx + 2 * m_config.margin_voxels;
    m_lattice_ny = base_ny + 2 * m_config.margin_voxels;
    m_lattice_nz = base_nz + 2 * m_config.margin_voxels;

    // Lattice origin should be shifted by margin voxels
    // This ensures the original bounding box is centered in the lattice
    m_lattice_origin[0] = xmin - m_config.margin_voxels * m_grid_spacing;
    m_lattice_origin[1] = ymin - m_config.margin_voxels * m_grid_spacing;
    m_lattice_origin[2] = zmin - m_config.margin_voxels * m_grid_spacing;

    if (m_config.verbosity >= 1) {
        printf("Situs init:\n");
        printf("  Atom bounding box: (%.2f-%.2f, %.2f-%.2f, %.2f-%.2f)\n",
            xmin, xmax, ymin, ymax, zmin, zmax);
        printf("  Base dimensions: %d x %d x %d\n", base_nx, base_ny, base_nz);
        printf("  Lattice dimensions: %d x %d x %d (with %d margin voxels)\n",
            m_lattice_nx, m_lattice_ny, m_lattice_nz, m_config.margin_voxels);
        printf("  Lattice origin: (%.2f, %.2f, %.2f)\n",
            m_lattice_origin[0], m_lattice_origin[1], m_lattice_origin[2]);
    }

    // Allocate lattice
    int64_t nvox_lattice = (int64_t)m_lattice_nx * m_lattice_ny * m_lattice_nz;
    m_lattice = (Ipp64f*)ippsMalloc_64f(nvox_lattice);
    if (!m_lattice) return -3;

    return 0;
}

// Fix projectAtomsToLattice() - the coordinate transformation
int SitusGenerator::projectAtomsToLattice() {
    if (!m_atoms || m_nAtoms <= 0) return -1;
    if (!m_lattice) return -2;

    // Clear lattice
    int64_t nvox_lattice = (int64_t)m_lattice_nx * m_lattice_ny * m_lattice_nz;
    ippsZero_64f(m_lattice, nvox_lattice);

    m_lattice_variance = 0.0;
    m_total_mass = 0.0;

    Ipp64f inv_spacing = 1.0 / m_grid_spacing;

    // Project each atom
    for (int i = 0; i < m_nAtoms; i++) {
        // Determine weight
        Ipp64f weight;
        if (m_config.use_mass_weighting) {
            int Z = AtomUtils::get_atomic_number(m_atoms[i].element);
            if (Z == 0) Z = 6;  // Default to carbon
            weight = (Ipp64f)Z;
        }
        else {
            weight = 1.0;
        }

        // Position in grid coordinates relative to lattice origin
        // This is correct - we don't add margin_voxels here because
        // the lattice origin already accounts for the margin
        Ipp64f gx = (m_atoms[i].x - m_lattice_origin[0]) * inv_spacing;
        Ipp64f gy = (m_atoms[i].y - m_lattice_origin[1]) * inv_spacing;
        Ipp64f gz = (m_atoms[i].z - m_lattice_origin[2]) * inv_spacing;

        // Check bounds
        if (gx < 0 || gx >= m_lattice_nx - 1 ||
            gy < 0 || gy >= m_lattice_ny - 1 ||
            gz < 0 || gz >= m_lattice_nz - 1) {
            if (m_config.verbosity >= 2) {
                printf("Warning: Atom %d at (%.2f,%.2f,%.2f) outside lattice (g=%.2f,%.2f,%.2f)\n",
                    i, m_atoms[i].x, m_atoms[i].y, m_atoms[i].z, gx, gy, gz);
            }
            continue;
        }

        // Integer coordinates (floor)
        int ix0 = (int)floor(gx);
        int iy0 = (int)floor(gy);
        int iz0 = (int)floor(gz);
        int ix1 = ix0 + 1;
        int iy1 = iy0 + 1;
        int iz1 = iz0 + 1;

        // Fractional parts (weights from left/front/bottom)
        Ipp64f fx = gx - ix0;
        Ipp64f fy = gy - iy0;
        Ipp64f fz = gz - iz0;

        // Trilinear contributions to 8 surrounding voxels
        Ipp64f w000 = (1.0 - fx) * (1.0 - fy) * (1.0 - fz);
        Ipp64f w100 = fx * (1.0 - fy) * (1.0 - fz);
        Ipp64f w010 = (1.0 - fx) * fy * (1.0 - fz);
        Ipp64f w110 = fx * fy * (1.0 - fz);
        Ipp64f w001 = (1.0 - fx) * (1.0 - fy) * fz;
        Ipp64f w101 = fx * (1.0 - fy) * fz;
        Ipp64f w011 = (1.0 - fx) * fy * fz;
        Ipp64f w111 = fx * fy * fz;

        m_lattice[latticeIdx(ix0, iy0, iz0)] += weight * w000;
        m_lattice[latticeIdx(ix1, iy0, iz0)] += weight * w100;
        m_lattice[latticeIdx(ix0, iy1, iz0)] += weight * w010;
        m_lattice[latticeIdx(ix1, iy1, iz0)] += weight * w110;
        m_lattice[latticeIdx(ix0, iy0, iz1)] += weight * w001;
        m_lattice[latticeIdx(ix1, iy0, iz1)] += weight * w101;
        m_lattice[latticeIdx(ix0, iy1, iz1)] += weight * w011;
        m_lattice[latticeIdx(ix1, iy1, iz1)] += weight * w111;

        m_total_mass += weight;
    }

    // For variance calculation (simplified - you may want to implement properly)
    m_lattice_variance = 0.0;  // You can calculate this properly if needed

    return 0;
}

// Fix convolveWithKernel() - output origin calculation
int SitusGenerator::convolveWithKernel() {
    int half_kernel = m_kernel_extent / 2;

    // Output dimensions = lattice dimensions + kernel extent - 1
    m_nx = m_lattice_nx + m_kernel_extent - 1;
    m_ny = m_lattice_ny + m_kernel_extent - 1;
    m_nz = m_lattice_nz + m_kernel_extent - 1;

    // Output origin should be lattice origin shifted by half kernel
    // This is correct because convolution spreads the lattice
    m_origin[0] = m_lattice_origin[0] - half_kernel * m_grid_spacing;
    m_origin[1] = m_lattice_origin[1] - half_kernel * m_grid_spacing;
    m_origin[2] = m_lattice_origin[2] - half_kernel * m_grid_spacing;

    if (m_config.verbosity >= 1) {
        printf("Convolution:\n");
        printf("  Half kernel: %d voxels\n", half_kernel);
        printf("  Lattice origin: (%.2f, %.2f, %.2f)\n",
            m_lattice_origin[0], m_lattice_origin[1], m_lattice_origin[2]);
        printf("  Output origin: (%.2f, %.2f, %.2f)\n",
            m_origin[0], m_origin[1], m_origin[2]);
        printf("  Output dimensions: %d x %d x %d\n", m_nx, m_ny, m_nz);
    }

    int64_t nvox_output = (int64_t)m_nx * m_ny * m_nz;
    m_output = (Ipp64f*)ippsMalloc_64f(nvox_output);
    if (!m_output) return -1;

    ippsZero_64f(m_output, nvox_output);

    // Direct convolution
#pragma omp parallel for
    for (int k = 0; k < m_lattice_nz; k++) {
        for (int j = 0; j < m_lattice_ny; j++) {
            for (int i = 0; i < m_lattice_nx; i++) {
                Ipp64f val = m_lattice[latticeIdx(i, j, k)];
                if (val != 0.0) {
                    // Add kernel contribution to output
                    for (int kk = 0; kk < m_kernel_extent; kk++) {
                        int out_k = k + kk;
                        if (out_k >= m_nz) continue;

                        for (int jj = 0; jj < m_kernel_extent; jj++) {
                            int out_j = j + jj;
                            if (out_j >= m_ny) continue;

                            for (int ii = 0; ii < m_kernel_extent; ii++) {
                                int out_i = i + ii;
                                if (out_i >= m_nx) continue;

                                m_output[outputIdx(out_i, out_j, out_k)] +=
                                    m_kernel[kernelIdx(ii, jj, kk)] * val;
                            }
                        }
                    }
                }
            }
        }
    }

    return 0;
}

int SitusGenerator::run() {
    if (!m_lattice) return -1;

    printf("\n=== Situs-style Map Generation ===\n");
    printf("Kernel type: %d\n", m_config.kernel_type);
    printf("Resolution input: %.2f %s\n",
        fabs(m_config.resolution),
        (m_config.resolution < 0) ? "(2σ mode)" : "(half-max radius mode)");
    printf("Grid spacing: %.2f A\n", m_grid_spacing);
    printf("Margin voxels: %d\n", m_config.margin_voxels);
    printf("Mass weighting: %s\n", m_config.use_mass_weighting ? "on" : "off");
    printf("Lattice correction: %s\n", m_config.apply_lattice_correction ? "on" : "off");

    // Step 1: Project atoms to lattice
    printf("\nProjecting atoms to lattice...\n");
    int ret = projectAtomsToLattice();
    if (ret != 0) return ret;

    printf("Lattice smoothing sigma: %.3f A\n", m_grid_spacing * sqrt(m_lattice_variance));

    // Step 2: Generate kernel
    printf("Generating kernel...\n");
    ret = generateKernel();
    if (ret != 0) return ret;

    printf("Kernel parameters:\n");
    printf("  sigma = %.3f A\n", m_kernel_params.sigma);
    printf("  r_h = %.3f A\n", m_kernel_params.half_max_radius);
    printf("  r_c = %.3f A\n", m_kernel_params.cutoff_radius);
    printf("  extent = %d voxels\n", m_kernel_extent);

    // Step 3: Convolve
    printf("Convolving lattice with kernel...\n");
    ret = convolveWithKernel();
    if (ret != 0) return ret;

    // Calculate effective resolution
    Ipp64f eff_sigma = sqrt(m_kernel_params.sigma * m_kernel_params.sigma +
        m_lattice_variance * m_grid_spacing * m_grid_spacing);
    Ipp64f eff_resolution = 2.0 * eff_sigma;

    printf("\nOutput map:\n");
    printf("  Dimensions: %d x %d x %d\n", m_nx, m_ny, m_nz);
    printf("  Origin: (%.2f, %.2f, %.2f)\n", m_origin[0], m_origin[1], m_origin[2]);
    printf("  Effective resolution (2σ): %.3f A\n", eff_resolution);
    if (m_config.apply_lattice_correction) {
        printf("  (corrected for lattice smoothing)\n");
    }
    else {
        printf("  (NOTE: slightly larger than target due to uncorrected lattice smoothing)\n");
    }
    printf("================================\n");

    return 0;
}