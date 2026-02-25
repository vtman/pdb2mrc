// src/chimerax_generator.cpp
#include "chimerax_generator.hpp"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <vector>

ChimeraXGenerator::ChimeraXGenerator() {
    m_output = nullptr;
    nx = ny = nz = 0;
    m_origin[0] = m_origin[1] = m_origin[2] = 0.0;
    m_step = 0.0;
    m_cutoff_range = 5.0f;
    m_sigma = 0.0f;
}

ChimeraXGenerator::~ChimeraXGenerator() {
    if (m_output) {
        ippsFree(m_output);
        m_output = nullptr;
    }
}

int ChimeraXGenerator::init(const Atom* atoms, int n_atoms,
    Ipp64f resolution, Ipp64f grid_spacing,
    Ipp64f padding, float cutoff_range) {
    if (!atoms || n_atoms <= 0 || resolution <= 0.0) return -1;

    m_cutoff_range = cutoff_range;
    m_step = (grid_spacing > 0) ? grid_spacing : resolution / 3.0;

    // ChimeraX formula: σ = resolution * sigma_factor where sigma_factor = 1/(π√2)
    m_sigma = (float)(resolution / (M_PI * M_SQRT2));

    // Calculate bounding box with padding
    Ipp64f xmin = atoms[0].x, xmax = atoms[0].x;
    Ipp64f ymin = atoms[0].y, ymax = atoms[0].y;
    Ipp64f zmin = atoms[0].z, zmax = atoms[0].z;

    for (int i = 1; i < n_atoms; i++) {
        if (atoms[i].x < xmin) xmin = atoms[i].x;
        if (atoms[i].x > xmax) xmax = atoms[i].x;
        if (atoms[i].y < ymin) ymin = atoms[i].y;
        if (atoms[i].y > ymax) ymax = atoms[i].y;
        if (atoms[i].z < zmin) zmin = atoms[i].z;
        if (atoms[i].z > zmax) zmax = atoms[i].z;
    }

    // Add padding (ChimeraX uses 3*resolution as default)
    Ipp64f pad = (padding > 0) ? padding : 3.0 * resolution;
    xmin -= pad;
    xmax += pad;
    ymin -= pad;
    ymax += pad;
    zmin -= pad;
    zmax += pad;

    // Calculate grid dimensions
    nx = (int)ceil((xmax - xmin) / m_step) + 1;
    ny = (int)ceil((ymax - ymin) / m_step) + 1;
    nz = (int)ceil((zmax - zmin) / m_step) + 1;

    // Make dimensions even if needed (ChimeraX often does this for FFT)
    if (nx % 2) nx++;
    if (ny % 2) ny++;
    if (nz % 2) nz++;

    // Adjust origin to exact voxel positions
    m_origin[0] = xmin;
    m_origin[1] = ymin;
    m_origin[2] = zmin;

    // Convert atom coordinates to grid indices (ijk space)
    m_centers.resize(n_atoms * 3);
    m_weights.resize(n_atoms);
    m_sdevs.resize(n_atoms * 3);

    for (int i = 0; i < n_atoms; i++) {
        // Convert to grid coordinates (like ChimeraX's xyz_to_ijk_transform)
        m_centers[i * 3 + 0] = (float)((atoms[i].x - m_origin[0]) / m_step);
        m_centers[i * 3 + 1] = (float)((atoms[i].y - m_origin[1]) / m_step);
        m_centers[i * 3 + 2] = (float)((atoms[i].z - m_origin[2]) / m_step);

        // Atomic number as weight (like ChimeraX's element_numbers)
        int Z = AtomUtils::get_atomic_number(atoms[i].element);
        m_weights[i] = (Z > 0) ? (float)Z : 6.0f;  // Default to carbon

        // Normalized sigma per dimension (σ/step)
        float inv_step = 1.0f / (float)m_step;
        m_sdevs[i * 3 + 0] = (float)(m_sigma * inv_step);
        m_sdevs[i * 3 + 1] = (float)(m_sigma * inv_step);
        m_sdevs[i * 3 + 2] = (float)(m_sigma * inv_step);
    }

    return 0;
}

// Direct implementation of ChimeraX's sum_of_gaussians algorithm
static void sum_of_gaussians_core(const std::vector<float>& centers,
    const std::vector<float>& weights,
    const std::vector<float>& sdevs,
    float maxrange,
    float* matrix,
    int nx, int ny, int nz,
    int64_t ms0, int64_t ms1, int64_t ms2) {

    int64_t n = weights.size();

    for (int64_t c = 0; c < n; ++c) {
        float sd[3] = { sdevs[c * 3], sdevs[c * 3 + 1], sdevs[c * 3 + 2] };
        if (sd[0] == 0 || sd[1] == 0 || sd[2] == 0)
            continue;

        float cijk[3] = { centers[c * 3], centers[c * 3 + 1], centers[c * 3 + 2] };

        // Calculate bounds - exactly like ChimeraX
        int ijk_min[3], ijk_max[3];
        ijk_min[0] = (int)ceil(cijk[0] - maxrange * sd[0]);
        ijk_max[0] = (int)floor(cijk[0] + maxrange * sd[0]);
        ijk_min[1] = (int)ceil(cijk[1] - maxrange * sd[1]);
        ijk_max[1] = (int)floor(cijk[1] + maxrange * sd[1]);
        ijk_min[2] = (int)ceil(cijk[2] - maxrange * sd[2]);
        ijk_max[2] = (int)floor(cijk[2] + maxrange * sd[2]);

        // Clamp to grid boundaries
        ijk_min[0] = (ijk_min[0] < 0) ? 0 : ijk_min[0];
        ijk_min[1] = (ijk_min[1] < 0) ? 0 : ijk_min[1];
        ijk_min[2] = (ijk_min[2] < 0) ? 0 : ijk_min[2];
        ijk_max[0] = (ijk_max[0] >= nx) ? nx - 1 : ijk_max[0];
        ijk_max[1] = (ijk_max[1] >= ny) ? ny - 1 : ijk_max[1];
        ijk_max[2] = (ijk_max[2] >= nz) ? nz - 1 : ijk_max[2];

        float cf = weights[c];

        // Nested loops ordered for cache efficiency (k,j,i) - exactly like ChimeraX
        for (int k = ijk_min[2]; k <= ijk_max[2]; ++k) {
            float dk = (k - cijk[2]) / sd[2];
            float k2 = dk * dk;

            for (int j = ijk_min[1]; j <= ijk_max[1]; ++j) {
                float dj = (j - cijk[1]) / sd[1];
                float jk2 = dj * dj + k2;

                int64_t base_idx = k * ms0 + j * ms1;

                for (int i = ijk_min[0]; i <= ijk_max[0]; ++i) {
                    float di = (i - cijk[0]) / sd[0];
                    float ijk2 = di * di + jk2;
                    matrix[base_idx + i * ms2] += cf * (float)exp(-0.5 * ijk2);
                }
            }
        }
    }
}

int ChimeraXGenerator::run() {
    if (m_centers.empty()) return -1;

    int64_t nvox = (int64_t)nx * ny * nz;

    // Allocate output grid
    m_output = (Ipp64f*)ippsMalloc_64f(nvox);
    if (!m_output) return -2;

    ippsZero_64f(m_output, nvox);

    // Strides for 3D array access (like ChimeraX's matrix.stride())
    int64_t ms0 = (int64_t)ny * nx;  // stride for k (slowest)
    int64_t ms1 = nx;                  // stride for j
    int64_t ms2 = 1;                   // stride for i (fastest)

    printf("ChimeraX map generation:\n");
    printf("  Grid: %d x %d x %d (%lld voxels)\n", nx, ny, nz, nvox);
    printf("  σ = %.3f Å, cutoff = %.1fσ\n", m_sigma, m_cutoff_range);
    printf("  Atoms: %zu\n", m_weights.size());

    // Convert to float for processing (ChimeraX uses float)
    std::vector<float> float_output(nvox, 0.0f);

    // Process in parallel (ChimeraX releases GIL but doesn't use OpenMP)
#pragma omp parallel for
    for (int64_t c = 0; c < (int64_t)m_weights.size(); ++c) {
        // Each thread needs its own copy of bounds calculation
        // But we're using the core function that's thread-safe
    }

    // For simplicity, call the core function (could be parallelized)
    sum_of_gaussians_core(m_centers, m_weights, m_sdevs, m_cutoff_range,
        float_output.data(), nx, ny, nz,
        ms0, ms1, ms2);

    // Apply normalization (exactly like ChimeraX)
    // normalization = pow(2*pi,-1.5) * pow(sdev,-3)
    float norm = (float)(pow(2.0 * M_PI, -1.5) * pow(m_sigma, -3.0));

#pragma omp parallel for
    for (int64_t i = 0; i < nvox; i++) {
        m_output[i] = float_output[i] * norm;
    }

    // Normalize to max = 1.0 (ChimeraX does this for display)
    Ipp64f max_val;
    ippsMax_64f(m_output, nvox, &max_val);
    if (max_val > 0) {
        ippsMulC_64f_I(1.0 / max_val, m_output, nvox);
    }

    return 0;
}