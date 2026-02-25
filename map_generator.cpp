// src/map_generator.cpp
#include "map_generator.hpp"
#include <stdio.h>
#include <string.h>
#include <math.h>

//=============================================================================
// MapParams implementation
//=============================================================================

MapParams::MapParams() {
    resolution = 6.0;
    criterion = CRITERION_RAYLEIGH;
    amplitude_mode = AMPLITUDE_PENG1996;  // Default to Peng 1996
    grid_spacing = 0.0;
    padding = 0.0;
    cutoff_level = 0.001;
    use_bfactors = 0;
    b_default = 20.0;
    cutoff_range = 5.0;
}

//=============================================================================
// ElementWorkspace implementation
//=============================================================================

ElementWorkspace::ElementWorkspace() {
    grid = nullptr;
    fft_buffer = nullptr;
    kernel_buffer = nullptr;
    fft_desc = nullptr;
}

ElementWorkspace::~ElementWorkspace() {
    cleanup();
}

void ElementWorkspace::cleanup() {
    if (grid) {
        ippsFree(grid);
        grid = nullptr;
    }
    if (fft_buffer) {
        ippsFree(fft_buffer);
        fft_buffer = nullptr;
    }
    if (kernel_buffer) {
        ippsFree(kernel_buffer);
        kernel_buffer = nullptr;
    }
    if (fft_desc) {
        DftiFreeDescriptor(&fft_desc);
        fft_desc = nullptr;
    }
}

//=============================================================================
// MapGenerator implementation
//=============================================================================

MapGenerator::MapGenerator() {
    m_profiles = nullptr;
    nProfiles = 0;
    m_element_to_profile = nullptr;

    m_atoms = nullptr;
    nAtoms = 0;
    vStartAtom = nullptr;
    nUnique = 0;

    m_elem_workspace = nullptr;
    m_n_elem_workspace = 0;

    vOutput = nullptr;

    m_out_origin[0] = 0.0;
    m_out_origin[1] = 0.0;
    m_out_origin[2] = 0.0;

    // Initialize FFT buffers to nullptr
    vData = nullptr;
    vcData = nullptr;
    vKernel = nullptr;
    vcKernel = nullptr;
    descFor = nullptr;
    descBack = nullptr;

    // Initialize buffer sizes
    nt = 0;
    nct = 0;

    // Initialize grid dimensions
    nx = ny = nz = 0;
    m_xyz_offset = 0;
    m_origin[0] = m_origin[1] = m_origin[2] = 0.0;
    m_spacing = 0.0;

    // Initialize scattering table
    Scattering::init(&m_table);
}

MapGenerator::~MapGenerator() {
    if (descFor != nullptr) { DftiFreeDescriptor(&descFor); descFor = nullptr; }
    if (descBack != nullptr) { DftiFreeDescriptor(&descBack); descBack = nullptr; }
    cleanup();
}

void MapGenerator::cleanup() {
    // Free element-to-profile mapping
    if (m_element_to_profile) {
        ippsFree(m_element_to_profile);
        m_element_to_profile = nullptr;
    }

    // Free unique elements array
    if (vStartAtom) {
        ippsFree(vStartAtom);
        vStartAtom = nullptr;
    }
    nUnique = 0;

    // Free profiles
    if (m_profiles) {
        Scattering::free_profiles(m_profiles, nProfiles);
        m_profiles = nullptr;
        nProfiles = 0;
    }

    // Free element workspace
    if (m_elem_workspace) {
        for (int i = 0; i < m_n_elem_workspace; i++) {
            m_elem_workspace[i].cleanup();
        }
        ippsFree(m_elem_workspace);
        m_elem_workspace = nullptr;
        m_n_elem_workspace = 0;
    }

    // Free FFT buffers
    if (vData != nullptr) { mkl_free(vData); vData = nullptr; }
    if (vcData != nullptr) { mkl_free(vcData); vcData = nullptr; }
    if (vKernel != nullptr) { mkl_free(vKernel); vKernel = nullptr; }
    if (vcKernel != nullptr) { mkl_free(vcKernel); vcKernel = nullptr; }
    if (descFor != nullptr) { DftiFreeDescriptor(&descFor); descFor = nullptr; }
    if (descBack != nullptr) { DftiFreeDescriptor(&descBack); descBack = nullptr; }

    // Free output map (if not already released)
    if (vOutput != nullptr) { ippsFree(vOutput); vOutput = nullptr; }

    // Reset atom pointers (we don't own them)
    m_atoms = nullptr;
    nAtoms = 0;
}

//=============================================================================
// Private helper methods
//=============================================================================

int MapGenerator::findUniqueElements() {
    if (!m_atoms || nAtoms <= 0) return -1;

    int element_present[MAX_ELEMENTS] = { 0 };
    int temp_unique[MAX_ELEMENTS];
    nUnique = 0;

    for (int i = 0; i < nAtoms; i++) {
        int idx = m_atoms[i].element_idx;
        if (idx >= 0 && idx < m_table.n_entries && !element_present[idx]) {
            element_present[idx] = 1;
            temp_unique[nUnique++] = idx;
        }
    }

    if (nUnique == 0) return -2;

    // Allocate and copy unique elements
    vStartAtom = (int*)ippsMalloc_32s(nUnique);
    if (!vStartAtom) return -3;

    for (int i = 0; i < nUnique; i++) {
        vStartAtom[i] = temp_unique[i];
    }

    return 0;
}

int MapGenerator::calculateGridDimensions() {
    if (!m_atoms || nAtoms <= 0) return -1;

    // Find global Rmax across all profiles
    Ipp64f global_rmax = 0.0;
    for (int i = 0; i < nProfiles; i++) {
        if (m_profiles[i].rmax > global_rmax) {
            global_rmax = m_profiles[i].rmax;
        }
    }

    // Add extra padding
    global_rmax += m_params.padding;

    // Calculate bounding box from atom positions only
    Ipp64f xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = xmax = m_atoms[0].x;
    ymin = ymax = m_atoms[0].y;
    zmin = zmax = m_atoms[0].z;

    for (int i = 1; i < nAtoms; i++) {
        if (m_atoms[i].x < xmin) xmin = m_atoms[i].x;
        if (m_atoms[i].x > xmax) xmax = m_atoms[i].x;
        if (m_atoms[i].y < ymin) ymin = m_atoms[i].y;
        if (m_atoms[i].y > ymax) ymax = m_atoms[i].y;
        if (m_atoms[i].z < zmin) zmin = m_atoms[i].z;
        if (m_atoms[i].z > zmax) zmax = m_atoms[i].z;
    }


    m_xyz_offset = (int)ceil(global_rmax / m_spacing) + 1;

    Ipp64f dOffset;

    dOffset = (double)m_xyz_offset * m_spacing;

    // Expand bounding box by global_rmax
    xmin -= dOffset;
    xmax += dOffset;
    ymin -= dOffset;
    ymax += dOffset;
    zmin -= dOffset;
    zmax += dOffset;

    // Calculate number of grid points
    nx = 4 * (((int)ceil((xmax - xmin) / m_spacing) + 3) / 4);
    ny = 4 * (((int)ceil((ymax - ymin) / m_spacing) + 3) / 4);
    nz = 4 * (((int)ceil((zmax - zmin) / m_spacing) + 3) / 4);

    // Adjust max to be exactly (n-1)*spacing from min
    xmax = xmin + (nx - 1) * m_spacing;
    ymax = ymin + (ny - 1) * m_spacing;
    zmax = zmin + (nz - 1) * m_spacing;

    m_origin[0] = xmin;
    m_origin[1] = ymin;
    m_origin[2] = zmin;

    // Calculate buffer sizes
    nt = (int64_t)nx * ny * nz;
    nct = (int64_t)nz * ny * (nx / 2 + 1);

    return 0;
}

int MapGenerator::allocateFFTBuffers() {
    // Allocate real buffer
    vData = (Ipp64f*)mkl_malloc(nt * sizeof(Ipp64f), 64);
    if (vData == nullptr) return -1;

    // Allocate complex buffer
    vcData = (Ipp64fc*)mkl_malloc(nct * sizeof(Ipp64fc), 64);
    if (vcData == nullptr) return -2;

    // Allocate kernel real buffer
    vKernel = (Ipp64f*)mkl_malloc(nt * sizeof(Ipp64f), 64);
    if (vKernel == nullptr) return -3;

    // Allocate kernel complex buffer
    vcKernel = (Ipp64fc*)mkl_malloc(nct * sizeof(Ipp64fc), 64);
    if (vcKernel == nullptr) return -4;

    // All allocations succeeded
    return 0;
}

int MapGenerator::setupFFTDescriptor() {
    MKL_LONG status;
    MKL_LONG dims[3] = { nz, ny, nx };

    // Strides for real data (row-major order: [z][y][x])
    MKL_LONG rstrides[4] = { 0, ny * nx, nx, 1 };

    // Strides for complex data (conjugate-even storage)
    MKL_LONG cstrides[4] = { 0, ny * (nx / 2 + 1), (nx / 2 + 1), 1 };


    status = DftiCreateDescriptor(&descFor, DFTI_DOUBLE, DFTI_REAL, 3, dims);
    if (status != 0) {
        printf("Error: DftiCreateDescriptor failed: %s\n", DftiErrorMessage(status));
        return -1;
    }
    status = DftiCreateDescriptor(&descBack, DFTI_DOUBLE, DFTI_REAL, 3, dims);
    if (status != 0) {
        printf("Error: DftiCreateDescriptor failed: %s\n", DftiErrorMessage(status));
        return -1;
    }

    status = DftiSetValue(descFor, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != 0) {
        printf("Error: DFTI_PLACEMENT failed: %s\n", DftiErrorMessage(status));
        return -2;
    }
    status = DftiSetValue(descBack, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    if (status != 0) {
        printf("Error: DFTI_PLACEMENT failed: %s\n", DftiErrorMessage(status));
        return -2;
    }

    status = DftiSetValue(descFor, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != 0) {
        printf("Error: DFTI_CONJUGATE_EVEN_STORAGE failed: %s\n", DftiErrorMessage(status));
        return -3;
    }
    status = DftiSetValue(descBack, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    if (status != 0) {
        printf("Error: DFTI_CONJUGATE_EVEN_STORAGE failed: %s\n", DftiErrorMessage(status));
        return -3;
    }

    status = DftiSetValue(descFor, DFTI_INPUT_STRIDES, rstrides);
    if (status != 0) {
        printf("Error: DFTI_INPUT_STRIDES failed: %s\n", DftiErrorMessage(status));
        return -4;
    }
    status = DftiSetValue(descBack, DFTI_OUTPUT_STRIDES, rstrides);
    if (status != 0) {
        printf("Error: DFTI_OUTPUT_STRIDES failed: %s\n", DftiErrorMessage(status));
        return -4;
    }

    status = DftiSetValue(descFor, DFTI_OUTPUT_STRIDES, cstrides);
    if (status != 0) {
        printf("Error: DFTI_OUTPUT_STRIDES failed: %s\n", DftiErrorMessage(status));
        return -5;
    }
    status = DftiSetValue(descBack, DFTI_INPUT_STRIDES, cstrides);
    if (status != 0) {
        printf("Error: DFTI_INPUT_STRIDES failed: %s\n", DftiErrorMessage(status));
        return -5;
    }

    status = DftiCommitDescriptor(descFor);
    if (status != 0) {
        printf("Error: DftiCommitDescriptor failed: %s\n", DftiErrorMessage(status));
        return -7;
    }
    status = DftiCommitDescriptor(descBack);
    if (status != 0) {
        printf("Error: DftiCommitDescriptor failed: %s\n", DftiErrorMessage(status));
        return -7;
    }

    return 0;
}

void MapGenerator::projectElementAtoms(int element_idx, Ipp64f* grid, int xyz_offset) {
    Ipp64f inv_spacing = 1.0 / m_spacing;

    for (int i = 0; i < nAtoms; i++) {
        if (m_atoms[i].element_idx != element_idx) continue;

        // Convert to grid coordinates
        Ipp64f gx = (m_atoms[i].x - m_origin[0]) * inv_spacing;
        Ipp64f gy = (m_atoms[i].y - m_origin[1]) * inv_spacing;
        Ipp64f gz = (m_atoms[i].z - m_origin[2]) * inv_spacing;

        int ix0 = (int)floor(gx);
        int iy0 = (int)floor(gy);
        int iz0 = (int)floor(gz);

        // Check bounds (we're in padded region so should be safe)
        if (ix0 < 0 || iy0 < 0 || iz0 < 0) continue;

        // Fractional parts
        Ipp64f fx = gx - ix0;
        Ipp64f fy = gy - iy0;
        Ipp64f fz = gz - iz0;

        // Trilinear weights
        Ipp64f w000 = (1.0 - fx) * (1.0 - fy) * (1.0 - fz);
        Ipp64f w100 = fx * (1.0 - fy) * (1.0 - fz);
        Ipp64f w010 = (1.0 - fx) * fy * (1.0 - fz);
        Ipp64f w110 = fx * fy * (1.0 - fz);
        Ipp64f w001 = (1.0 - fx) * (1.0 - fy) * fz;
        Ipp64f w101 = fx * (1.0 - fy) * fz;
        Ipp64f w011 = (1.0 - fx) * fy * fz;
        Ipp64f w111 = fx * fy * fz;

        Ipp64f weight = m_atoms[i].occupancy;

        // Add to 8 surrounding voxels (with offset for padding)
        int64_t base = (iz0 * ny + iy0) * nx + ix0;

        grid[base] += weight * w000;
        grid[base + 1] += weight * w100;

        int64_t row = base + nx;
        grid[row] += weight * w010;
        grid[row + 1] += weight * w110;

        int64_t slice = base + (int64_t)nx * ny;
        grid[slice] += weight * w001;
        grid[slice + 1] += weight * w101;

        int64_t slice_row = slice + nx;
        grid[slice_row] += weight * w011;
        grid[slice_row + 1] += weight * w111;
    }
}

void MapGenerator::createElementKernel(int element_idx, Ipp64f* kernel, int kernel_radius_voxels) {
    // Find profile index for this element
    int profile_idx = m_element_to_profile[element_idx];
    if (profile_idx < 0) return;

    ElementProfile* profile = &m_profiles[profile_idx];

    for (int k = -kernel_radius_voxels; k <= kernel_radius_voxels; k++) {
        int kk = k;
        if (kk < 0) kk += nz;
        Ipp64f dz2 = (double)(k * k);

        for (int j = -kernel_radius_voxels; j <= kernel_radius_voxels; j++) {
            int jj = j;
            if (jj < 0) jj += ny;
            Ipp64f dyz2 = dz2 + (double)(j * j);

            for (int i = -kernel_radius_voxels; i <= kernel_radius_voxels; i++) {
                int ii = i;
                if (ii < 0) ii += nx;
                Ipp64f r = m_spacing * sqrt((double)(i * i) + dyz2);

                int64_t idx = ((int64_t)kk * ny + jj) * nx + ii;
                kernel[idx] = Scattering::profile_at_r(profile, r);
            }
        }
    }
}

void MapGenerator::addToFinalMap(const Ipp64f* src) {
    ippsAdd_64f_I(src, vOutput, nt);
}

void MapGenerator::normalizeMap() {
    Ipp64f max_val = 0.0;
    ippsMax_64f(vOutput, nt, &max_val);

    if (max_val > 0) {
        Ipp64f inv_max = 1.0 / max_val;
        ippsMulC_64f_I(inv_max, vOutput, nt);
    }
}

void MapGenerator::printMapStatistics() {
    Ipp64f min_val, max_val, mean_val, std_val;

    ippsMinMax_64f(vOutput, nt, &min_val, &max_val);
    ippsMean_64f(vOutput, nt, &mean_val);
    ippsStdDev_64f(vOutput, nt, &std_val);

    printf("Map statistics: min=%.6f, max=%.6f, mean=%.6f, std=%.6f\n",
        min_val, max_val, mean_val, std_val);
}

int MapGenerator::processElement(int element_idx) {

    MKL_LONG status;

    // Clear real buffer for this element
    ippsZero_64f(vData, nt);

    // Project atoms of this element to grid
    projectElementAtoms(element_idx, vData, m_xyz_offset);

    // Clear kernel real buffer
    ippsZero_64f(vKernel, nt);

    // Create kernel
    int profile_idx = m_element_to_profile[element_idx];
    if (profile_idx >= 0) {
        int kernel_radius_voxels = (int)ceil(m_profiles[profile_idx].rmax / m_spacing);
        createElementKernel(element_idx, vKernel, kernel_radius_voxels);
    }

    // Forward FFT of grid (real-to-complex) - out-of-place
    status = DftiComputeForward(descFor, vData, vcData);
    if (status != 0) {
        printf("Warning: Forward FFT failed for grid: %s\n", DftiErrorMessage(status));
        return -1;
    }

    // Forward FFT of kernel (real-to-complex) - out-of-place
    status = DftiComputeForward(descFor, vKernel, vcKernel);
    if (status != 0) {
        printf("Warning: Forward FFT failed for kernel: %s\n", DftiErrorMessage(status));
        return -2;
    }

    ippsMul_64fc_I(vcKernel, vcData, nct);

    // Inverse FFT (complex-to-real)
    status = DftiComputeBackward(descBack, vcData, vData);
    if (status != 0) {
        printf("Warning: Inverse FFT failed: %s\n", DftiErrorMessage(status));
        return -6;
    }

    // Add result to final map
    addToFinalMap(vData);
    return 0;
}

//=============================================================================
// Public methods
//=============================================================================

int MapGenerator::init(const MapParams* params, const Atom* atoms, int n_atoms) {
    if (!params || !atoms || n_atoms <= 0) return -1;

    // Store parameters and atoms
    m_params = *params;
    m_atoms = atoms;
    nAtoms = n_atoms;
    m_spacing = params->grid_spacing;

    // Precompute profiles ONLY for elements present in atoms
    // Pass amplitude mode to precompute_profiles
    int ret = Scattering::precompute_profiles(&m_table,
        params->resolution, params->grid_spacing, params->cutoff_level,
        params->amplitude_mode,  // New parameter
        atoms, n_atoms,
        &m_profiles, &nProfiles);

    if (ret != 0) return ret;

    // Create mapping from element index to profile index
    m_element_to_profile = (int*)ippsMalloc_32s(m_table.n_entries);
    if (!m_element_to_profile) {
        Scattering::free_profiles(m_profiles, nProfiles);
        return -5;
    }

    // Initialize all to -1 (not present)
    ippsSet_32s(-1, m_element_to_profile, m_table.n_entries);

    // Fill mapping for profiles that exist
    for (int i = 0; i < nProfiles; i++) {
        int elem_idx = m_profiles[i].element_idx;
        if (elem_idx >= 0 && elem_idx < m_table.n_entries) {
            m_element_to_profile[elem_idx] = i;
        }
    }

    return 0;
}

int MapGenerator::run() {
    if (!m_atoms || nAtoms <= 0) return -1;

    // Step 1: Find unique elements present in the atoms
    int ret = findUniqueElements();
    if (ret != 0 || nUnique == 0) return -2;

    printf("Found %d unique elements\n", nUnique);
    printf("Amplitude mode: %s\n", amplitude_mode_name(m_params.amplitude_mode));

    // Step 2: Calculate grid dimensions and FFT parameters
    ret = calculateGridDimensions();
    if (ret != 0) return -3;

    // Store output dimensions and origin
    m_out_origin[0] = m_origin[0];
    m_out_origin[1] = m_origin[1];
    m_out_origin[2] = m_origin[2];

    printf("Grid: %d x %d x %d (voxel size: %.2f A)\n",
        nx, ny, nz, m_params.grid_spacing);
    printf("FFT grid: %d x %d x %d\n", nx, ny, nz);
    printf("Total voxels: %lld\n", nt);

    // Step 3: Allocate output map
    vOutput = (Ipp64f*)ippsMalloc_64f(nt);
    if (vOutput == nullptr) return -4;

    ippsZero_64f(vOutput, nt);

    // Step 4: Allocate FFT buffers
    ret = allocateFFTBuffers();
    if (ret != 0) return -5;

    // Step 5: Setup FFT descriptor
    ret = setupFFTDescriptor();
    if (ret != 0) {
        printf("Error: Failed to setup FFT descriptor (error %d)\n", ret);
        return -6;
    }

    // Step 6: Process each element
    for (int e = 0; e < nUnique; e++) {
        int elem_idx = vStartAtom[e];

        printf("Processing element %s...\n", m_table.entries[elem_idx].name);

        ret = processElement(elem_idx);

        if (ret != 0) {
            printf("Warning: Failed to process element %d (error %d)\n", elem_idx, ret);
        }
    }

    // Step 7: Normalize if requested
    normalizeMap();

    // Step 8: Print statistics
    printMapStatistics();

    return 0;
}

//=============================================================================
// Utility functions
//=============================================================================

Ipp64f resolution_to_sigma(Ipp64f resolution, ResolutionCriterion criterion) {
    switch (criterion) {
    case CRITERION_RAYLEIGH: return resolution / 1.665;
    case CRITERION_CHIMERAX: return resolution / (M_PI * M_SQRT2);
    case CRITERION_EMAN2:    return resolution / (M_PI * sqrt(8.0));  // R/(π√8)
    case CRITERION_FSC_0143: return resolution / (1.1 * 1.665);
    case CRITERION_FSC_05:   return resolution / (1.3 * 1.665);
    default: return resolution / 1.665;
    }
}



