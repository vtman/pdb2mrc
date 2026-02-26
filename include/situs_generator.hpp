// include/situs_generator.hpp
#pragma once

#include "types.hpp"
#include "atom.hpp"
#include <vector>
#include <cmath>

enum SitusKernelType {
    SITUS_KERNEL_GAUSSIAN = 1,
    SITUS_KERNEL_TRIANGULAR = 2,
    SITUS_KERNEL_SEMI_EPANECHNIKOV = 3,
    SITUS_KERNEL_EPANECHNIKOV = 4,
    SITUS_KERNEL_HARD_SPHERE = 5,
    SITUS_KERNEL_DEFAULT = SITUS_KERNEL_GAUSSIAN
};


struct SitusGeneratorConfig {
    SitusKernelType kernel_type;
    Ipp64f resolution;           // Target resolution in Å (positive = half-max radius, negative = 2σ)
    Ipp64f grid_spacing;         // Voxel size in Å
    Ipp64f padding;              // Extra margin in voxels (default 2)
    int use_mass_weighting;       // Use atomic masses vs unit weights
    int apply_lattice_correction; // Apply variance correction
    int margin_voxels;            // Extra voxels around structure (default 2)
    Ipp64f kernel_amplitude;      // Scaling factor for kernel
    int verbosity;                // Add this for debug output

    SitusGeneratorConfig() {
        kernel_type = SITUS_KERNEL_GAUSSIAN;
        resolution = 10.0;
        grid_spacing = 2.0;
        padding = 0.0;
        use_mass_weighting = 1;
        apply_lattice_correction = 1;
        margin_voxels = 2;
        kernel_amplitude = 1.0;
        verbosity = 0;  // Default to 0
    }
};

// Kernel parameters for each type
struct SitusKernelParams {
    Ipp64f sigma;          // Gaussian sigma (1D)
    Ipp64f half_max_radius; // r_h
    Ipp64f cutoff_radius;   // r_c
};

class SitusGenerator {
private:
    // Configuration
    SitusGeneratorConfig m_config;

    // Atom data
    const Atom* m_atoms;
    int m_nAtoms;

    // Grid parameters
    int m_nx, m_ny, m_nz;           // Output grid dimensions
    int m_lattice_nx, m_lattice_ny, m_lattice_nz; // Lattice projection dimensions
    Ipp64f m_origin[3];              // Output origin
    Ipp64f m_lattice_origin[3];       // Lattice projection origin
    Ipp64f m_grid_spacing;

    // Data buffers
    Ipp64f* m_lattice;               // Projection lattice
    Ipp64f* m_kernel;                 // Smoothing kernel
    Ipp64f* m_output;                 // Final output map

    // Kernel parameters
    SitusKernelParams m_kernel_params;
    int m_kernel_extent;               // Kernel size in voxels (K)

    // Statistics
    Ipp64f m_lattice_variance;
    Ipp64f m_total_mass;

    // Private helper methods
    int computeBoundingBox(const Atom* atoms, int n_atoms,
        Ipp64f& xmin, Ipp64f& xmax,
        Ipp64f& ymin, Ipp64f& ymax,
        Ipp64f& zmin, Ipp64f& zmax);

    void computeKernelParams();
    int generateKernel();
    int projectAtomsToLattice();  // No parameters - uses m_atoms, m_nAtoms
    int convolveWithKernel();

    // Kernel generation functions
    void generateGaussianKernel(Ipp64f* kernel, int extent, Ipp64f sigma, Ipp64f amplitude);
    void generateTriangularKernel(Ipp64f* kernel, int extent, Ipp64f rh, Ipp64f amplitude);
    void generateSemiEpanechnikovKernel(Ipp64f* kernel, int extent, Ipp64f rh, Ipp64f amplitude);
    void generateEpanechnikovKernel(Ipp64f* kernel, int extent, Ipp64f rh, Ipp64f amplitude);
    void generateHardSphereKernel(Ipp64f* kernel, int extent, Ipp64f rh, Ipp64f amplitude);

    // Index calculations
    inline int64_t latticeIdx(int i, int j, int k) const {
        return ((int64_t)k * m_lattice_ny + j) * m_lattice_nx + i;
    }

    inline int64_t kernelIdx(int i, int j, int k) const {
        return ((int64_t)k * m_kernel_extent + j) * m_kernel_extent + i;
    }

    inline int64_t outputIdx(int i, int j, int k) const {
        return ((int64_t)k * m_ny + j) * m_nx + i;
    }

public:
    SitusGenerator();
    ~SitusGenerator();

    // Initialize with configuration
    int init(const SitusGeneratorConfig* config, const Atom* atoms, int n_atoms);

    // Run map generation
    int run();

    // Get results
    int getOutputNx() const { return m_nx; }
    int getOutputNy() const { return m_ny; }
    int getOutputNz() const { return m_nz; }
    void getOutputOrigin(Ipp64f origin[3]) const {
        origin[0] = m_origin[0];
        origin[1] = m_origin[1];
        origin[2] = m_origin[2];
    }
    Ipp64f getGridSpacing() const { return m_grid_spacing; }

    // Transfer ownership of output map
    Ipp64f* releaseOutput() {
        Ipp64f* map = m_output;
        m_output = nullptr;
        return map;
    }

    // Cleanup
    void cleanup();

    // Disable copy
    SitusGenerator(const SitusGenerator&) = delete;
    SitusGenerator& operator=(const SitusGenerator&) = delete;
};