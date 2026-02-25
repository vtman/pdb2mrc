// include/emmer_generator.hpp
#pragma once

#include "types.hpp"
#include "atom.hpp"
#include "enums.hpp"
#include <mkl.h>  // For DFTI_DESCRIPTOR_HANDLE

#define EMMER_MAX_ELEMENTS 98
#define EMMER_N_GAUSS 5

// EMmer-specific configuration
struct EmmerGeneratorConfig {
    int model_index;           // Which model to use (for multi-model PDB)
    int align_output;           // Apply MRC alignment (flip+rotate)
    int set_refmac_blur;        // Apply Refmac-compatible blur
    double blur;                 // Manual blur value (if >0 overrides auto)
    int symmetry_expansion;      // Apply space group symmetry
    double cutoff_level;         // Density cutoff for radius determination
    double rate;                 // Shannon rate for grid spacing (default 1.5)

    EmmerGeneratorConfig();
};

// 5-Gaussian coefficient structure
struct EmmerGaussianCoeff {
    char element[4];
    double a[EMMER_N_GAUSS];
    double b[EMMER_N_GAUSS];
    int atomic_number;

    double fe0() const;
};

// Precalculated Gaussian for an element at given B-factor
struct EmmerPrecalcGaussian {
    double amplitude;
    double width;
    double radius;
};

// Per-element data (no STL vectors - use arrays)
struct EmmerElementData {
    int element_idx;           // Index in COEFFICIENTS
    int atom_count;            // Number of atoms of this element
    int* atom_indices;         // Dynamically allocated array
    EmmerPrecalcGaussian precalc[EMMER_N_GAUSS];
    double avg_b_factor;
};

// Helper class for FFT operations
struct EmmerFFT {
    int nx, ny, nz;
    int64_t nreals, ncomplex;
    Ipp64f* real_in;
    Ipp64fc* complex_out;
    Ipp64fc* kernel;
    DFTI_DESCRIPTOR_HANDLE fft_desc;

    void init();
    void cleanup();
    int setup(int nx_, int ny_, int nz_);
    int forward(Ipp64f* input);
    int applyUnblur(double blur, double grid_spacing);
    int inverse(Ipp64f* output);
};

class EmmerGenerator {
private:
    // Configuration
    EmmerGeneratorConfig m_config;

    // Atom data
    const Atom* m_atoms;
    int m_nAtoms;

    // Grid parameters
    int m_nx, m_ny, m_nz;
    Ipp64f m_origin[3];
    Ipp64f m_grid_spacing;
    Ipp64f m_apix[3];
    double m_d_min;

    // Output map
    Ipp64f* m_output;
    EmmerFFT m_fft;

    // Precomputed Gaussian coefficients
    static const EmmerGaussianCoeff COEFFICIENTS[];
    static const int N_COEFFICIENTS;

    // Per-element data (fixed array, m_nElements tracks count)
    EmmerElementData m_element_data[EMMER_MAX_ELEMENTS];
    int m_nElements;

    // Working buffer
    Ipp64f* m_work_grid;

    // Private methods
    int findUniqueElements();
    int resolveParameters(const Ipp64f* cell_dimensions);
    double calculateRefmacBlur();
    void precalculateForElement(EmmerElementData* elem, double B_eff);
    double estimateRadius(const EmmerPrecalcGaussian* gaussians);
    void addAtomToGrid(int atom_idx, const EmmerElementData* elem, Ipp64f* grid);
    void applyOutputAlignment();
    void symmetrizeGrid();

public:
    EmmerGenerator();
    ~EmmerGenerator();

    int init(const EmmerGeneratorConfig* config,
        const Atom* atoms,
        int n_atoms,
        Ipp64f resolution,
        Ipp64f grid_spacing,
        const Ipp64f* cell_dimensions);

    int run();

    // Get results
    int getOutputNx() const { return m_nx; }
    int getOutputNy() const { return m_ny; }
    int getOutputNz() const { return m_nz; }
    void getOutputOrigin(Ipp64f origin[3]) const;
    Ipp64f getGridSpacing() const { return m_grid_spacing; }

    // Transfer ownership
    Ipp64f* releaseOutput();

    // Cleanup
    void cleanup();

    // Disable copy
    EmmerGenerator(const EmmerGenerator&) = delete;
    EmmerGenerator& operator=(const EmmerGenerator&) = delete;
};