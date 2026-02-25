// include/map_generator.hpp
#pragma once

#include "types.hpp"
#include "atom.hpp"
#include "enums.hpp"  // New header for enums
#include "scattering.hpp"  // Now this includes enums.hpp, not map_generator.hpp
#include "grid_projection.hpp"
#include <mkl.h>

// Map generation parameters
struct MapParams {
    Ipp64f resolution;           // Target resolution in Å
    ResolutionCriterion criterion;
    AmplitudeMode amplitude_mode; // Amplitude scaling mode
    Ipp64f grid_spacing;         // Voxel size (0 = auto)
    Ipp64f cutoff_level;         // Drop to this fraction of max (default 0.001)
    Ipp64f padding;
    Ipp64f cutoff_range;
    int use_bfactors;            // Apply B-factors from PDB
    Ipp64f b_default;            // Default B-factor

    MapParams();
};

// Element workspace for FFT processing
struct ElementWorkspace {
    Ipp64f* grid;                // Grid for this element
    Ipp64fc* fft_buffer;         // FFT buffer for grid
    Ipp64fc* kernel_buffer;      // Separate buffer for kernel FFT
    DFTI_DESCRIPTOR_HANDLE fft_desc;  // Single descriptor reused for forward/backward

    ElementWorkspace();
    ~ElementWorkspace();
    void cleanup();
};

// Main map generator class
class MapGenerator {
private:
    // Scattering data
    ScatteringTable m_table;
    ElementProfile* m_profiles;      // Profiles for elements present
    int nProfiles;                 // Number of profiles
    int* m_element_to_profile;        // Maps element_idx -> profile index

    // Atom data (now as members)
    const Atom* m_atoms;              // Pointer to atoms (not owned)
    int nAtoms;                     // Number of atoms
    int* vStartAtom;            // Array of unique element indices
    int nUnique;                     // Number of unique elements

    // Working spaces
    ElementWorkspace* m_elem_workspace;
    int m_n_elem_workspace;

    // Grid parameters
    Ipp64f m_spacing;
    int nx, ny, nz;  // Local grid dimensions
    int m_xyz_offset;                         // Offset for padding
    Ipp64f m_origin[3];

    // Buffer sizes (now as members)
    int64_t nt;                     // Size of real FFT buffer
    int64_t nct;                   // Size of complex FFT buffer

    // Output data (now as members)
    Ipp64f m_out_origin[3];                    // Output origin

    // FFT buffers (as class members)
    Ipp64f* vData, * vKernel, * vOutput;
    Ipp64fc* vcData, * vcKernel;
    DFTI_DESCRIPTOR_HANDLE descFor, descBack;

    // Parameters (stored for internal use)
    MapParams m_params;

    // Private helper methods
    int findUniqueElements();
    int calculateGridDimensions();
    int allocateFFTBuffers();
    int setupFFTDescriptor();
    int processElement(int element_idx);
    void projectElementAtoms(int element_idx, Ipp64f* grid, int xyz_offset);
    void createElementKernel(int element_idx, Ipp64f* kernel, int kernel_radius_voxels);
    void addToFinalMap(const Ipp64f* src);
    void normalizeMap();
    void printMapStatistics();

public:
    MapGenerator();
    ~MapGenerator();

    // Initialize with atoms to determine which elements are present
    int init(const MapParams* params, const Atom* atoms, int n_atoms);

    // Generate map - now no output parameters, results are stored internally
    int run();

    // Get results after run() completes
    int getOutputNx() const { return nx; }
    int getOutputNy() const { return ny; }
    int getOutputNz() const { return nz; }
    void getOutputOrigin(Ipp64f origin[3]) const {
        origin[0] = m_out_origin[0];
        origin[1] = m_out_origin[1];
        origin[2] = m_out_origin[2];
    }

    // Transfer ownership of output map (caller must free with ippsFree)
    Ipp64f* releaseOutputMap() {
        Ipp64f* map = vOutput;
        vOutput = nullptr;
        return map;
    }

    // Cleanup (called automatically by destructor)
    void cleanup();

    // Disable copy
    MapGenerator(const MapGenerator&) = delete;
    MapGenerator& operator=(const MapGenerator&) = delete;
};

// Utility functions
Ipp64f resolution_to_sigma(Ipp64f resolution, ResolutionCriterion criterion);