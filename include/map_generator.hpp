// include/map_generator.hpp
#pragma once

#include "types.hpp"
#include "atom.hpp"
#include "enums.hpp"
#include "scattering.hpp"
#include "grid_projection.hpp"
#include <mkl.h>

/**
 * @file map_generator.hpp
 * @brief Default map generation using Peng1996 method with convolution chain
 *
 * Implements the complete convolution chain described in the README:
 * ρ_final = (∑ δ(r - r_n)) * K_element * G_Bn * G_res
 *
 * The method separates per-atom B-factor blurring from element-specific
 * and resolution-dependent blurring for computational efficiency.
 */

 // Map generation parameters
struct MapParams {
    Ipp64f resolution;           ///< Target resolution in Å
    ResolutionCriterion criterion; ///< Resolution to sigma conversion
    AmplitudeMode amplitude_mode; ///< Amplitude scaling mode
    Ipp64f grid_spacing;         ///< Voxel size in Å (0 = auto)
    Ipp64f padding;              ///< Extra padding around atoms in Å
    Ipp64f cutoff_level;         ///< Density cutoff for kernel (default 0.001)
    int use_bfactors;            ///< Apply B-factors from PDB
    Ipp64f b_default;            ///< Default B-factor if not in PDB

    MapParams();
};

/**
 * @struct ElementWorkspace
 * @brief Working space for processing a single element type
 *
 * Contains per-element grids and FFT buffers to avoid reallocation
 * when processing multiple elements.
 */
struct ElementWorkspace {
    Ipp64f* grid;                 ///< Real-space grid for this element
    Ipp64fc* fft_buffer;          ///< FFT buffer for this element's grid
    Ipp64fc* kernel_buffer;       ///< FFT buffer for this element's kernel
    DFTI_DESCRIPTOR_HANDLE fft_desc;  ///< FFT descriptor (optional per-element)

    ElementWorkspace();
    ~ElementWorkspace();
    void cleanup();
};

/**
 * @class MapGenerator
 * @brief Implements the default Peng1996 method with convolution chain
 *
 * The generation follows a two-step approach:
 * 1. B-factor blurring: Each atom is distributed to nearby grid points using
 *    analytical integration (error functions) for exact occupancy conservation
 * 2. Combined element+resolution kernel applied via FFT convolution
 *
 * Proper normalization ensures each atom's total contribution equals fe0(0)
 * and B-factor handling follows the formulas in the README.
 */
class MapGenerator {
private:
    // Scattering data
    ScatteringTable m_table;
    ElementProfile* m_profiles;      ///< Profiles for elements present
    int nProfiles;                    ///< Number of profiles
    int* m_element_to_profile;        ///< Maps element_idx -> profile index

    // Atom data
    const Atom* m_atoms;              ///< Pointer to atoms (not owned)
    int nAtoms;                        ///< Number of atoms
    int* vStartAtom;                   ///< Array of unique element indices
    int nUnique;                       ///< Number of unique elements

    // Working spaces
    ElementWorkspace* m_elem_workspace;
    int m_n_elem_workspace;

    // Grid parameters
    Ipp64f m_spacing;                  ///< Grid spacing in Å
    int nx, ny, nz;                     ///< Local grid dimensions
    int m_xyz_offset;                   ///< Offset for padding
    Ipp64f m_origin[3];                 ///< Grid origin

    // Buffer sizes
    int64_t nt;                         ///< Size of real FFT buffer (nx*ny*nz)
    int64_t nct;                        ///< Size of complex FFT buffer (nz*ny*(nx/2+1))

    // Output data
    Ipp64f m_out_origin[3];              ///< Output origin for final map

    // FFT buffers
    Ipp64f* vData, * vKernel, * vOutput;
    Ipp64fc* vcData, * vcKernel;
    DFTI_DESCRIPTOR_HANDLE descFor, descBack;

    // Parameters
    MapParams m_params;

    /**
     * @brief Find unique elements present in atoms
     * @return 0 on success, negative on error
     */
    int findUniqueElements();

    /**
     * @brief Calculate grid dimensions based on atom positions and profiles
     * @return 0 on success, negative on error
     */
    int calculateGridDimensions();

    /**
     * @brief Allocate FFT buffers
     * @return 0 on success, negative on error
     */
    int allocateFFTBuffers();

    /**
     * @brief Setup MKL FFT descriptors
     * @return 0 on success, negative on error
     */
    int setupFFTDescriptor();

    /**
     * @brief Process a single element type
     * @param element_idx Index of element to process
     * @return 0 on success, negative on error
     */
    int processElement(int element_idx);

    /**
     * @brief Project atoms of a specific element with B-factor blurring
     *
     * Implements Step 1 from README: delta function + B-factor blurring combined.
     * Uses analytical integration (error functions) to distribute atoms to grid.
     *
     * @param element_idx Element index to project
     * @param grid Output grid for this element
     * @param xyz_offset Padding offset
     */
    void projectElementAtoms(int element_idx, Ipp64f* grid, int xyz_offset);

    /**
     * @brief Create combined kernel for element + resolution blur
     *
     * Implements combined kernel from README:
     * K_combined = K_element * G_res
     *
     * @param element_idx Element index
     * @param kernel Output kernel array
     * @param kernel_radius_voxels Kernel radius in voxels
     */
    void createElementKernel(int element_idx, Ipp64f* kernel, int kernel_radius_voxels);

    /**
     * @brief Add processed element map to final output
     * @param src Source element map
     */
    void addToFinalMap(const Ipp64f* src);

    /**
     * @brief Normalize final map
     *
     * Ensures proper scaling according to amplitude mode:
     * - Peng1996: Each atom's total = fe0(0)
     * - Atomic Number: Each atom's total = Z
     */
    void normalizeMap();

    /**
     * @brief Print map statistics for debugging
     */
    void printMapStatistics();

public:
    MapGenerator();
    ~MapGenerator();

    /**
     * @brief Initialize generator with parameters and atoms
     * @param params Map generation parameters
     * @param atoms Array of atoms
     * @param n_atoms Number of atoms
     * @return 0 on success, negative on error
     */
    int init(const MapParams* params, const Atom* atoms, int n_atoms);

    /**
     * @brief Run map generation
     * @return 0 on success, negative on error
     */
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

    /**
     * @brief Transfer ownership of output map
     * @return Pointer to output map (caller must free with ippsFree)
     */
    Ipp64f* releaseOutputMap() {
        Ipp64f* map = vOutput;
        vOutput = nullptr;
        return map;
    }

    /**
     * @brief Clean up all allocated resources
     */
    void cleanup();

    // Disable copy
    MapGenerator(const MapGenerator&) = delete;
    MapGenerator& operator=(const MapGenerator&) = delete;
};

/**
 * @brief Convert resolution to Gaussian sigma based on criterion
 * @param resolution Target resolution in Å
 * @param criterion Resolution criterion to use
 * @return Gaussian sigma in Å
 */
Ipp64f resolution_to_sigma(Ipp64f resolution, ResolutionCriterion criterion);