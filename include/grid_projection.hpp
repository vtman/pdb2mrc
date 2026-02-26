// include/grid_projection.hpp
#pragma once

#include "types.hpp"
#include "atom.hpp"
#include "scattering.hpp"

/**
 * Projection grid with double precision
 */
struct ProjectionGrid {
    int nx, ny, nz;           // Grid dimensions
    Ipp64f sx, sy, sz;        // Grid spacing (Å)
    Ipp64f ox, oy, oz;        // Origin (min coordinates)
    int64_t nvox;             // Total voxels (nx*ny*nz)
    Ipp64f* grid;             // Projected grid (weights)
    Ipp64f* work;             // Working buffer

    ProjectionGrid();
    ~ProjectionGrid();
    void cleanup();
};

class GridProjection {
private:
    static char m_error_msg[256];

public:
    /**
     * Initialize projection grid
     * @param grid Grid structure
     * @param bbox Bounding box
     * @param spacing Grid spacing in Å
     * @param pad_fft If 1, pad to powers of 2 for FFT
     * @return 0 on success, negative on error
     */
    static int init(ProjectionGrid* grid,
        const BoundingBox* bbox,
        Ipp64f spacing,
        int pad_fft);

    /**
     * Free grid memory
     * @param grid Grid structure
     */
    static void free_grid(ProjectionGrid* grid);

    /**
     * Clear grid (set to zero)
     * @param grid Grid structure
     */
    static void clear(ProjectionGrid* grid);

    /**
     * Calculate bounding box from atoms with padding
     * @param atoms Array of atoms
     * @param n_atoms Number of atoms
     * @param padding Padding to add (Å)
     * @param bbox Output bounding box
     */
    static void calc_bounding_box(const Atom* atoms,
        int n_atoms,
        Ipp64f padding,
        BoundingBox* bbox);

    /**
     * Project atoms to grid using trilinear interpolation
     * Each atom contributes to 8 surrounding voxels with weights
     *
     * @param grid Projection grid
     * @param atoms Array of atoms
     * @param n_atoms Number of atoms
     * @param table Scattering table (for fe0 values)
     * @return 0 on success, negative on error
     */
    static int project_atoms(ProjectionGrid* grid,
        const Atom* atoms,
        int n_atoms,
        const ScatteringTable* table);

    /**
     * Get error message from last operation
     * @return Error string
     */
    static const char* get_error();
};

// Inline helpers
static inline Ipp64f grid_get(const ProjectionGrid* grid,
    int i, int j, int k) {
    return grid->grid[((int64_t)k * grid->ny + j) * grid->nx + i];
}

static inline void grid_set(ProjectionGrid* grid,
    int i, int j, int k,
    Ipp64f value) {
    grid->grid[((int64_t)k * grid->ny + j) * grid->nx + i] = value;
}

static inline void grid_add(ProjectionGrid* grid,
    int i, int j, int k,
    Ipp64f value) {
    grid->grid[((int64_t)k * grid->ny + j) * grid->nx + i] += value;
}

// For backward compatibility
inline int grid_init(ProjectionGrid* grid,
    const BoundingBox* bbox,
    Ipp64f spacing,
    int pad_fft) {
    return GridProjection::init(grid, bbox, spacing, pad_fft);
}

inline void grid_free(ProjectionGrid* grid) {
    GridProjection::free_grid(grid);
}

inline void grid_clear(ProjectionGrid* grid) {
    GridProjection::clear(grid);
}

inline void grid_calc_bounding_box(const Atom* atoms,
    int n_atoms,
    Ipp64f padding,
    BoundingBox* bbox) {
    GridProjection::calc_bounding_box(atoms, n_atoms, padding, bbox);
}

inline int grid_project_atoms(ProjectionGrid* grid,
    const Atom* atoms,
    int n_atoms,
    const ScatteringTable* table) {
    return GridProjection::project_atoms(grid, atoms, n_atoms, table);
}