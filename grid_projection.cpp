// src/grid_projection.cpp
#include "grid_projection.hpp"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

char GridProjection::m_error_msg[256] = { 0 };

ProjectionGrid::ProjectionGrid() {
    nx = ny = nz = 0;
    sx = sy = sz = 0.0;
    ox = oy = oz = 0.0;
    nvox = 0;
    grid = nullptr;
    work = nullptr;
}

ProjectionGrid::~ProjectionGrid() {
    cleanup();
}

void ProjectionGrid::cleanup() {
    if (grid) {
        ippsFree(grid);
        grid = nullptr;
    }
    if (work) {
        ippsFree(work);
        work = nullptr;
    }
}

const char* GridProjection::get_error() {
    return m_error_msg;
}

int GridProjection::init(ProjectionGrid* grid,
    const BoundingBox* bbox,
    Ipp64f spacing,
    int pad_fft) {

    if (!grid || !bbox) return -1;

    grid->sx = spacing;
    grid->sy = spacing;
    grid->sz = spacing;
    grid->ox = bbox->xmin;
    grid->oy = bbox->ymin;
    grid->oz = bbox->zmin;

    // Calculate base dimensions
    int nx_base = (int)((bbox->xmax - bbox->xmin) / spacing) + 1;
    int ny_base = (int)((bbox->ymax - bbox->ymin) / spacing) + 1;
    int nz_base = (int)((bbox->zmax - bbox->zmin) / spacing) + 1;

    if (pad_fft) {
        // Pad to powers of 2 for FFT efficiency
        grid->nx = 1;
        while (grid->nx < nx_base) grid->nx <<= 1;

        grid->ny = 1;
        while (grid->ny < ny_base) grid->ny <<= 1;

        grid->nz = 1;
        while (grid->nz < nz_base) grid->nz <<= 1;
    }
    else {
        grid->nx = nx_base;
        grid->ny = ny_base;
        grid->nz = nz_base;
    }

    grid->nvox = (int64_t)grid->nx * grid->ny * grid->nz;

    // Allocate grid
    grid->grid = (Ipp64f*)ippsMalloc_64f(grid->nvox);
    if (!grid->grid) return -2;

    grid->work = (Ipp64f*)ippsMalloc_64f(grid->nvox);
    if (!grid->work) {
        ippsFree(grid->grid);
        grid->grid = nullptr;
        return -3;
    }

    ippsZero_64f(grid->grid, grid->nvox);
    ippsZero_64f(grid->work, grid->nvox);

    return 0;
}

void GridProjection::free_grid(ProjectionGrid* grid) {
    if (!grid) return;
    grid->cleanup();
}

void GridProjection::clear(ProjectionGrid* grid) {
    if (grid && grid->grid) {
        ippsZero_64f(grid->grid, grid->nvox);
    }
}

void GridProjection::calc_bounding_box(const Atom* atoms,
    int n_atoms,
    Ipp64f padding,
    BoundingBox* bbox) {

    if (n_atoms == 0) {
        bbox->xmin = bbox->xmax = 0.0;
        bbox->ymin = bbox->ymax = 0.0;
        bbox->zmin = bbox->zmax = 0.0;
        bbox->padding = padding;
        return;
    }

    bbox->xmin = bbox->xmax = atoms[0].x;
    bbox->ymin = bbox->ymax = atoms[0].y;
    bbox->zmin = bbox->zmax = atoms[0].z;

    for (int i = 1; i < n_atoms; i++) {
        if (atoms[i].x < bbox->xmin) bbox->xmin = atoms[i].x;
        if (atoms[i].x > bbox->xmax) bbox->xmax = atoms[i].x;
        if (atoms[i].y < bbox->ymin) bbox->ymin = atoms[i].y;
        if (atoms[i].y > bbox->ymax) bbox->ymax = atoms[i].y;
        if (atoms[i].z < bbox->zmin) bbox->zmin = atoms[i].z;
        if (atoms[i].z > bbox->zmax) bbox->zmax = atoms[i].z;
    }

    // Add padding
    bbox->xmin -= padding;
    bbox->xmax += padding;
    bbox->ymin -= padding;
    bbox->ymax += padding;
    bbox->zmin -= padding;
    bbox->zmax += padding;

    bbox->padding = padding;
}

int GridProjection::project_atoms(ProjectionGrid* grid,
    const Atom* atoms,
    int n_atoms,
    const ScatteringTable* table) {

    if (!grid || !grid->grid || !atoms || n_atoms <= 0) return -1;

    clear(grid);

    Ipp64f inv_sx = 1.0 / grid->sx;
    Ipp64f inv_sy = 1.0 / grid->sy;
    Ipp64f inv_sz = 1.0 / grid->sz;

    // For each atom
    for (int i = 0; i < n_atoms; i++) {
        const Atom* atom = &atoms[i];

        // Get fe0 (scattering power at zero angle)
        Ipp64f fe0;
        if (table && atom->element_idx >= 0 &&
            atom->element_idx < table->n_entries) {
            fe0 = table->entries[atom->element_idx].fe0;
        }
        else {
            // Default to carbon if unknown
            fe0 = 2.808;
        }

        Ipp64f weight = fe0 * atom->occupancy;

        // Position in grid coordinates
        Ipp64f gx = (atom->x - grid->ox) * inv_sx;
        Ipp64f gy = (atom->y - grid->oy) * inv_sy;
        Ipp64f gz = (atom->z - grid->oz) * inv_sz;

        // Integer coordinates (floor)
        int ix0 = (int)gx;
        int iy0 = (int)gy;
        int iz0 = (int)gz;

        // Check bounds
        if (ix0 < 0 || ix0 >= grid->nx - 1 ||
            iy0 < 0 || iy0 >= grid->ny - 1 ||
            iz0 < 0 || iz0 >= grid->nz - 1) {
            // Atom outside grid - shouldn't happen with padding
            continue;
        }

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

        // Add to 8 surrounding voxels
        int64_t base_idx = ((int64_t)iz0 * grid->ny + iy0) * grid->nx + ix0;

        grid->grid[base_idx] += weight * w000;
        grid->grid[base_idx + 1] += weight * w100;

        int64_t next_row = base_idx + grid->nx;
        grid->grid[next_row] += weight * w010;
        grid->grid[next_row + 1] += weight * w110;

        int64_t next_slice = base_idx + (int64_t)grid->nx * grid->ny;
        grid->grid[next_slice] += weight * w001;
        grid->grid[next_slice + 1] += weight * w101;

        int64_t next_slice_row = next_slice + grid->nx;
        grid->grid[next_slice_row] += weight * w011;
        grid->grid[next_slice_row + 1] += weight * w111;
    }

    return 0;
}