// include/chimerax_generator.hpp
#pragma once

#include "types.hpp"
#include "atom.hpp"
#include <vector>

class ChimeraXGenerator {
private:
    std::vector<float> m_centers;  // i,j,k coordinates in grid units
    std::vector<float> m_weights;   // atomic numbers (Z)
    std::vector<float> m_sdevs;     // normalized sigmas (σ/step) for each dimension
    float m_cutoff_range;            // typically 5.0
    float m_sigma;                   // Gaussian width in Å

    Ipp64f* m_output;
    int nx, ny, nz;
    Ipp64f m_origin[3];
    Ipp64f m_step;                    // grid spacing

    // Transform atom coordinates to grid indices
    void xyz_to_ijk(const Atom* atoms, int n_atoms,
        const Ipp64f origin[3], Ipp64f step,
        std::vector<float>& ijk);

public:
    ChimeraXGenerator();
    ~ChimeraXGenerator();

    int init(const Atom* atoms, int n_atoms,
        Ipp64f resolution, Ipp64f grid_spacing,
        Ipp64f padding, float cutoff_range = 5.0f);

    int run();

    Ipp64f* releaseOutput() {
        Ipp64f* map = m_output;
        m_output = nullptr;
        return map;
    }

    void getDimensions(int& x, int& y, int& z) const { x = nx; y = ny; z = nz; }
    void getOrigin(Ipp64f origin[3]) const {
        origin[0] = m_origin[0];
        origin[1] = m_origin[1];
        origin[2] = m_origin[2];
    }
    Ipp64f getStep() const { return m_step; }

    ChimeraXGenerator(const ChimeraXGenerator&) = delete;
    ChimeraXGenerator& operator=(const ChimeraXGenerator&) = delete;
};