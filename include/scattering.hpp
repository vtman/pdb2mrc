// include/scattering.hpp
#pragma once

#include "types.hpp"
#include "atom.hpp"
#include "enums.hpp"
#include <math.h>

#define MAX_ELEMENTS 98
#define N_GAUSS 5
#define PROFILE_POINTS 1000  // Number of points in radial profile

// Peng 1996 5-Gaussian coefficients
struct ScatteringEntry {
    char name[8];             // Element symbol
    Ipp64f a[N_GAUSS];        // Amplitudes
    Ipp64f b[N_GAUSS];        // Falloffs (A^2)
    Ipp64f fe0;               // f_e(0) = sum(a)
    int atomic_number;        // Atomic number Z
};

// Scattering factor table
struct ScatteringTable {
    ScatteringEntry entries[MAX_ELEMENTS];
    int n_entries;
};

// Precomputed radial profile for an element at given resolution
struct ElementProfile {
    int element_idx;          // Index in scattering table
    Ipp64f* profile;          // Radial profile values (1D) - size PROFILE_POINTS
    int profile_size;         // Number of points in profile (should be PROFILE_POINTS)
    Ipp64f step;              // Step size between profile points
    Ipp64f rmax;              // Maximum radius
    Ipp64f spacing;           // Grid spacing used (for reference)
};

class Scattering {
private:
    static const ScatteringEntry PENG1996_NEUTRAL[];
    static int strieq(const char* a, const char* b);
    static Ipp64f atom_density_at_r(const ScatteringEntry* entry,
        Ipp64f r, Ipp64f resolution, ResolutionCriterion criterion);

public:
    // Initialize Peng 1996 table with neutral atoms
    static void init(ScatteringTable* table);

    // Get element index (case-insensitive)
    static int get_index(const ScatteringTable* table, const char* element);

    // Precompute profiles for elements present in atoms at given resolution
    static int precompute_profiles(const ScatteringTable* table,
        Ipp64f resolution,
        ResolutionCriterion criterion,
        Ipp64f grid_spacing,
        Ipp64f cutoff_level,
        AmplitudeMode amplitude_mode,
        const Atom* atoms,
        int n_atoms,
        ElementProfile** profiles,
        int* n_profiles);

    // Free precomputed profiles
    static void free_profiles(ElementProfile* profiles, int n_profiles);

    // Get interpolated profile value at distance r
    static inline Ipp64f profile_at_r(const ElementProfile* profile, Ipp64f r) {
        if (!profile || !profile->profile) return 0.0;

        if (r <= 0.0) return profile->profile[0];
        if (r >= profile->rmax) return profile->profile[profile->profile_size - 1];

        Ipp64f idx = r / profile->step;
        int i0 = (int)idx;
        int i1 = i0 + 1;

        if (i1 >= profile->profile_size) {
            return profile->profile[profile->profile_size - 1];
        }

        Ipp64f frac = idx - i0;
        return profile->profile[i0] * (1.0 - frac) + profile->profile[i1] * frac;
    }
};

// For backward compatibility with existing code
inline void scattering_init(ScatteringTable* table) {
    Scattering::init(table);
}

inline int scattering_get_index(const ScatteringTable* table, const char* element) {
    return Scattering::get_index(table, element);
}

// UPDATED: Backward compatibility wrapper with default criterion
inline int scattering_precompute_profiles(const ScatteringTable* table,
    Ipp64f resolution,
    Ipp64f grid_spacing,
    Ipp64f cutoff_level,
    AmplitudeMode amplitude_mode,
    const Atom* atoms,
    int n_atoms,
    ElementProfile** profiles,
    int* n_profiles) {

    // Use Rayleigh criterion as default for backward compatibility
    return Scattering::precompute_profiles(table, resolution, CRITERION_RAYLEIGH,
        grid_spacing, cutoff_level, amplitude_mode,
        atoms, n_atoms, profiles, n_profiles);
}

inline Ipp64f scattering_profile_at_r(const ElementProfile* profile, Ipp64f r) {
    return Scattering::profile_at_r(profile, r);
}