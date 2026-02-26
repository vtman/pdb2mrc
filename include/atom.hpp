// include/atom.hpp
#pragma once

#include "types.hpp"

/**
 * Atom structure with double precision
 */
struct Atom {
    Ipp64f x, y, z;           // Coordinates in Å
    Ipp64f occupancy;         // Occupancy (0-1)
    Ipp64f b_factor;          // Isotropic B-factor (Å²)
    int element_idx;          // Index into scattering table
    int serial;               // PDB serial number
    char element[4];          // Element symbol (e.g., "C", "FE")
    char chain[2];            // Chain identifier
    char resname[5];          // Residue name
    int resseq;               // Residue sequence number
    int hetatm;               // 1 for HETATM, 0 for ATOM
};

/**
 * Bounding box with double precision
 */
struct BoundingBox {
    Ipp64f xmin, xmax;
    Ipp64f ymin, ymax;
    Ipp64f zmin, zmax;
    Ipp64f padding;
};

/**
 * Atom statistics
 */
struct AtomStats {
    Ipp64f center_of_mass[3];
    Ipp64f geometric_center[3];
    Ipp64f radius;            // Maximum distance from center
    int n_atoms;
    int n_elements;           // Number of unique elements
};

class AtomUtils {
private:
    // Atomic numbers for common elements
    struct AtomicEntry {
        const char* symbol;
        int z;
    };

    static const AtomicEntry ATOMIC_NUMBERS[];
    static int strieq(const char* a, const char* b);

public:
    /**
     * Calculate atom statistics
     * @param atoms Array of atoms
     * @param n_atoms Number of atoms
     * @param stats Output statistics
     */
    static void calc_stats(const Atom* atoms, int n_atoms, AtomStats* stats);

    /**
     * Get atomic number from element symbol
     * @param element Element symbol
     * @return Atomic number or 0 if unknown
     */
    static int get_atomic_number(const char* element);
};

// For backward compatibility with existing code
inline void atom_calc_stats(const Atom* atoms, int n_atoms, AtomStats* stats) {
    AtomUtils::calc_stats(atoms, n_atoms, stats);
}

inline int atom_get_atomic_number(const char* element) {
    return AtomUtils::get_atomic_number(element);
}