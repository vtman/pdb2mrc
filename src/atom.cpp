// src/atom.cpp
#include "atom.hpp"
#include <string.h>
#include <ctype.h>
#include <math.h>

// Atomic numbers for common elements
const AtomUtils::AtomicEntry AtomUtils::ATOMIC_NUMBERS[] = {
    {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6},
    {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12},
    {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18},
    {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
    {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
    {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36},
    {"Rb", 37}, {"Sr", 38}, {"Y", 39}, {"Zr", 40}, {"Nb", 41}, {"Mo", 42},
    {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48},
    {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53}, {"Xe", 54},
    {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60},
    {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66},
    {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72},
    {"Ta", 73}, {"W", 74}, {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78},
    {"Au", 79}, {"Hg", 80}, {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84},
    {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
    {"Pa", 91}, {"U", 92}, {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96},
    {"Bk", 97}, {"Cf", 98}, {nullptr, 0}
};

// Case-insensitive string comparison
int AtomUtils::strieq(const char* a, const char* b) {
    while (*a && *b) {
        if (toupper(*a) != toupper(*b)) return 0;
        a++; b++;
    }
    return *a == *b;
}

int AtomUtils::get_atomic_number(const char* element) {
    if (!element) return 0;

    for (int i = 0; ATOMIC_NUMBERS[i].symbol != nullptr; i++) {
        if (strieq(ATOMIC_NUMBERS[i].symbol, element)) {
            return ATOMIC_NUMBERS[i].z;
        }
    }
    return 0;  // Unknown
}

void AtomUtils::calc_stats(const Atom* atoms, int n_atoms, AtomStats* stats) {
    if (!atoms || n_atoms == 0 || !stats) return;

    // Initialize
    stats->center_of_mass[0] = 0.0;
    stats->center_of_mass[1] = 0.0;
    stats->center_of_mass[2] = 0.0;

    stats->geometric_center[0] = 0.0;
    stats->geometric_center[1] = 0.0;
    stats->geometric_center[2] = 0.0;

    stats->radius = 0.0;
    stats->n_atoms = n_atoms;

    // Find bounds for geometric center
    Ipp64f xmin = atoms[0].x, xmax = atoms[0].x;
    Ipp64f ymin = atoms[0].y, ymax = atoms[0].y;
    Ipp64f zmin = atoms[0].z, zmax = atoms[0].z;

    // Calculate center of mass
    Ipp64f mass_sum = 0.0;
    Ipp64f com_x = 0.0, com_y = 0.0, com_z = 0.0;

    for (int i = 0; i < n_atoms; i++) {
        int z = get_atomic_number(atoms[i].element);
        if (z == 0) z = 6;  // Default to carbon

        com_x += z * atoms[i].x;
        com_y += z * atoms[i].y;
        com_z += z * atoms[i].z;
        mass_sum += z;

        if (atoms[i].x < xmin) xmin = atoms[i].x;
        if (atoms[i].x > xmax) xmax = atoms[i].x;
        if (atoms[i].y < ymin) ymin = atoms[i].y;
        if (atoms[i].y > ymax) ymax = atoms[i].y;
        if (atoms[i].z < zmin) zmin = atoms[i].z;
        if (atoms[i].z > zmax) zmax = atoms[i].z;
    }

    if (mass_sum > 0) {
        stats->center_of_mass[0] = com_x / mass_sum;
        stats->center_of_mass[1] = com_y / mass_sum;
        stats->center_of_mass[2] = com_z / mass_sum;
    }

    stats->geometric_center[0] = (xmin + xmax) * 0.5;
    stats->geometric_center[1] = (ymin + ymax) * 0.5;
    stats->geometric_center[2] = (zmin + zmax) * 0.5;

    // Calculate radius (max distance from geometric center)
    for (int i = 0; i < n_atoms; i++) {
        Ipp64f dx = atoms[i].x - stats->geometric_center[0];
        Ipp64f dy = atoms[i].y - stats->geometric_center[1];
        Ipp64f dz = atoms[i].z - stats->geometric_center[2];
        Ipp64f dist = sqrt(dx * dx + dy * dy + dz * dz);
        if (dist > stats->radius) stats->radius = dist;
    }

    // Count unique elements (simplified)
    stats->n_elements = 0;
    // This would require a more complex implementation to count unique elements
}