// include/pdb_reader.hpp
#pragma once

#include "atom.hpp"
#include <stdio.h>

// PDB reader configuration
struct PDBReaderConfig {
    int filter_hydrogen;      // Exclude hydrogen atoms
    float bfactor_cutoff;     // Exclude atoms with B > cutoff (0 = no cutoff)
    int use_hetatm;           // Include HETATM records
    int verbosity;            // 0 = silent, 1 = basic, 2 = verbose

    PDBReaderConfig();
};

class PDBReader {
private:
    static char m_error_msg[256];

    // Helper methods
    static int is_hydrogen(const char* element);
    static int is_water(const char* resname);
    static void parse_element(const char* line, char* element);
    static float parse_float(const char* line, int start, int len);
    static int first_pass(const char* buffer, const PDBReaderConfig* config,
        int* element_counts, int max_elements);
    static int second_pass(const char* buffer, const PDBReaderConfig* config,
        Atom* atoms, int* element_starts, int* element_counts,
        int max_elements);

public:
    static int read_file(const char* filename,
        const PDBReaderConfig* config,
        Atom** atoms,
        int* n_atoms);

    static const char* get_error();

    // Use Ipp64f (double) for coordinates - no unnecessary conversions
    static void center_of_mass(const Atom* atoms, int n_atoms,
        Ipp64f* cx, Ipp64f* cy, Ipp64f* cz);

    static void geometric_center(const Atom* atoms, int n_atoms,
        Ipp64f* cx, Ipp64f* cy, Ipp64f* cz);
};