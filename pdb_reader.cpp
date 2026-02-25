// src/pdb_reader.cpp
#include "pdb_reader.hpp"
#include "atom.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

char PDBReader::m_error_msg[256] = { 0 };

PDBReaderConfig::PDBReaderConfig() {
    filter_hydrogen = 1;   // Filter H by default
    bfactor_cutoff = 0.0f;  // No cutoff
    use_hetatm = 1;         // Include HETATM
    verbosity = 1;
}

const char* PDBReader::get_error() {
    return m_error_msg;
}

// Check if element is hydrogen
int PDBReader::is_hydrogen(const char* element) {
    if (!element) return 0;
    // H, H1, H2, etc. but not He, Hg, etc.
    if (element[0] == 'H' || element[0] == 'h') {
        if (element[1] == '\0' || isdigit(element[1])) {
            return 1;
        }
        // Check for common hydrogen-like
        if (element[1] == 'D' || element[1] == 'd') return 1;  // Deuterium
    }
    return 0;
}

// Check if residue is water
int PDBReader::is_water(const char* resname) {
    if (!resname) return 0;
    const char* waters[] = { "HOH", "H2O", "TIP3", "TIP4", "WAT", nullptr };
    for (int i = 0; waters[i] != nullptr; i++) {
        if (strncmp(resname, waters[i], 3) == 0) return 1;
    }
    return 0;
}

// Parse element from PDB line (columns 77-78 or 13-14)
void PDBReader::parse_element(const char* line, char* element) {
    element[0] = '\0';

    // Try columns 77-78 first (standard PDB format)
    if (strlen(line) >= 78) {
        if (line[76] != ' ') {  // Column 77 (0-based index 76)
            element[0] = line[76];
            if (line[77] != ' ' && line[77] != '\0') {  // Column 78
                element[1] = line[77];
                element[2] = '\0';
            }
            else {
                element[1] = '\0';
            }
            return;
        }
    }

    // Try columns 13-14 (old PDB format)
    if (strlen(line) >= 14) {
        if (line[12] != ' ') {  // Column 13
            element[0] = line[12];
            if (line[13] != ' ' && line[13] != '\0') {  // Column 14
                element[1] = line[13];
                element[2] = '\0';
            }
            else {
                element[1] = '\0';
            }
            return;
        }
    }

    // Try to parse from atom name (some PDB files put element in atom name)
    // Atom name is in columns 13-16 (0-based 12-15)
    if (strlen(line) >= 16) {
        char atom_name[5] = { line[12], line[13], line[14], line[15], '\0' };

        // Common patterns: " C  ", " CA ", " N  ", etc.
        // Element is often the first non-space character
        for (int i = 0; i < 4; i++) {
            if (atom_name[i] != ' ') {
                element[0] = atom_name[i];
                // Check if next char is lowercase (part of element name like "Ca")
                if (i + 1 < 4 && atom_name[i + 1] >= 'a' && atom_name[i + 1] <= 'z') {
                    element[1] = atom_name[i + 1];
                    element[2] = '\0';
                }
                else {
                    element[1] = '\0';
                }
                return;
            }
        }
    }

    // Default to unknown
    strcpy(element, "X");
}

// Parse float from fixed columns
float PDBReader::parse_float(const char* line, int start, int len) {
    char buf[32];
    int i, j;

    // Skip spaces
    for (i = 0; i < len && line[start + i] == ' '; i++);

    // Copy number
    for (j = 0; i < len; i++, j++) {
        buf[j] = line[start + i];
    }
    buf[j] = '\0';

    return (float)atof(buf);
}

// First pass: count atoms and count elements
int PDBReader::first_pass(const char* buffer, const PDBReaderConfig* config,
    int* element_counts, int max_elements) {
    char line[82];
    const char* p = buffer;
    int atom_count = 0;

    // Initialize element counts to zero
    for (int i = 0; i < max_elements; i++) {
        element_counts[i] = 0;
    }

    while (*p) {
        // Copy line
        int i = 0;
        while (*p && *p != '\n' && *p != '\r' && i < 80) {
            line[i++] = *p++;
        }
        line[i] = '\0';

        // Skip newlines
        while (*p == '\n' || *p == '\r') p++;

        // Check record type
        if (strncmp(line, "ATOM  ", 6) == 0 ||
            (config->use_hetatm && strncmp(line, "HETATM", 6) == 0)) {

            char element[4];
            parse_element(line, element);

            // Apply filters
            if (config->filter_hydrogen && is_hydrogen(element)) {
                continue;
            }

            char resname[4] = { line[17], line[18], line[19], '\0' };

            float b = parse_float(line, 60, 6);
            if (config->bfactor_cutoff > 0 && b > config->bfactor_cutoff) {
                continue;
            }

            // Get atomic number
            int z = atom_get_atomic_number(element);
            if (z == 0) z = 6;  // Default to carbon if unknown

            if (z > 0 && z <= max_elements) {
                element_counts[z - 1]++;  // Increment count for this element
            }

            atom_count++;
        }
    }

    return atom_count;
}

// Second pass: read atoms with element-based indexing
int PDBReader::second_pass(const char* buffer, const PDBReaderConfig* config,
    Atom* atoms, int* element_starts, int* element_counts,
    int max_elements) {
    char line[82];
    const char* p = buffer;
    int total_atoms_read = 0;

    // Create working copy of starts that we'll increment
    int* current_pos = (int*)malloc(max_elements * sizeof(int));
    if (!current_pos) return -1;

    for (int i = 0; i < max_elements; i++) {
        current_pos[i] = element_starts[i];
    }

    while (*p) {
        int i = 0;
        while (*p && *p != '\n' && *p != '\r' && i < 80) {
            line[i++] = *p++;
        }
        line[i] = '\0';

        while (*p == '\n' || *p == '\r') p++;

        if (strncmp(line, "ATOM  ", 6) == 0 ||
            (config->use_hetatm && strncmp(line, "HETATM", 6) == 0)) {

            char element[4];
            parse_element(line, element);

            // Apply filters again
            if (config->filter_hydrogen && is_hydrogen(element)) {
                continue;
            }

            char resname[4] = { line[17], line[18], line[19], '\0' };

            float b = parse_float(line, 60, 6);
            if (config->bfactor_cutoff > 0 && b > config->bfactor_cutoff) {
                continue;
            }

            // Get atomic number
            int z = atom_get_atomic_number(element);
            if (z == 0) z = 6;  // Default to carbon if unknown

            if (z > 0 && z <= max_elements) {
                // Check if we haven't exceeded the expected count for this element
                int element_idx = z - 1;
                if (current_pos[element_idx] < element_starts[element_idx] + element_counts[element_idx]) {
                    // Place atom at current position for this element
                    int pos = current_pos[element_idx];
                    Atom* atom = &atoms[pos];

                    // Parse coordinates
                    atom->x = parse_float(line, 30, 8);
                    atom->y = parse_float(line, 38, 8);
                    atom->z = parse_float(line, 46, 8);

                    // Parse occupancy and B-factor
                    atom->occupancy = parse_float(line, 54, 6);
                    if (atom->occupancy < 0.01f) atom->occupancy = 1.0f;

                    atom->b_factor = parse_float(line, 60, 6);
                    if (atom->b_factor < 0.01f) atom->b_factor = 20.0f;  // Default

                    // Parse serial, chain, residue
                    char serial_str[6] = { line[6], line[7], line[8], line[9], line[10], '\0' };
                    atom->serial = atoi(serial_str);

                    atom->chain[0] = line[21];
                    atom->chain[1] = '\0';

                    strncpy(atom->resname, resname, 3);
                    atom->resname[3] = '\0';

                    char resseq_str[5] = { line[22], line[23], line[24], line[25], '\0' };
                    atom->resseq = atoi(resseq_str);

                    strcpy(atom->element, element);
                    atom->element_idx = -1;  // Will be set later
                    atom->hetatm = (strncmp(line, "HETATM", 6) == 0) ? 1 : 0;

                    current_pos[element_idx]++;
                    total_atoms_read++;

                    if (config->verbosity >= 2) {
                        printf("Atom %d: %s at (%.2f, %.2f, %.2f) B=%.1f (Z=%d, pos=%d)\n",
                            atom->serial, atom->element, atom->x, atom->y, atom->z, atom->b_factor, z, pos);
                    }
                }
            }
        }
    }

    free(current_pos);
    return total_atoms_read;
}

int PDBReader::read_file(const char* filename,
    const PDBReaderConfig* config,
    Atom** atoms,
    int* n_atoms) {

    FILE* fp = fopen(filename, "rb");
    if (!fp) {
        snprintf(m_error_msg, sizeof(m_error_msg),
            "Cannot open file: %s", filename);
        return -1;
    }

    // Get file size
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    // Read entire file
    char* buffer = (char*)malloc(fsize + 1);
    if (!buffer) {
        fclose(fp);
        snprintf(m_error_msg, sizeof(m_error_msg),
            "Memory allocation failed");
        return -2;
    }

    size_t bytes_read = fread(buffer, 1, fsize, fp);
    buffer[bytes_read] = '\0';
    fclose(fp);

    // First pass: count atoms and count elements
    int element_counts[118] = { 0 };  // Elements up to Oganesson
    int total_atoms = first_pass(buffer, config, element_counts, 118);

    if (total_atoms == 0) {
        free(buffer);
        snprintf(m_error_msg, sizeof(m_error_msg),
            "No atoms found after filtering");
        return -3;
    }

    if (config->verbosity >= 1) {
        printf("First pass: Found %d total atoms\n", total_atoms);
        printf("Element distribution:\n");
        const char* symbols[] = { "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
                                 "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
                                 "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
                                 "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
                                 "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac",
                                 "Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh",
                                 "Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og" };
        for (int z = 0; z < 118; z++) {
            if (element_counts[z] > 0) {
                if (z < sizeof(symbols) / sizeof(symbols[0])) {
                    printf("  %s (Z=%d): %d atoms\n", symbols[z], z + 1, element_counts[z]);
                }
                else {
                    printf("  Element Z=%d: %d atoms\n", z + 1, element_counts[z]);
                }
            }
        }
    }

    // Calculate start positions for each element
    int* element_starts = (int*)malloc(119 * sizeof(int));  // One extra for end
    if (!element_starts) {
        free(buffer);
        snprintf(m_error_msg, sizeof(m_error_msg), "Memory allocation failed");
        return -4;
    }

    element_starts[0] = 0;
    for (int z = 1; z <= 118; z++) {
        element_starts[z] = element_starts[z - 1] + element_counts[z - 1];
    }

    if (config->verbosity >= 2) {
        printf("Element start positions:\n");
        for (int z = 0; z < 118; z++) {
            if (element_counts[z] > 0) {
                printf("  Z=%d: starts at %d, count %d\n", z + 1, element_starts[z], element_counts[z]);
            }
        }
        printf("Total atoms: %d (last+1: %d)\n", total_atoms, element_starts[118]);
    }

    // Allocate atoms array
    *atoms = (Atom*)malloc(total_atoms * sizeof(Atom));
    if (!*atoms) {
        free(buffer);
        free(element_starts);
        snprintf(m_error_msg, sizeof(m_error_msg), "Memory allocation failed");
        return -4;
    }

    // Initialize the atoms array to zero
    memset(*atoms, 0, total_atoms * sizeof(Atom));

    // Second pass: read atoms
    int atoms_read = second_pass(buffer, config, *atoms, element_starts, element_counts, 118);

    free(element_starts);
    free(buffer);

    if (atoms_read != total_atoms) {
        fprintf(stderr, "Warning: Read %d atoms but expected %d\n", atoms_read, total_atoms);
        // Adjust n_atoms to actual number read
        total_atoms = atoms_read;
    }

    *n_atoms = atoms_read;

    if (config->verbosity >= 1) {
        printf("Second pass: Read %d atoms from %s\n", atoms_read, filename);
    }

    return 0;
}

void PDBReader::center_of_mass(const Atom* atoms, int n_atoms,
    Ipp64f* cx, Ipp64f* cy, Ipp64f* cz) {
    Ipp64f sx = 0.0, sy = 0.0, sz = 0.0;
    Ipp64f mass = 0.0;

    for (int i = 0; i < n_atoms; i++) {
        int z = AtomUtils::get_atomic_number(atoms[i].element);
        if (z == 0) z = 6;  // Default to carbon if unknown

        sx += z * atoms[i].x;
        sy += z * atoms[i].y;
        sz += z * atoms[i].z;
        mass += z;
    }

    if (mass > 0) {
        *cx = sx / mass;
        *cy = sy / mass;
        *cz = sz / mass;
    }
    else {
        *cx = *cy = *cz = 0.0;
    }
}

void PDBReader::geometric_center(const Atom* atoms, int n_atoms,
    Ipp64f* cx, Ipp64f* cy, Ipp64f* cz) {
    if (n_atoms == 0) {
        *cx = *cy = *cz = 0.0;
        return;
    }

    Ipp64f xmin = atoms[0].x, xmax = atoms[0].x;
    Ipp64f ymin = atoms[0].y, ymax = atoms[0].y;
    Ipp64f zmin = atoms[0].z, zmax = atoms[0].z;

    for (int i = 1; i < n_atoms; i++) {
        if (atoms[i].x < xmin) xmin = atoms[i].x;
        if (atoms[i].x > xmax) xmax = atoms[i].x;
        if (atoms[i].y < ymin) ymin = atoms[i].y;
        if (atoms[i].y > ymax) ymax = atoms[i].y;
        if (atoms[i].z < zmin) zmin = atoms[i].z;
        if (atoms[i].z > zmax) zmax = atoms[i].z;
    }

    *cx = (xmin + xmax) * 0.5;
    *cy = (ymin + ymax) * 0.5;
    *cz = (zmin + zmax) * 0.5;
}