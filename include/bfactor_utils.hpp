// include/bfactor_utils.hpp
#pragma once

#include "atom.hpp"
#include "enums.hpp"

#define MAX_UNIQUE_ELEMENTS 50

struct ElementBStats {
    char element[4];           // Element symbol
    int count;                  // Number of atoms
    double mean;                 // Mean B-factor
    double median;               // Median B-factor
    double rms;                  // Root mean square
    double min;                  // Minimum
    double max;                  // Maximum
    double eq_resolution;        // Equivalent resolution
};

struct GlobalBStats {
    int count;                   // Total atoms
    double mean;                 // Mean B-factor
    double median;               // Median B-factor
    double rms;                  // Root mean square
    double min;                  // Minimum
    double max;                  // Maximum
    double eq_res_mean;          // Equivalent resolution from mean
    double eq_res_rms;           // Equivalent resolution from RMS
};

class BFactorAnalyzer {
private:
    static int strieq(const char* a, const char* b);
    static int find_element_index(ElementBStats* stats, int n_stats, const char* element);
    static int compare_double(const void* a, const void* b);

public:
    /**
     * Calculate per-element B-factor statistics
     * @param atoms Array of atoms
     * @param n_atoms Number of atoms
     * @param elem_stats Output array for element statistics (caller must free)
     * @param n_elem Output number of unique elements
     * @return 0 on success, negative on error
     */
    static int calculate_per_element_stats(const Atom* atoms, int n_atoms,
        ElementBStats** elem_stats, int* n_elem);

    /**
     * Calculate global B-factor statistics
     * @param atoms Array of atoms
     * @param n_atoms Number of atoms
     * @param stats Output global statistics
     * @return 0 on success
     */
    static void calculate_global_stats(const Atom* atoms, int n_atoms, GlobalBStats* stats);

    /**
     * Convert B-factor to equivalent resolution based on criterion
     */
    static double bfactor_to_resolution(double B, ResolutionCriterion criterion);

    /**
     * Print statistics to console
     */
    static void print_stats(const ElementBStats* elem_stats, int n_elem,
        const GlobalBStats* global_stats,
        ResolutionCriterion criterion,
        MapGenerationMethod method);

    /**
     * Free allocated element statistics
     */
    static void free_stats(ElementBStats* elem_stats);
};