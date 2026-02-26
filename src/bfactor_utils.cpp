// src/bfactor_utils.cpp
#include "bfactor_utils.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

// Case-insensitive string comparison
int BFactorAnalyzer::strieq(const char* a, const char* b) {
    while (*a && *b) {
        if (toupper(*a) != toupper(*b)) return 0;
        a++; b++;
    }
    return *a == *b;
}

// Find element index in stats array
int BFactorAnalyzer::find_element_index(ElementBStats* stats, int n_stats, const char* element) {
    for (int i = 0; i < n_stats; i++) {
        if (strieq(stats[i].element, element)) {
            return i;
        }
    }
    return -1;
}

// Comparison function for qsort (doubles)
int BFactorAnalyzer::compare_double(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

int BFactorAnalyzer::calculate_per_element_stats(const Atom* atoms, int n_atoms,
    ElementBStats** elem_stats, int* n_elem) {
    if (!atoms || n_atoms <= 0 || !elem_stats || !n_elem) return -1;

    // First pass: count unique elements
    ElementBStats temp_stats[MAX_UNIQUE_ELEMENTS];
    int n_unique = 0;

    for (int i = 0; i < n_atoms; i++) {
        const char* element = atoms[i].element;
        int idx = find_element_index(temp_stats, n_unique, element);

        if (idx < 0 && n_unique < MAX_UNIQUE_ELEMENTS) {
            // New element
            strcpy(temp_stats[n_unique].element, element);
            temp_stats[n_unique].count = 1;
            n_unique++;
        }
        else if (idx >= 0) {
            // Existing element
            temp_stats[idx].count++;
        }
    }

    if (n_unique == 0) {
        *elem_stats = NULL;
        *n_elem = 0;
        return 0;
    }

    // Allocate final stats array
    *elem_stats = (ElementBStats*)malloc(n_unique * sizeof(ElementBStats));
    if (!*elem_stats) return -2;

    // Initialize
    for (int i = 0; i < n_unique; i++) {
        strcpy((*elem_stats)[i].element, temp_stats[i].element);
        (*elem_stats)[i].count = temp_stats[i].count;
        (*elem_stats)[i].mean = 0.0;
        (*elem_stats)[i].median = 0.0;
        (*elem_stats)[i].rms = 0.0;
        (*elem_stats)[i].min = 1e10;
        (*elem_stats)[i].max = -1e10;
    }

    // Allocate temporary arrays for each element's B-factors
    double** b_arrays = (double**)malloc(n_unique * sizeof(double*));
    if (!b_arrays) {
        free(*elem_stats);
        return -3;
    }

    for (int i = 0; i < n_unique; i++) {
        b_arrays[i] = (double*)malloc((*elem_stats)[i].count * sizeof(double));
        if (!b_arrays[i]) {
            // Clean up already allocated arrays
            for (int j = 0; j < i; j++) free(b_arrays[j]);
            free(b_arrays);
            free(*elem_stats);
            return -4;
        }
    }

    // Reset counts for filling arrays
    int* current_pos = (int*)calloc(n_unique, sizeof(int));
    if (!current_pos) {
        for (int i = 0; i < n_unique; i++) free(b_arrays[i]);
        free(b_arrays);
        free(*elem_stats);
        return -5;
    }

    // Second pass: fill B-factor arrays
    for (int i = 0; i < n_atoms; i++) {
        int idx = find_element_index(*elem_stats, n_unique, atoms[i].element);
        if (idx >= 0) {
            int pos = current_pos[idx];
            if (pos < (*elem_stats)[idx].count) {
                b_arrays[idx][pos] = atoms[i].b_factor;
                current_pos[idx]++;

                // Update min/max
                if (atoms[i].b_factor < (*elem_stats)[idx].min) {
                    (*elem_stats)[idx].min = atoms[i].b_factor;
                }
                if (atoms[i].b_factor > (*elem_stats)[idx].max) {
                    (*elem_stats)[idx].max = atoms[i].b_factor;
                }
            }
        }
    }

    // Calculate statistics for each element
    for (int i = 0; i < n_unique; i++) {
        double sum = 0.0;
        double sum_sq = 0.0;

        for (int j = 0; j < (*elem_stats)[i].count; j++) {
            double b = b_arrays[i][j];
            sum += b;
            sum_sq += b * b;
        }

        (*elem_stats)[i].mean = sum / (*elem_stats)[i].count;
        (*elem_stats)[i].rms = sqrt(sum_sq / (*elem_stats)[i].count);

        // Calculate median (sort the array)
        qsort(b_arrays[i], (*elem_stats)[i].count, sizeof(double), compare_double);

        int mid = (*elem_stats)[i].count / 2;
        if ((*elem_stats)[i].count % 2 == 1) {
            (*elem_stats)[i].median = b_arrays[i][mid];
        }
        else {
            (*elem_stats)[i].median = (b_arrays[i][mid - 1] + b_arrays[i][mid]) / 2.0;
        }
    }

    // Clean up temporary arrays
    for (int i = 0; i < n_unique; i++) free(b_arrays[i]);
    free(b_arrays);
    free(current_pos);

    *n_elem = n_unique;
    return 0;
}

void BFactorAnalyzer::calculate_global_stats(const Atom* atoms, int n_atoms, GlobalBStats* stats) {
    if (!atoms || n_atoms <= 0 || !stats) return;

    stats->count = n_atoms;
    stats->min = 1e10;
    stats->max = -1e10;

    double sum = 0.0;
    double sum_sq = 0.0;

    for (int i = 0; i < n_atoms; i++) {
        double b = atoms[i].b_factor;
        sum += b;
        sum_sq += b * b;

        if (b < stats->min) stats->min = b;
        if (b > stats->max) stats->max = b;
    }

    stats->mean = sum / n_atoms;
    stats->rms = sqrt(sum_sq / n_atoms);

    // Create array for median calculation
    double* values = (double*)malloc(n_atoms * sizeof(double));
    if (values) {
        for (int i = 0; i < n_atoms; i++) values[i] = atoms[i].b_factor;
        qsort(values, n_atoms, sizeof(double), compare_double);

        int mid = n_atoms / 2;
        if (n_atoms % 2 == 1) {
            stats->median = values[mid];
        }
        else {
            stats->median = (values[mid - 1] + values[mid]) / 2.0;
        }
        free(values);
    }
    else {
        stats->median = stats->mean;  // Fallback
    }
}

double BFactorAnalyzer::bfactor_to_resolution(double B, ResolutionCriterion criterion) {
    // σ_B² = B/(8π²)
    double sigma_B = sqrt(B / (8.0 * M_PI * M_PI));

    switch (criterion) {
    case CRITERION_RAYLEIGH:
        return sigma_B * 1.665;
    case CRITERION_CHIMERAX:
        return sigma_B * M_PI * M_SQRT2;
    case CRITERION_EMAN2:
        return sigma_B * M_PI * sqrt(8.0);
    default:
        return sigma_B * 1.665;
    }
}


void BFactorAnalyzer::print_stats(const ElementBStats* elem_stats, int n_elem,
    const GlobalBStats* global_stats,
    ResolutionCriterion criterion,
    MapGenerationMethod method) {

    printf("\n=== B-factor Analysis ===\n");
    printf("Method: %s\n", method_name(method));

    // Show how B-factors are used by this method
    switch (method) {
    case METHOD_EMMER:
        printf("B-factors: WILL BE USED (averaged per element + optional Refmac blur)\n");
        break;
    case METHOD_PENG1996:
        printf("B-factors: %s (use --use-bfac to enable per-element averaging)\n",
            "NOT USED by default");
        break;
    case METHOD_CHIMERAX:
        printf("B-factors: NOT USED (ChimeraX uses atomic numbers only)\n");
        break;
    case METHOD_SITUS:
        printf("B-factors: NOT USED (Situs uses mass/unit weighting only)\n");
        break;
    default:
        printf("B-factors: IGNORED\n");
    }

    printf("\nNOTE: B-factor statistics are provided for information only.\n");
    printf("      They do NOT affect map generation unless explicitly enabled.\n");

    // Per-element statistics
    printf("\nPer-element B-factor statistics:\n");
    printf("%-6s %9s %12s %12s %12s %12s %12s %16s\n",
        "Element", "Count", "Mean", "Median", "RMS", "Min", "Max", "Eq. Resolution (A)");
    printf("------------------------------------------------------------------------------------------------\n");

    for (int i = 0; i < n_elem; i++) {
        const ElementBStats* s = &elem_stats[i];
        double eq_res = bfactor_to_resolution(s->mean, criterion);

        printf("%-6s %10d %12.2f %12.2f %12.2f %12.2f %12.2f %16.2f\n",
            s->element, s->count, s->mean, s->median, s->rms, s->min, s->max, eq_res);
    }

    // Global statistics
    printf("\nGlobal B-factor statistics (all atoms):\n");
    printf("  Count: %d\n", global_stats->count);
    printf("  Mean:  %.2f A^2\n", global_stats->mean);
    printf("  Median: %.2f A^2\n", global_stats->median);
    printf("  RMS:   %.2f A^2\n", global_stats->rms);
    printf("  Min:   %.2f A^2\n", global_stats->min);
    printf("  Max:   %.2f A^2\n", global_stats->max);

    // Equivalent resolutions
    double eq_res_mean = bfactor_to_resolution(global_stats->mean, criterion);
    double eq_res_rms = bfactor_to_resolution(global_stats->rms, criterion);
    double eq_res_mean_chx = bfactor_to_resolution(global_stats->mean, CRITERION_CHIMERAX);
    double eq_res_mean_eman2 = bfactor_to_resolution(global_stats->mean, CRITERION_EMAN2);

    printf("\nEquivalent resolution based on mean B-factor:\n");
    printf("  Rayleigh criterion:  %.2f A\n", eq_res_mean);
    printf("  ChimeraX criterion:  %.2f A\n", eq_res_mean_chx);
    printf("  EMAN2 criterion:     %.2f A\n", eq_res_mean_eman2);

    printf("\nEquivalent resolution based on RMS B-factor:\n");
    printf("  Rayleigh criterion:  %.2f A\n", eq_res_rms);

    printf("============================\n\n");
}

void BFactorAnalyzer::free_stats(ElementBStats* elem_stats) {
    if (elem_stats) free(elem_stats);
}