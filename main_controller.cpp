// src/main_controller.cpp
#include "main_controller.hpp"
#include "bfactor_utils.hpp"
#include "enums.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef _WIN32
#include <windows.h>
#include <omp.h>
#else
#include <unistd.h>
#include <omp.h>
#endif

//=============================================================================
// CommandLineConfig implementation
//=============================================================================

CommandLineConfig::CommandLineConfig() {
    memset(input_file, 0, sizeof(input_file));
    memset(output_file, 0, sizeof(output_file));

    // Initialize map parameters with defaults
    map_params.resolution = 6.0;
    map_params.criterion = CRITERION_RAYLEIGH;
    map_params.amplitude_mode = AMPLITUDE_PENG1996;
    map_params.grid_spacing = 0.0;
    map_params.padding = 3.0;
    map_params.cutoff_level = 0.001;
    map_params.use_bfactors = 0;
    map_params.b_default = 20.0;

    // Initialize PDB reader config with defaults
    pdb_config.filter_hydrogen = 1;
    pdb_config.bfactor_cutoff = 0.0f;
    pdb_config.use_hetatm = 1;
    pdb_config.verbosity = 1;

    // Initialize other parameters
    cutoff_range = 5.0f;  // For ChimeraX only
    threads = 0;
    verbose = 1;
    no_normalize = 0;
    filter_water = 0;
}

void CommandLineConfig::print_usage(const char* progname) {
    printf("\npdb2mrc - Convert PDB to electron scattering map\n");
    printf("==================================================\n\n");
    printf("Usage: %s -i input.pdb -o output.mrc [options]\n\n", progname);

    printf("Required:\n");
    printf("  -i FILE    Input PDB file\n");
    printf("  -o FILE    Output MRC file\n\n");

    printf("Method Selection:\n");
    printf("  --method STR   Generation method: peng1996 (default), chimerax, situs, emmer\n\n");

    printf("Common Options:\n");
    printf("  -r FLOAT   Target resolution in A (default: 6.0)\n");
    printf("  -c STR     Resolution criterion: rayleigh (default), chimerax, eman2\n");
    printf("  -s FLOAT   Voxel size in A (default: auto = r/3)\n");
    printf("  -p FLOAT   Padding around atoms in A (default: 3.0)\n");
    printf("  -t INT     Number of OpenMP threads (0 = auto, default: 0)\n");
    printf("  -v         Verbose output\n");
    printf("  -q         Quiet mode\n");
    printf("  -h         Show this help\n\n");

    printf("Peng1996/AtomicNumber Options:\n");
    printf("  -a STR     Amplitude mode: peng1996 (default), atomic-number\n");
    printf("  --use-bfac Apply B-factors from PDB (default: off)\n");
    printf("  --b-default FLOAT Default B-factor if not in PDB (default: 20.0)\n");
    printf("  --bfactor-cutoff FLOAT   Exclude atoms with B > cutoff\n\n");

    printf("ChimeraX-Specific Options:\n");
    printf("  --cutoff FLOAT   Cutoff range in sigma (default: 5.0)\n");

    printf("Situs-Specific Options:\n");
    printf("  --situs-kernel NUM    Kernel type (1-5):\n");
    printf("                          1 = Gaussian\n");
    printf("                          2 = Triangular\n");
    printf("                          3 = Semi-Epanechnikov\n");
    printf("                          4 = Epanechnikov\n");
    printf("                          5 = Hard Sphere\n");
    printf("  --situs-kernel-type TYPE  Kernel type by name:\n");
    printf("                          gaussian, triangular, semi-epanechnikov,\n");
    printf("                          epanechnikov, hard-sphere\n");
    printf("  --situs-halfmax FLOAT   Set resolution as half-max radius (positive)\n");
    printf("  --situs-2sigma FLOAT    Set resolution as 2sigma (negative input)\n");
    printf("  --situs-margin INT      Margin voxels around structure (default: 2)\n");
    printf("  --situs-mass            Use mass weighting (default)\n");
    printf("  --situs-no-mass         Use unit weights\n");
    printf("  --situs-correction      Apply lattice smoothing correction (default)\n");
    printf("  --situs-no-correction   Skip lattice smoothing correction\n");
    printf("  --situs-amplitude FLOAT Kernel amplitude scaling (default: 1.0)\n\n");

    printf("EMmer-Specific Options:\n");
    printf("  --emmer-align         Apply MRC output alignment (flip+rotate) (default)\n");
    printf("  --emmer-no-align      Skip output alignment\n");
    printf("  --emmer-refmac-blur   Apply Refmac-compatible blur (default)\n");
    printf("  --emmer-no-blur       Skip Refmac blur\n");
    printf("  --emmer-blur FLOAT    Manual blur value in A^2 (overrides auto)\n");
    printf("  --emmer-symmetry      Apply space group symmetry (default)\n");
    printf("  --emmer-no-symmetry   Skip symmetry expansion\n");
    printf("  --emmer-cutoff FLOAT  Density cutoff for radius (default: 1e-5)\n");
    printf("  --emmer-rate FLOAT    Shannon rate for grid (default: 1.5)\n\n");

    printf("Filtering Options (all methods):\n");
    printf("  --filter-h       Filter hydrogen atoms (default: on)\n");
    printf("  --no-filter-h    Keep hydrogen atoms\n");
    printf("  --bfactor-cutoff FLOAT   Exclude atoms with B > cutoff\n\n");
    printf("Note: Water molecules (HOH, WAT, H2O, TIP) are ALWAYS filtered out.\n\n");

    printf("Examples:\n");
    printf("  # Default Peng1996 mode\n");
    printf("  pdb2mrc -i 1ake.pdb -o 1ake.mrc -r 8.0\n\n");
    printf("  # ChimeraX mode\n");
    printf("  pdb2mrc -i 1ake.pdb -o 1ake_chx.mrc -r 8.0 --method chimerax\n\n");
    printf("  # Situs mode with Gaussian kernel, half-max radius 8.0 A\n");
    printf("  pdb2mrc -i 1ake.pdb -o 1ake_situs.mrc --method situs --situs-halfmax 8.0\n\n");
    printf("  # EMmer mode with International Tables coefficients\n");
    printf("  pdb2mrc -i 1ake.pdb -o 1ake_emmer.mrc --method emmer -r 8.0\n\n");
}

static ResolutionCriterion parse_criterion(const char* s) {
    if (strcmp(s, "chimerax") == 0) return CRITERION_CHIMERAX;
    if (strcmp(s, "eman2") == 0) return CRITERION_EMAN2;
    if (strcmp(s, "rayleigh") == 0) return CRITERION_RAYLEIGH;
    // Default to Rayleigh if unknown
    fprintf(stderr, "Warning: Unknown criterion '%s', using Rayleigh\n", s);
    return CRITERION_RAYLEIGH;
}

static AmplitudeMode parse_amplitude_mode(const char* s) {
    if (strcmp(s, "atomic-number") == 0) return AMPLITUDE_ATOMIC_NUMBER;
    if (strcmp(s, "z") == 0) return AMPLITUDE_ATOMIC_NUMBER;
    if (strcmp(s, "chimerax") == 0) return AMPLITUDE_CHIMERAX;
    if (strcmp(s, "peng1996") == 0) return AMPLITUDE_PENG1996;
    return AMPLITUDE_PENG1996;
}

int CommandLineConfig::parse(int argc, char** argv) {
    // Default method
    method = METHOD_PENG1996;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 1;
        }
        else if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            STRCPY(input_file, argv[++i]);
        }
        else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            STRCPY(output_file, argv[++i]);
        }
        else if (strcmp(argv[i], "--method") == 0 && i + 1 < argc) {
            const char* m = argv[++i];
            if (strcmp(m, "peng1996") == 0 || strcmp(m, "peng") == 0) {
                method = METHOD_PENG1996;
            }
            else if (strcmp(m, "chimerax") == 0 || strcmp(m, "chx") == 0) {
                method = METHOD_CHIMERAX;
            }
            else if (strcmp(m, "situs") == 0) {
                method = METHOD_SITUS;
            }
            else if (strcmp(m, "emmer") == 0 || strcmp(m, "gemmi") == 0) {
                method = METHOD_EMMER;
            }
            else {
                fprintf(stderr, "Error: Unknown method '%s'\n", m);
                fprintf(stderr, "  Available methods: peng1996, chimerax, situs, emmer\n");
                return -1;
            }
        }
        else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
            map_params.resolution = atof(argv[++i]);
            situs_config.resolution = map_params.resolution;
            if (map_params.resolution <= 0) {
                fprintf(stderr, "Error: Resolution must be positive\n");
                return -1;
            }
        }
        else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
            map_params.criterion = parse_criterion(argv[++i]);
        }
        else if (strcmp(argv[i], "-a") == 0 && i + 1 < argc) {
            map_params.amplitude_mode = parse_amplitude_mode(argv[++i]);
        }
        else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            map_params.grid_spacing = atof(argv[++i]);
            situs_config.grid_spacing = map_params.grid_spacing;
            if (map_params.grid_spacing <= 0) {
                fprintf(stderr, "Error: Grid spacing must be positive\n");
                return -1;
            }
        }
        else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) {
            map_params.padding = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--cutoff") == 0 && i + 1 < argc) {
            cutoff_range = (float)atof(argv[++i]);
            if (cutoff_range <= 0) {
                fprintf(stderr, "Error: Cutoff range must be positive\n");
                return -1;
            }
        }
        else if (strcmp(argv[i], "--use-bfac") == 0) {
            map_params.use_bfactors = 1;
        }
        else if (strcmp(argv[i], "--b-default") == 0 && i + 1 < argc) {
            map_params.b_default = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--filter-h") == 0) {
            pdb_config.filter_hydrogen = 1;
        }
        else if (strcmp(argv[i], "--no-filter-h") == 0) {
            pdb_config.filter_hydrogen = 0;
        }
        else if (strcmp(argv[i], "--bfactor-cutoff") == 0 && i + 1 < argc) {
            pdb_config.bfactor_cutoff = (float)atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            // For backward compatibility
            pdb_config.bfactor_cutoff = (float)atof(argv[++i]);
        }
        else if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            threads = atoi(argv[++i]);
            if (threads < 0) threads = 0;
        }
        else if (strcmp(argv[i], "-v") == 0) {
            verbose = 2;
            pdb_config.verbosity = 2;
            situs_config.verbosity = 2;
        }
        else if (strcmp(argv[i], "-q") == 0) {
            verbose = 0;
            pdb_config.verbosity = 0;
            situs_config.verbosity = 0;
        }

        // Situs-specific options
        else if (strcmp(argv[i], "--situs-kernel") == 0 && i + 1 < argc) {
            int k = atoi(argv[++i]);
            if (k >= 1 && k <= 5) {
                situs_config.kernel_type = (SitusKernelType)k;
            }
            else {
                fprintf(stderr, "Error: Situs kernel type must be 1-5\n");
                return -1;
            }
        }
        else if (strcmp(argv[i], "--situs-kernel-type") == 0 && i + 1 < argc) {
            const char* type = argv[++i];
            if (strcmp(type, "gaussian") == 0)
                situs_config.kernel_type = SITUS_KERNEL_GAUSSIAN;
            else if (strcmp(type, "triangular") == 0)
                situs_config.kernel_type = SITUS_KERNEL_TRIANGULAR;
            else if (strcmp(type, "semi-epanechnikov") == 0 || strcmp(type, "semi") == 0)
                situs_config.kernel_type = SITUS_KERNEL_SEMI_EPANECHNIKOV;
            else if (strcmp(type, "epanechnikov") == 0 || strcmp(type, "epan") == 0)
                situs_config.kernel_type = SITUS_KERNEL_EPANECHNIKOV;
            else if (strcmp(type, "hard-sphere") == 0 || strcmp(type, "sphere") == 0)
                situs_config.kernel_type = SITUS_KERNEL_HARD_SPHERE;
            else {
                fprintf(stderr, "Error: Unknown kernel type '%s'\n", type);
                return -1;
            }
        }
        else if (strcmp(argv[i], "--situs-halfmax") == 0 && i + 1 < argc) {
            situs_config.resolution = fabs(atof(argv[++i]));
        }
        else if (strcmp(argv[i], "--situs-2sigma") == 0 && i + 1 < argc) {
            situs_config.resolution = -fabs(atof(argv[++i]));
        }
        else if (strcmp(argv[i], "--situs-margin") == 0 && i + 1 < argc) {
            situs_config.margin_voxels = atoi(argv[++i]);
            if (situs_config.margin_voxels < 0) {
                fprintf(stderr, "Error: Margin voxels must be >= 0\n");
                return -1;
            }
        }
        else if (strcmp(argv[i], "--situs-no-mass") == 0) {
            situs_config.use_mass_weighting = 0;
        }
        else if (strcmp(argv[i], "--situs-mass") == 0) {
            situs_config.use_mass_weighting = 1;
        }
        else if (strcmp(argv[i], "--situs-no-correction") == 0) {
            situs_config.apply_lattice_correction = 0;
        }
        else if (strcmp(argv[i], "--situs-correction") == 0) {
            situs_config.apply_lattice_correction = 1;
        }
        else if (strcmp(argv[i], "--situs-amplitude") == 0 && i + 1 < argc) {
            situs_config.kernel_amplitude = atof(argv[++i]);
            if (situs_config.kernel_amplitude <= 0) {
                fprintf(stderr, "Error: Kernel amplitude must be positive\n");
                return -1;
            }
        }

        // EMmer-specific options
        else if (strcmp(argv[i], "--emmer-align") == 0) {
            emmer_config.align_output = 1;
        }
        else if (strcmp(argv[i], "--emmer-no-align") == 0) {
            emmer_config.align_output = 0;
        }
        else if (strcmp(argv[i], "--emmer-refmac-blur") == 0) {
            emmer_config.set_refmac_blur = 1;
        }
        else if (strcmp(argv[i], "--emmer-no-blur") == 0) {
            emmer_config.set_refmac_blur = 0;
        }
        else if (strcmp(argv[i], "--emmer-blur") == 0 && i + 1 < argc) {
            emmer_config.blur = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--emmer-symmetry") == 0) {
            emmer_config.symmetry_expansion = 1;
        }
        else if (strcmp(argv[i], "--emmer-no-symmetry") == 0) {
            emmer_config.symmetry_expansion = 0;
        }
        else if (strcmp(argv[i], "--emmer-cutoff") == 0 && i + 1 < argc) {
            emmer_config.cutoff_level = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--emmer-rate") == 0 && i + 1 < argc) {
            emmer_config.rate = atof(argv[++i]);
        }
        else {
            fprintf(stderr, "Error: Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return -1;
        }
    }

    // Check required arguments
    if (input_file[0] == 0 || output_file[0] == 0) {
        fprintf(stderr, "Error: Input and output files are required\n");
        print_usage(argv[0]);
        return -1;
    }

    // Auto grid spacing if not specified
    if (map_params.grid_spacing <= 0) {
        map_params.grid_spacing = map_params.resolution / 3.0;
        situs_config.grid_spacing = map_params.grid_spacing;
        if (map_params.grid_spacing < 0.5) map_params.grid_spacing = 0.5;
    }

    // For Situs method, if resolution mode not explicitly set, use half-max radius mode
    if (method == METHOD_SITUS && situs_config.resolution == 0) {
        situs_config.resolution = map_params.resolution;
    }

    // Set OpenMP threads
    if (threads > 0) {
        omp_set_num_threads(threads);
    }

    // Print configuration if verbose
    if (verbose >= 2) {
        printf("\n=== Configuration ===\n");
        printf("Method: ");
        switch (method) {
        case METHOD_PENG1996: printf("Peng1996\n"); break;
        case METHOD_CHIMERAX: printf("ChimeraX\n"); break;
        case METHOD_SITUS: printf("Situs\n"); break;
        case METHOD_EMMER: printf("EMmer/GEMMI\n"); break;
        }
        printf("Input file: %s\n", input_file);
        printf("Output file: %s\n", output_file);
        printf("Resolution: %.2f A\n", map_params.resolution);
        printf("Criterion: %s\n", criterion_name(map_params.criterion));
        printf("Amplitude mode: %s\n", amplitude_mode_name(map_params.amplitude_mode));
        printf("Grid spacing: %.2f A\n", map_params.grid_spacing);
        printf("Padding: %.2f A\n", map_params.padding);
        printf("Use B-factors: %s\n", map_params.use_bfactors ? "yes" : "no");
        printf("Filter hydrogen: %s\n", pdb_config.filter_hydrogen ? "yes" : "no");

        if (method == METHOD_CHIMERAX) {
            printf("\nChimeraX-specific settings:\n");
            printf("  Cutoff range: %.1f sigma\n", cutoff_range);
        }

        if (method == METHOD_SITUS) {
            printf("\nSitus-specific settings:\n");
            printf("  Kernel type: ");
            switch (situs_config.kernel_type) {
            case SITUS_KERNEL_GAUSSIAN: printf("Gaussian\n"); break;
            case SITUS_KERNEL_TRIANGULAR: printf("Triangular\n"); break;
            case SITUS_KERNEL_SEMI_EPANECHNIKOV: printf("Semi-Epanechnikov\n"); break;
            case SITUS_KERNEL_EPANECHNIKOV: printf("Epanechnikov\n"); break;
            case SITUS_KERNEL_HARD_SPHERE: printf("Hard Sphere\n"); break;
            }
            printf("  Resolution mode: %s (%.2f A)\n",
                situs_config.resolution > 0 ? "half-max radius" : "2 sigma",
                fabs(situs_config.resolution));
            printf("  Margin voxels: %d\n", situs_config.margin_voxels);
            printf("  Mass weighting: %s\n", situs_config.use_mass_weighting ? "on" : "off");
            printf("  Lattice correction: %s\n", situs_config.apply_lattice_correction ? "on" : "off");
            printf("  Kernel amplitude: %.3f\n", situs_config.kernel_amplitude);
        }

        if (method == METHOD_EMMER) {
            printf("\nEMmer-specific settings:\n");
            printf("  Output alignment: %s\n", emmer_config.align_output ? "yes" : "no");
            printf("  Refmac blur: %s\n", emmer_config.set_refmac_blur ? "enabled" : "disabled");
            if (emmer_config.blur > 0)
                printf("  Manual blur: %.2f A^2\n", emmer_config.blur);
            printf("  Symmetry expansion: %s\n", emmer_config.symmetry_expansion ? "yes" : "no");
            printf("  Cutoff level: %.2e\n", emmer_config.cutoff_level);
            printf("  Shannon rate: %.2f\n", emmer_config.rate);
        }

        // B-factor usage message by method
        if (method == METHOD_EMMER) {
            printf("B-factors: WILL BE USED (averaged per element)\n");
            if (emmer_config.set_refmac_blur) {
                printf("  + Refmac blur: enabled");
                if (emmer_config.blur > 0)
                    printf(" (manual: %.2f A^2)", emmer_config.blur);
                printf("\n");
            }
        }
        else {
            printf("B-factors: IGNORED (only EMmer method uses B-factors)\n");
        }

        printf("Filter water: ALWAYS (HOH, WAT, H2O, TIP removed)\n");

        printf("====================\n\n");
    }

    return 0;
}

//=============================================================================
// PDB2MRCApp implementation
//=============================================================================

PDB2MRCApp::PDB2MRCApp() {
    m_atoms = nullptr;
    nAtoms = 0;
    vOutput = nullptr;
    m_generator = nullptr;
    m_chimerax_generator = nullptr;
    m_situs_generator = nullptr;
    m_emmer_generator = nullptr;
    nx = ny = nz = 0;
    m_origin[0] = m_origin[1] = m_origin[2] = 0.0;

    // Initialize scattering table
    Scattering::init(&m_table);
}

PDB2MRCApp::~PDB2MRCApp() {
    cleanup();
}

void PDB2MRCApp::cleanup() {
    if (m_atoms) {
        free(m_atoms);
        m_atoms = nullptr;
    }
    if (m_generator) {
        delete m_generator;
        m_generator = nullptr;
    }
    if (m_chimerax_generator) {
        delete m_chimerax_generator;
        m_chimerax_generator = nullptr;
    }
    if (m_situs_generator) {
        delete m_situs_generator;
        m_situs_generator = nullptr;
    }
    if (m_emmer_generator) {
        delete m_emmer_generator;
        m_emmer_generator = nullptr;
    }
    if (vOutput) {
        ippsFree(vOutput);
        vOutput = nullptr;
    }
}

int PDB2MRCApp::read_pdb_file() {
    if (m_config.verbose >= 1) {
        printf("Reading PDB file: %s\n", m_config.input_file);
    }
    clock_t start = clock();

    int ret = PDBReader::read_file(m_config.input_file, &m_config.pdb_config,
        &m_atoms, &nAtoms);

    clock_t end = clock();
    double cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;

    if (ret != 0) {
        fprintf(stderr, "Error reading PDB file: %s\n", PDBReader::get_error());
        return ret;
    }

    if (m_config.verbose >= 1) {
        printf("Read %d atoms in %.2f seconds\n", nAtoms, cpu_time);
    }

    if (nAtoms == 0) {
        fprintf(stderr, "Error: No atoms found after filtering\n");
        return -3;
    }

    // Add B-factor analysis for all methods
    ElementBStats* elem_stats = NULL;
    int n_elem_stats = 0;
    GlobalBStats global_stats;

    int ret_bf = BFactorAnalyzer::calculate_per_element_stats(m_atoms, nAtoms,
        &elem_stats, &n_elem_stats);
    if (ret_bf == 0 && elem_stats) {
        BFactorAnalyzer::calculate_global_stats(m_atoms, nAtoms, &global_stats);

        // Print B-factor statistics
        BFactorAnalyzer::print_stats(elem_stats, n_elem_stats, &global_stats,
            m_config.map_params.criterion,
            m_config.method);

        BFactorAnalyzer::free_stats(elem_stats);
    }

    return 0;
}

int PDB2MRCApp::map_elements_to_indices() {
    int unknown_elements = 0;

    for (int i = 0; i < nAtoms; i++) {
        if (m_config.map_params.amplitude_mode != AMPLITUDE_CHIMERAX) {
            // Only needed for Peng1996/AtomicNumber modes
            m_atoms[i].element_idx = Scattering::get_index(&m_table, m_atoms[i].element);
            if (m_atoms[i].element_idx < 0) {
                if (m_config.verbose >= 1) {
                    fprintf(stderr, "Warning: Unknown element %s, using carbon\n",
                        m_atoms[i].element);
                }
                m_atoms[i].element_idx = Scattering::get_index(&m_table, "C");
                unknown_elements++;
            }
        }
    }

    if (unknown_elements > 0 && m_config.verbose >= 1) {
        printf("Warning: %d atoms with unknown elements (set to C)\n", unknown_elements);
    }

    return 0;
}

int PDB2MRCApp::initialize_generator() {
    // Clean up any existing generators
    if (m_generator) {
        delete m_generator;
        m_generator = nullptr;
    }
    if (m_chimerax_generator) {
        delete m_chimerax_generator;
        m_chimerax_generator = nullptr;
    }
    if (m_situs_generator) {
        delete m_situs_generator;
        m_situs_generator = nullptr;
    }
    if (m_emmer_generator) {
        delete m_emmer_generator;
        m_emmer_generator = nullptr;
    }

    // Use the method to decide which generator to create
    if (m_config.method == METHOD_CHIMERAX) {
        // Initialize ChimeraX generator
        m_chimerax_generator = new ChimeraXGenerator();
        if (!m_chimerax_generator) {
            fprintf(stderr, "Error: Failed to create ChimeraX generator\n");
            return -1;
        }

        int ret = m_chimerax_generator->init(
            m_atoms, nAtoms,
            m_config.map_params.resolution,
            m_config.map_params.grid_spacing,
            m_config.map_params.padding,
            m_config.cutoff_range
        );

        if (ret != 0) {
            fprintf(stderr, "Error: Failed to initialize ChimeraX generator (error %d)\n", ret);
            delete m_chimerax_generator;
            m_chimerax_generator = nullptr;
            return ret;
        }

        if (m_config.verbose >= 1) {
            printf("Initialized ChimeraX generator\n");
        }
    }
    else if (m_config.method == METHOD_SITUS) {
        // Initialize Situs generator
        m_situs_generator = new SitusGenerator();
        if (!m_situs_generator) {
            fprintf(stderr, "Error: Failed to create Situs generator\n");
            return -1;
        }

        int ret = m_situs_generator->init(&m_config.situs_config, m_atoms, nAtoms);
        if (ret != 0) {
            fprintf(stderr, "Error: Failed to initialize Situs generator (error %d)\n", ret);
            delete m_situs_generator;
            m_situs_generator = nullptr;
            return ret;
        }

        if (m_config.verbose >= 1) {
            printf("Initialized Situs generator\n");
        }
    }
    else if (m_config.method == METHOD_EMMER) {
        // Initialize EMmer generator
        m_emmer_generator = new EmmerGenerator();
        if (!m_emmer_generator) {
            fprintf(stderr, "Error: Failed to create EMmer generator\n");
            return -1;
        }

        int ret = m_emmer_generator->init(
            &m_config.emmer_config,
            m_atoms,
            nAtoms,
            m_config.map_params.resolution,
            m_config.map_params.grid_spacing,
            nullptr  // cell_dimensions
        );

        if (ret != 0) {
            fprintf(stderr, "Error: Failed to initialize EMmer generator (error %d)\n", ret);
            delete m_emmer_generator;
            m_emmer_generator = nullptr;
            return ret;
        }

        if (m_config.verbose >= 1) {
            printf("Initialized EMmer/GEMMI generator\n");
        }
    }
    else {
        // Initialize Peng1996/AtomicNumber generator (default)
        m_generator = new MapGenerator();
        if (!m_generator) {
            fprintf(stderr, "Error: Failed to create map generator\n");
            return -1;
        }

        int ret = m_generator->init(&m_config.map_params, m_atoms, nAtoms);
        if (ret != 0) {
            fprintf(stderr, "Error: Failed to initialize map generator (error %d)\n", ret);
            delete m_generator;
            m_generator = nullptr;
            return ret;
        }

        if (m_config.verbose >= 1) {
            printf("Initialized Peng1996/AtomicNumber generator\n");
        }
    }

    return 0;
}

int PDB2MRCApp::generate_map() {
    clock_t start = clock();
    int ret = 0;

    if (m_config.method == METHOD_CHIMERAX) {
        if (!m_chimerax_generator) {
            fprintf(stderr, "Error: ChimeraX generator not initialized\n");
            return -1;
        }

        if (m_config.verbose >= 1) {
            printf("\nGenerating ChimeraX-style map at %.2f A resolution\n",
                m_config.map_params.resolution);
            printf("  Amplitude mode: %s\n", amplitude_mode_name(m_config.map_params.amplitude_mode));
            printf("  Criterion: %s\n", criterion_name(m_config.map_params.criterion));
            printf("  Grid spacing: %.2f A, Padding: %.2f A\n",
                m_config.map_params.grid_spacing, m_config.map_params.padding);
            printf("  sigma = R/(pi sqrt(2)) = %.3f A\n",
                m_config.map_params.resolution / (M_PI * M_SQRT2));
            printf("  Cutoff range: %.1f sigma\n", m_config.cutoff_range);
        }

        ret = m_chimerax_generator->run();

        if (ret == 0) {
            vOutput = m_chimerax_generator->releaseOutput();
            m_chimerax_generator->getDimensions(nx, ny, nz);
            m_chimerax_generator->getOrigin(m_origin);
            m_config.map_params.grid_spacing = m_chimerax_generator->getStep();
        }
    }
    else if (m_config.method == METHOD_SITUS) {
        if (!m_situs_generator) {
            fprintf(stderr, "Error: Situs generator not initialized\n");
            return -1;
        }

        if (m_config.verbose >= 1) {
            printf("\nGenerating Situs-style map\n");
            printf("  Kernel type: %d\n", m_config.situs_config.kernel_type);
            printf("  Resolution: %.2f %s\n",
                fabs(m_config.situs_config.resolution),
                (m_config.situs_config.resolution < 0) ? "(2sigma mode)" : "(half-max radius mode)");
            printf("  Grid spacing: %.2f A\n", m_config.situs_config.grid_spacing);
        }

        ret = m_situs_generator->run();

        if (ret == 0) {
            vOutput = m_situs_generator->releaseOutput();
            nx = m_situs_generator->getOutputNx();
            ny = m_situs_generator->getOutputNy();
            nz = m_situs_generator->getOutputNz();
            m_situs_generator->getOutputOrigin(m_origin);
            m_config.map_params.grid_spacing = m_situs_generator->getGridSpacing();
        }
    }
    else if (m_config.method == METHOD_EMMER) {
        if (!m_emmer_generator) {
            fprintf(stderr, "Error: EMmer generator not initialized\n");
            return -1;
        }

        if (m_config.verbose >= 1) {
            printf("\nGenerating EMmer-style map at %.2f A resolution\n",
                m_config.map_params.resolution);
            printf("  Using International Tables c4322 coefficients\n");
            printf("  Grid spacing: %.2f A\n", m_config.map_params.grid_spacing);
            printf("  Refmac blur: %s\n",
                m_config.emmer_config.set_refmac_blur ? "enabled" : "disabled");
            if (m_config.emmer_config.blur > 0)
                printf("  Manual blur: %.2f A^2\n", m_config.emmer_config.blur);
            printf("  Output alignment: %s\n",
                m_config.emmer_config.align_output ? "yes" : "no");
            printf("  Symmetry expansion: %s\n",
                m_config.emmer_config.symmetry_expansion ? "yes" : "no");
        }

        ret = m_emmer_generator->run();

        if (ret == 0) {
            vOutput = m_emmer_generator->releaseOutput();
            nx = m_emmer_generator->getOutputNx();
            ny = m_emmer_generator->getOutputNy();
            nz = m_emmer_generator->getOutputNz();
            m_emmer_generator->getOutputOrigin(m_origin);
            m_config.map_params.grid_spacing = m_emmer_generator->getGridSpacing();
        }
    }
    else {
        if (!m_generator) {
            fprintf(stderr, "Error: Map generator not initialized\n");
            return -1;
        }

        if (m_config.verbose >= 1) {
            printf("\nGenerating map at %.2f A resolution (%s)\n",
                m_config.map_params.resolution,
                criterion_name(m_config.map_params.criterion));
            printf("  Amplitude mode: %s\n", amplitude_mode_name(m_config.map_params.amplitude_mode));
            printf("  Grid spacing: %.2f A, Padding: %.2f A\n",
                m_config.map_params.grid_spacing, m_config.map_params.padding);
            printf("  Use B-factors: %s\n", m_config.map_params.use_bfactors ? "yes" : "no");
        }

        ret = m_generator->run();

        if (ret == 0) {
            vOutput = m_generator->releaseOutputMap();
            nx = m_generator->getOutputNx();
            ny = m_generator->getOutputNy();
            nz = m_generator->getOutputNz();
            m_generator->getOutputOrigin(m_origin);
        }
    }

    clock_t end = clock();
    double cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;

    if (ret != 0) {
        fprintf(stderr, "Error: Map generation failed with error %d\n", ret);
        return ret;
    }

    if (m_config.verbose >= 1) {
        printf("Map generated in %.2f seconds\n", cpu_time);
        printf("Grid size: %d x %d x %d (%lld voxels)\n",
            nx, ny, nz, (long long)nx * ny * nz);
        printf("Origin: (%.2f, %.2f, %.2f)\n", m_origin[0], m_origin[1], m_origin[2]);
    }

    return 0;
}

int PDB2MRCApp::write_output() {
    int64_t nvox = (int64_t)nx * ny * nz;

    // Convert double to float for MRC output
    Ipp32f* float_output = ippsMalloc_32f(nvox);
    if (!float_output) {
        fprintf(stderr, "Error: Memory allocation failed for output conversion\n");
        return -6;
    }

    IppStatus status = ippsConvert_64f32f(vOutput, float_output, nvox);
    if (status != ippStsNoErr) {
        fprintf(stderr, "Error: Failed to convert map to float: %d\n", status);
        ippsFree(float_output);
        return -7;
    }

    // Write MRC file
    float origin_float[3] = { (float)m_origin[0], (float)m_origin[1], (float)m_origin[2] };
    int ret = MRCWriter::write(m_config.output_file, float_output,
        nx, ny, nz,
        origin_float, (float)m_config.map_params.grid_spacing);

    ippsFree(float_output);

    if (ret != 0) {
        fprintf(stderr, "Error: Failed to write MRC file: %s\n", MRCWriter::get_error());
        return ret;
    }

    if (m_config.verbose >= 1) {
        printf("Map written to %s\n", m_config.output_file);
    }

    return 0;
}

void PDB2MRCApp::print_summary() {
    if (m_config.verbose >= 1) {
        printf("\n=== Processing Complete ===\n");
        printf("Input: %s\n", m_config.input_file);
        printf("Output: %s\n", m_config.output_file);
        printf("Atoms: %d\n", nAtoms);
        printf("Resolution: %.2f A (%s)\n",
            m_config.map_params.resolution,
            criterion_name(m_config.map_params.criterion));
        printf("Mode: %s\n", amplitude_mode_name(m_config.map_params.amplitude_mode));
        printf("Grid: %d x %d x %d voxels (%.2f A spacing)\n",
            nx, ny, nz, m_config.map_params.grid_spacing);
        printf("===========================\n");
    }
}

int PDB2MRCApp::run(int argc, char** argv) {
    // Parse command line
    int parse_result = m_config.parse(argc, argv);
    if (parse_result != 0) {
        return (parse_result == 1) ? 0 : 1;  // 1 = help requested, exit gracefully
    }

    // Set verbosity from config
    if (m_config.verbose >= 2) {
        printf("Command line parsed successfully\n");
        printf("Input: %s\n", m_config.input_file);
        printf("Output: %s\n", m_config.output_file);
        printf("Resolution: %.2f A\n", m_config.map_params.resolution);
        printf("Mode: %s\n", amplitude_mode_name(m_config.map_params.amplitude_mode));
    }

    // Read PDB file
    int ret = read_pdb_file();
    if (ret != 0) return 2;

    // Map elements to indices (not needed for ChimeraX but harmless)
    ret = map_elements_to_indices();
    if (ret != 0) return 3;

    // Initialize appropriate generator
    ret = initialize_generator();
    if (ret != 0) return 4;

    // Generate map
    ret = generate_map();
    if (ret != 0) return 5;

    // Write output
    ret = write_output();
    if (ret != 0) return 6;

    // Print summary
    print_summary();

    return 0;
}