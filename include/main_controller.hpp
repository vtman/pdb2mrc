// include/main_controller.hpp
#pragma once

#include "types.hpp"
#include "atom.hpp"
#include "pdb_reader.hpp"
#include "mrc_writer.hpp"
#include "scattering.hpp"
#include "map_generator.hpp"
#include "chimerax_generator.hpp"
#include "situs_generator.hpp"
#include "emmer_generator.hpp"  // New header

// Command line configuration
struct CommandLineConfig {
    char input_file[1024];
    char output_file[1024];
    MapParams map_params;
    PDBReaderConfig pdb_config;
    float cutoff_range;
    int threads;
    int verbose;
    int no_normalize;
    int filter_water;
    MapGenerationMethod method;
    SitusGeneratorConfig situs_config;
    EmmerGeneratorConfig emmer_config;  // New

    CommandLineConfig();
    void print_usage(const char* progname);
    int parse(int argc, char** argv);
};

// Main application controller
class PDB2MRCApp {
private:
    CommandLineConfig m_config;
    ScatteringTable m_table;
    Atom* m_atoms;
    int nAtoms;
    Ipp64f* vOutput;
    int nx, ny, nz;
    Ipp64f m_origin[3];

    // Map generators
    MapGenerator* m_generator;
    ChimeraXGenerator* m_chimerax_generator;
    SitusGenerator* m_situs_generator;
    EmmerGenerator* m_emmer_generator;   // New

    int read_pdb_file();
    int map_elements_to_indices();
    int initialize_generator();
    int generate_map();
    int write_output();
    void print_summary();
    void cleanup();

public:
    PDB2MRCApp();
    ~PDB2MRCApp();

    int run(int argc, char** argv);

    PDB2MRCApp(const PDB2MRCApp&) = delete;
    PDB2MRCApp& operator=(const PDB2MRCApp&) = delete;
};