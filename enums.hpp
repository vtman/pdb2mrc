// include/enums.hpp
#pragma once

// Resolution criterion definitions - matching README
enum ResolutionCriterion {
    CRITERION_RAYLEIGH = 0,    // σ = R/1.665
    CRITERION_CHIMERAX,        // σ = R/(π√2)
    CRITERION_EMAN2,           // σ = R/(π√8)
    CRITERION_DEFAULT = CRITERION_RAYLEIGH
};

enum MapGenerationMethod {
    METHOD_PENG1996 = 0,       // Default method with convolution chain
    METHOD_CHIMERAX,            // ChimeraX molmap algorithm
    METHOD_SITUS,               // Situs multi-kernel approach
    METHOD_EMMER                // EMmer/GEMMI with International Tables
};

// Amplitude mode definitions
enum AmplitudeMode {
    AMPLITUDE_PENG1996 = 0,     // Peng1996 fe0 values
    AMPLITUDE_ATOMIC_NUMBER,     // Atomic number Z (EMAN2 style)
    AMPLITUDE_CHIMERAX,          // ChimeraX molmap (always Z)
    AMPLITUDE_EMMER              // EMmer (International Tables)
};

const char* criterion_name(ResolutionCriterion criterion);
const char* amplitude_mode_name(AmplitudeMode mode);
const char* method_name(MapGenerationMethod method);