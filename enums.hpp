// include/enums.hpp
#pragma once

// Resolution criterion definitions
enum ResolutionCriterion {
    CRITERION_RAYLEIGH = 0,
    CRITERION_CHIMERAX,
    CRITERION_EMAN2,
    CRITERION_FSC_0143,
    CRITERION_FSC_05,
    CRITERION_DEFAULT = CRITERION_RAYLEIGH
};

enum MapGenerationMethod {
    METHOD_PENG1996 = 0,
    METHOD_CHIMERAX,
    METHOD_SITUS,
    METHOD_EMMER      // New method
};

// Amplitude mode definitions
enum AmplitudeMode {
    AMPLITUDE_PENG1996 = 0,
    AMPLITUDE_ATOMIC_NUMBER,
    AMPLITUDE_CHIMERAX,
    AMPLITUDE_EMMER    // New mode
};

const char* criterion_name(ResolutionCriterion criterion);
const char* amplitude_mode_name(AmplitudeMode mode);
const char* method_name(MapGenerationMethod method);