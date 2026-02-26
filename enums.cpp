// src/enums.cpp
#include "enums.hpp"

const char* criterion_name(ResolutionCriterion criterion) {
    switch (criterion) {
    case CRITERION_RAYLEIGH: return "Rayleigh (sigma = R/1.665)";
    case CRITERION_CHIMERAX: return "ChimeraX (sigma = R/(π√2))";
    case CRITERION_EMAN2:    return "EMAN2 (sigma = R/(π√8))";
    default:                  return "Unknown";
    }
}

const char* amplitude_mode_name(AmplitudeMode mode) {
    switch (mode) {
    case AMPLITUDE_PENG1996:      return "Peng1996 (fe0)";
    case AMPLITUDE_ATOMIC_NUMBER: return "Atomic Number (Z)";
    case AMPLITUDE_CHIMERAX:      return "ChimeraX molmap (Z)";
    case AMPLITUDE_EMMER:         return "EMmer (IT c4322)";
    default:                       return "Unknown";
    }
}

const char* method_name(MapGenerationMethod method) {
    switch (method) {
    case METHOD_PENG1996: return "Peng1996 (Default)";
    case METHOD_CHIMERAX: return "ChimeraX molmap";
    case METHOD_SITUS:    return "Situs";
    case METHOD_EMMER:    return "EMmer/GEMMI";
    default:              return "Unknown";
    }
}