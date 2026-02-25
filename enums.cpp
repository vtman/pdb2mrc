// src/enums.cpp
#include "enums.hpp"
#include <math.h>

const char* criterion_name(ResolutionCriterion criterion) {
    switch (criterion) {
    case CRITERION_RAYLEIGH: return "Rayleigh";
    case CRITERION_CHIMERAX: return "ChimeraX";
    case CRITERION_EMAN2:    return "EMAN2";
    case CRITERION_FSC_0143: return "FSC=0.143";
    case CRITERION_FSC_05:   return "FSC=0.5";
    default:                  return "Unknown";
    }
}

const char* amplitude_mode_name(AmplitudeMode mode) {
    switch (mode) {
    case AMPLITUDE_PENG1996:      return "Peng1996 (fe0)";
    case AMPLITUDE_ATOMIC_NUMBER: return "Atomic Number (Z) - EMAN2 style";
    case AMPLITUDE_CHIMERAX:       return "ChimeraX molmap";
    case AMPLITUDE_EMMER:        return "EMmer (IT c4322)";
    default:                       return "Unknown";
    }
}