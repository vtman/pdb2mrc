// include/mrc_writer.hpp
#pragma once

#include "types.hpp"
#include <stdint.h>
#include <stdio.h>

// MRC/CCP4 header structure
// Based on MRC2014 standard
struct MRHeader {
    int32_t nx, ny, nz;           // Grid dimensions
    int32_t mode;                  // 2 = 32-bit float
    int32_t nx_start, ny_start, nz_start;  // Starting point
    int32_t mx, my, mz;            // Grid sampling (same as nx,ny,nz)
    float cella[3];                 // Cell dimensions (Å)
    float cellb[3];                 // Cell angles (degrees)
    int32_t mapc, mapr, maps;       // Column, row, section axes
    float dmin, dmax, dmean;        // Density statistics
    int32_t ispg;                    // Space group (1 = P1)
    int32_t nsymbt;                  // Extended header size
    int32_t extra[25];               // Extra space
    float origin[3];                 // Origin in Å
    char map[4];                     // Must be "MAP "
    int32_t machstamp[2];            // Machine stamp
    float rms;                       // RMS deviation
    int32_t nlabl;                    // Number of labels
    char labels[800];                 // Text labels
    int32_t version;                  // Version (20140 for MRC2014)
};

class MRCWriter {
private:
    static char m_error_msg[256];

public:
    static void init_header(MRHeader* header,
        int nx, int ny, int nz,
        const float origin[3],
        float pixel_size);

    static const char* get_error();

    static int write(const char* filename,
        const Ipp32f* data,  // Use Ipp32f consistently
        int nx, int ny, int nz,
        const float origin[3],
        float pixel_size);

    static int write_with_header(const char* filename,
        const Ipp32f* data,  // Use Ipp32f consistently
        const MRHeader* header);
};