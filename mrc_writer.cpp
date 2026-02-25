// src/mrc_writer.cpp
#include "mrc_writer.hpp"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

char MRCWriter::m_error_msg[256] = { 0 };

void MRCWriter::init_header(MRHeader* header,
    int nx, int ny, int nz,
    const float origin[3],
    float pixel_size) {

    if (!header) return;

    memset(header, 0, sizeof(MRHeader));

    header->nx = nx;
    header->ny = ny;
    header->nz = nz;
    header->mode = 2;  // 32-bit float

    header->nx_start = 0;
    header->ny_start = 0;
    header->nz_start = 0;

    header->mx = nx;
    header->my = ny;
    header->mz = nz;

    header->cella[0] = nx * pixel_size;
    header->cella[1] = ny * pixel_size;
    header->cella[2] = nz * pixel_size;

    header->cellb[0] = 90.0f;
    header->cellb[1] = 90.0f;
    header->cellb[2] = 90.0f;

    header->mapc = 1;
    header->mapr = 2;
    header->maps = 3;

    header->ispg = 1;  // P1
    header->nsymbt = 0;

    // Clear extra space
    memset(header->extra, 0, sizeof(header->extra));

    header->origin[0] = origin[0];
    header->origin[1] = origin[1];
    header->origin[2] = origin[2];

    header->map[0] = 'M';
    header->map[1] = 'A';
    header->map[2] = 'P';
    header->map[3] = ' ';

    // Machine stamp: 0x4444 = little-endian
    header->machstamp[0] = 0x4444;
    header->machstamp[1] = 0x0000;

    header->version = 20140;  // MRC2014

    // Add default label
    strcpy(header->labels, "Created by pdb2mrc");
    header->nlabl = 1;
}

const char* MRCWriter::get_error() {
    return m_error_msg;
}

int MRCWriter::write(const char* filename,
    const Ipp32f* data,
    int nx, int ny, int nz,
    const float origin[3],
    float pixel_size) {

    MRHeader header;
    init_header(&header, nx, ny, nz, origin, pixel_size);

    // Compute statistics
    IppStatus status;
    status = ippsMinMax_32f(data, nx * ny * nz, &header.dmin, &header.dmax);
    if (status != ippStsNoErr) {
        snprintf(m_error_msg, sizeof(m_error_msg),
            "Failed to compute min/max: %d", status);
        return -1;
    }

    status = ippsMean_32f(data, nx * ny * nz, &header.dmean, ippAlgHintAccurate);
    if (status != ippStsNoErr) {
        snprintf(m_error_msg, sizeof(m_error_msg),
            "Failed to compute mean: %d", status);
        return -2;
    }

    Ipp32f std;
    status = ippsStdDev_32f(data, nx * ny * nz, &std, ippAlgHintAccurate);
    if (status != ippStsNoErr) {
        snprintf(m_error_msg, sizeof(m_error_msg),
            "Failed to compute std dev: %d", status);
        return -3;
    }
    header.rms = std;

    return write_with_header(filename, data, &header);
}

int MRCWriter::write_with_header(const char* filename,
    const Ipp32f* data,
    const MRHeader* header) {

    FILE* fp = fopen(filename, "wb");
    if (!fp) {
        snprintf(m_error_msg, sizeof(m_error_msg),
            "Cannot open file: %s", filename);
        return -1;
    }

    // Write header
    size_t written = fwrite(header, 1, 1024, fp);
    if (written != 1024) {
        fclose(fp);
        snprintf(m_error_msg, sizeof(m_error_msg),
            "Failed to write header (wrote %zu of 1024 bytes)", written);
        return -2;
    }

    // Write data
    int nvox = header->nx * header->ny * header->nz;
    written = fwrite(data, sizeof(float), nvox, fp);
    if (written != (size_t)nvox) {
        fclose(fp);
        snprintf(m_error_msg, sizeof(m_error_msg),
            "Failed to write data (wrote %zu of %d elements)", written, nvox);
        return -3;
    }

    fclose(fp);
    return 0;
}