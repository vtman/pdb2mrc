#pragma once
#pragma once

/**
 * Cross-platform definitions for pdb2mrc
 * Supports Windows and Linux/Unix
 */

#ifdef _WIN32
 // Windows-specific includes
#define NOMINMAX
#include <windows.h>

// File I/O
#define FOPEN_READ "rb"
#define FOPEN_WRITE "wb"
#define FOPEN_READ_TEXT "r"
#define FOPEN_WRITE_TEXT "w"

// Path separator
#define PATH_SEPARATOR '\\'
#define PATH_SEPARATOR_STR "\\"

// Safe string functions
#define STRCPY(dest, src) strcpy_s(dest, sizeof(dest), src)
#define STRNCAT(dest, src, n) strncat_s(dest, sizeof(dest), src, n)
#define SPRINTF(dest, ...) sprintf_s(dest, sizeof(dest), __VA_ARGS__)
#define STRTOK_R(str, delim, saveptr) strtok_s(str, delim, saveptr)

// File offset
typedef __int64 file_offset_t;
#define FSEEK _fseeki64
#define FTELL _ftelli64

#else
 // Linux/Unix includes
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

// File I/O
#define FOPEN_READ "r"
#define FOPEN_WRITE "w"
#define FOPEN_READ_TEXT "r"
#define FOPEN_WRITE_TEXT "w"

// Path separator
#define PATH_SEPARATOR '/'
#define PATH_SEPARATOR_STR "/"

// Safe string functions
#define STRCPY(dest, src) strcpy(dest, src)
#define STRNCAT(dest, src, n) strncat(dest, src, n)
#define SPRINTF(dest, ...) sprintf(dest, __VA_ARGS__)
#define STRTOK_R(str, delim, saveptr) strtok_r(str, delim, saveptr)

// File offset
typedef off_t file_offset_t;
#define FSEEK fseeko
#define FTELL ftello

#endif

// Common integer types
#include <stdint.h>

// Endian detection
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define IS_BIG_ENDIAN 1
#elif defined(_WIN32)
#define IS_BIG_ENDIAN 0  // Windows is little-endian
#else
#include <endian.h>
#if __BYTE_ORDER == __BIG_ENDIAN
#define IS_BIG_ENDIAN 1
#else
#define IS_BIG_ENDIAN 0
#endif
#endif

// Alignment macros
#ifdef _WIN32
#define ALIGNED_MALLOC(size, align) _aligned_malloc(size, align)
#define ALIGNED_FREE(ptr) _aligned_free(ptr)
#else
#define ALIGNED_MALLOC(size, align) aligned_alloc(align, size)
#define ALIGNED_FREE(ptr) free(ptr)
#endif

// DLL export/import
#ifdef _WIN32
#ifdef PDB2MRC_EXPORTS
#define PDB2MRC_API __declspec(dllexport)
#else
#define PDB2MRC_API __declspec(dllimport)
#endif
#else
#define PDB2MRC_API __attribute__((visibility("default")))
#endif