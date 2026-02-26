# cmake/FindMKL.cmake
# Find Intel Math Kernel Library (MKL)

include(FindPackageHandleStandardArgs)

# Set default search paths
set(MKL_ROOT_DIR "" CACHE PATH "Root directory for Intel MKL")

# Check environment variables
if(NOT MKL_ROOT_DIR)
    if(DEFINED ENV{MKLROOT})
        set(MKL_ROOT_DIR $ENV{MKLROOT})
    elseif(DEFINED ENV{INTEL_MKL_ROOT})
        set(MKL_ROOT_DIR $ENV{INTEL_MKL_ROOT})
    endif()
endif()

# Common MKL library names
set(MKL_COMMON_LIBRARIES
    mkl_intel_lp64
    mkl_sequential
    mkl_core
)

set(MKL_COMMON_LIBRARIES_ILP64
    mkl_intel_ilp64
    mkl_sequential
    mkl_core
)

set(MKL_COMMON_LIBRARIES_OPENMP
    mkl_intel_lp64
    mkl_intel_thread
    mkl_core
)

set(MKL_COMMON_LIBRARIES_OPENMP_ILP64
    mkl_intel_ilp64
    mkl_intel_thread
    mkl_core
)

set(MKL_COMMON_LIBRARIES_TBB
    mkl_intel_lp64
    mkl_tbb_thread
    mkl_core
)

# Find MKL libraries
find_path(MKL_INCLUDE_DIR
    NAMES mkl.h
    PATHS
        ${MKL_ROOT_DIR}/include
        ${MKL_ROOT_DIR}/include/intel64/lp64
        /opt/intel/mkl/include
        /opt/intel/oneapi/mkl/latest/include
        /usr/local/intel/mkl/include
        /usr/include/mkl
    PATH_SUFFIXES include
)

# Try to find libraries in different configurations
foreach(lib IN LISTS MKL_COMMON_LIBRARIES)
    find_library(MKL_LIBRARY_${lib}
        NAMES ${lib} lib${lib}
        PATHS
            ${MKL_ROOT_DIR}/lib
            ${MKL_ROOT_DIR}/lib/intel64
            ${MKL_ROOT_DIR}/lib/intel64_lin
            /opt/intel/mkl/lib/intel64
            /opt/intel/oneapi/mkl/latest/lib/intel64
            /usr/local/intel/mkl/lib/intel64
            /usr/lib/x86_64-linux-gnu
        PATH_SUFFIXES lib
    )
    
    if(MKL_LIBRARY_${lib})
        list(APPEND MKL_LIBRARIES ${MKL_LIBRARY_${lib}})
    endif()
endforeach()

# Handle Windows platform
if(WIN32)
    # Add Windows-specific paths and libraries
    foreach(lib IN LISTS MKL_COMMON_LIBRARIES)
        find_library(MKL_LIBRARY_${lib}
            NAMES ${lib} lib${lib}
            PATHS
                ${MKL_ROOT_DIR}/lib/intel64_win
                ${MKL_ROOT_DIR}/lib/intel64
                "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64_win"
                "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/intel64_win"
            PATH_SUFFIXES lib
        )
    endforeach()
endif()

# Handle macOS platform
if(APPLE)
    # Add macOS-specific paths
    foreach(lib IN LISTS MKL_COMMON_LIBRARIES)
        find_library(MKL_LIBRARY_${lib}
            NAMES ${lib} lib${lib}
            PATHS
                ${MKL_ROOT_DIR}/lib
                /opt/intel/mkl/lib
                /opt/intel/oneapi/mkl/latest/lib
            PATH_SUFFIXES lib
        )
    endforeach()
endif()

# Set include directories
set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})

# Handle standard arguments
find_package_handle_standard_args(MKL
    REQUIRED_VARS MKL_INCLUDE_DIR MKL_LIBRARIES
    VERSION_VAR MKL_VERSION
)

if(MKL_FOUND AND NOT TARGET MKL::MKL)
    add_library(MKL::MKL UNKNOWN IMPORTED)
    set_target_properties(MKL::MKL PROPERTIES
        IMPORTED_LOCATION "${MKL_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIRS}"
    )
endif()

mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARIES)