# cmake/FindIPP.cmake
# Find Intel Integrated Performance Primitives (IPP)

include(FindPackageHandleStandardArgs)

# Set default search paths
set(IPP_ROOT_DIR "" CACHE PATH "Root directory for Intel IPP")

# Check environment variables
if(NOT IPP_ROOT_DIR)
    if(DEFINED ENV{IPPROOT})
        set(IPP_ROOT_DIR $ENV{IPPROOT})
    elseif(DEFINED ENV{INTEL_IPP_ROOT})
        set(IPP_ROOT_DIR $ENV{INTEL_IPP_ROOT})
    endif()
endif()

# Architecture detection
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(IPP_ARCH "intel64")
else()
    set(IPP_ARCH "ia32")
endif()

# Common IPP components
set(IPP_COMPONENTS
    ippcore
    ipps
    ippi
    ippvm
    ippcc
    ippcv
    ippch
)

# Find IPP include directory
find_path(IPP_INCLUDE_DIR
    NAMES ipp.h
    PATHS
        ${IPP_ROOT_DIR}/include
        ${IPP_ROOT_DIR}/include/${IPP_ARCH}
        /opt/intel/ipp/include
        /opt/intel/oneapi/ipp/latest/include
        /opt/intel/compilers_and_libraries/linux/ipp/include
        /usr/local/intel/ipp/include
    PATH_SUFFIXES include
)

# Find IPP libraries
foreach(comp IN LISTS IPP_COMPONENTS)
    # Try different naming conventions
    find_library(IPP_LIBRARY_${comp}
        NAMES
            ${comp}
            lib${comp}
            ${comp}mt        # Multi-threaded
            lib${comp}mt
            ${comp}${IPP_ARCH}
            lib${comp}${IPP_ARCH}
        PATHS
            ${IPP_ROOT_DIR}/lib
            ${IPP_ROOT_DIR}/lib/${IPP_ARCH}
            ${IPP_ROOT_DIR}/lib/intel64
            /opt/intel/ipp/lib/${IPP_ARCH}
            /opt/intel/oneapi/ipp/latest/lib/${IPP_ARCH}
            /opt/intel/compilers_and_libraries/linux/ipp/lib/${IPP_ARCH}
            /usr/local/intel/ipp/lib/${IPP_ARCH}
        PATH_SUFFIXES lib
    )
    
    if(IPP_LIBRARY_${comp})
        list(APPEND IPP_LIBRARIES ${IPP_LIBRARY_${comp}})
    endif()
endforeach()

# Handle Windows platform
if(WIN32)
    foreach(comp IN LISTS IPP_COMPONENTS)
        find_library(IPP_LIBRARY_${comp}
            NAMES
                ${comp}mt
                ${comp}
                lib${comp}mt
                lib${comp}
            PATHS
                ${IPP_ROOT_DIR}/lib/${IPP_ARCH}
                ${IPP_ROOT_DIR}/lib/intel64_win
                "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/ipp/lib/${IPP_ARCH}"
                "C:/Program Files (x86)/Intel/oneAPI/ipp/latest/lib/${IPP_ARCH}"
            PATH_SUFFIXES lib
        )
    endforeach()
endif()

# Handle macOS platform
if(APPLE)
    foreach(comp IN LISTS IPP_COMPONENTS)
        find_library(IPP_LIBRARY_${comp}
            NAMES
                ${comp}
                lib${comp}
            PATHS
                ${IPP_ROOT_DIR}/lib
                /opt/intel/ipp/lib
                /opt/intel/oneapi/ipp/latest/lib
            PATH_SUFFIXES lib
        )
    endforeach()
endif()

set(IPP_INCLUDE_DIRS ${IPP_INCLUDE_DIR})

# Handle standard arguments
find_package_handle_standard_args(IPP
    REQUIRED_VARS IPP_INCLUDE_DIR IPP_LIBRARIES
)

if(IPP_FOUND AND NOT TARGET IPP::IPP)
    add_library(IPP::IPP UNKNOWN IMPORTED)
    set_target_properties(IPP::IPP PROPERTIES
        IMPORTED_LOCATION "${IPP_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${IPP_INCLUDE_DIRS}"
    )
endif()

mark_as_advanced(IPP_INCLUDE_DIR IPP_LIBRARIES)