# CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(pdb2mrc VERSION 1.0.0 LANGUAGES CXX C)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Configuration options
option(PDB2MRC_BUILD_SHARED "Build shared libraries" OFF)
option(PDB2MRC_ENABLE_OPENMP "Enable OpenMP parallelization" ON)
option(PDB2MRC_ENABLE_TESTS "Build tests" OFF)
option(PDB2MRC_ENABLE_EXAMPLES "Build examples" OFF)
option(PDB2MRC_VERBOSE_CONFIG "Print verbose configuration info" OFF)

# Set build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Include CMake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Find required packages
find_package(OpenMP REQUIRED)
find_package(MKL REQUIRED)
find_package(IPP REQUIRED)

# Print configuration
if(PDB2MRC_VERBOSE_CONFIG)
    message(STATUS "=== pdb2mrc Configuration ===")
    message(STATUS "Version: ${PROJECT_VERSION}")
    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")
    message(STATUS "MKL found: ${MKL_FOUND}")
    message(STATUS "IPP found: ${IPP_FOUND}")
    message(STATUS "==============================")
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${IPP_INCLUDE_DIRS})
include_directories(${MKL_INCLUDE_DIRS})

# Collect source files
file(GLOB LIB_SOURCES
    "src/*.cpp"
)

# Exclude main.cpp from library sources
list(REMOVE_ITEM LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

# Create library
if(PDB2MRC_BUILD_SHARED)
    add_library(pdb2mrc SHARED ${LIB_SOURCES})
else()
    add_library(pdb2mrc STATIC ${LIB_SOURCES})
endif()

# Set library properties
set_target_properties(pdb2mrc PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    OUTPUT_NAME pdb2mrc
    PUBLIC_HEADER "include/atom.hpp;include/bfactor_utils.hpp;include/chimerax_generator.hpp;include/emmer_generator.hpp;include/enums.hpp;include/grid_projection.hpp;include/main_controller.hpp;include/map_generator.hpp;include/mrc_writer.hpp;include/pdb_reader.hpp;include/platform.h;include/scattering.hpp;include/situs_generator.hpp;include/types.hpp"
)

# Link libraries
target_link_libraries(pdb2mrc
    PUBLIC
        ${IPP_LIBRARIES}
        ${MKL_LIBRARIES}
        OpenMP::OpenMP_CXX
    PRIVATE
        ${CMAKE_DL_LIBS}
)

# Create executable
add_executable(pdb2mrc_exe src/main.cpp)
target_link_libraries(pdb2mrc_exe PRIVATE pdb2mrc)

# Rename executable to pdb2mrc
set_target_properties(pdb2mrc_exe PROPERTIES OUTPUT_NAME pdb2mrc)

# Installation
install(TARGETS pdb2mrc pdb2mrc_exe
    EXPORT pdb2mrcTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include/pdb2mrc
)

install(DIRECTORY include/ DESTINATION include/pdb2mrc
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)

# Export package configuration
install(EXPORT pdb2mrcTargets
    FILE pdb2mrcTargets.cmake
    NAMESPACE pdb2mrc::
    DESTINATION lib/cmake/pdb2mrc
)

# Generate package config file
include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/pdb2mrcConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/pdb2mrcConfig.cmake
    INSTALL_DESTINATION lib/cmake/pdb2mrc
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/pdb2mrcConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/pdb2mrcConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/pdb2mrcConfigVersion.cmake
    DESTINATION lib/cmake/pdb2mrc
)

# Optionally build tests
if(PDB2MRC_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Optionally build examples
if(PDB2MRC_ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif()