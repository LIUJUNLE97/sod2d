# Set minimum CMake version
cmake_minimum_required(VERSION 3.15)

# Set project versions: major, minor and patch
set(TEST_VERSION_MAJOR "0")
set(TEST_VERSION_MINOR "0")
set(TEST_VERSION_PATCH "0")
set(TEST_VERSION "${TEST_VERSION_MAJOR}.${TEST_VERSION_MINOR}.${TEST_VERSION_PATCH}")

# Options
#option(UNIT_TESTING "Enable unit tests" ON)
option(USE_MPI "Compile using MPI" ON)
option(USE_HDF5 "Compile using HDF5" ON)
option(USE_GEMPA "Compile using GEMPA" ON)

option(USE_PREDEFINED_COMPILERS "Use predefined compilers" ON)
option(USE_PCPOWER "Compiling in PC-POWER" OFF)
option(USE_MN "Compiling in MareNostrum" OFF)

option(USE_GPU "Compile using GPU" OFF)
option(USE_MEM_MANAGED "Compile using Mem managed" OFF)

# Folder with files configuring extra CMake options
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Configure compiler options
include(compilerOps)

# Set project name, version and languages
project(sod2d LANGUAGES C CXX Fortran)

# Configure MPI module
include(mpi)

# Configure HDF5 module
include(hdf5)

# Include external projects
include(ExternalProject)

# Set library types to be shared
set(LIBRARY_TYPE SHARED)
#set(LIBRARY_TYPE STATIC)

# Set a default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set compilation config for certain compilers
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wsign-conversion -Wold-style-cast")
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

# Set C/C++ standard
set(CMAKE_CXX_STANDARD 11)

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# Enable GEMPA if requested
include(gempa)

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(tool_postProc)

#enable folder unitt
#add_subdirectory(unitt)

# Enable testing
#if(UNIT_TESTING)
#    enable_testing()
#    add_subdirectory(tests)
#endif()
