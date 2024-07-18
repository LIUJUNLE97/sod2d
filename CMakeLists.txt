# Set minimum CMake version
cmake_minimum_required(VERSION 3.15)

# Set project versions: major, minor and patch
set(TEST_VERSION_MAJOR "0")
set(TEST_VERSION_MINOR "0")
set(TEST_VERSION_PATCH "0")
set(TEST_VERSION "${TEST_VERSION_MAJOR}.${TEST_VERSION_MINOR}.${TEST_VERSION_PATCH}")

# Options
option(UNIT_TESTING "Enable unit tests" OFF)
option(USE_MPI "Compile using MPI" ON)
option(USE_HDF5 "Compile using HDF5" ON)
option(USE_SMARTREDIS "Compile using SmartRedis" ON)

option(TOOL_MESHPART "Compiles tool MeshPart" OFF)
option(TOOL_COPYRES  "Compiles tool CopyResults" OFF)
option(TOOL_COMMPERF "Compiles tool CommPerf" OFF)

option(USE_PREDEFINED_COMPILERS "Use predefined compilers" OFF)
option(USE_PCPOWER "Compiling in PC-POWER" OFF)
option(USE_ALVIS "Compiling in Alvis" OFF)
option(USE_BERZELIUS "Compiling in Berzelius" OFF)
option(USE_MN4 "Compiling in MareNostrum" OFF)
option(USE_MN5 "Compiling in MareNostrum5" OFF)
option(USE_CRAY "Compiling with CrayMpi" OFF)
option(SHARED_LIBS "Compile using shared libraries" OFF)

option(USE_GPU "Compile using GPU" ON)
option(USE_MEM_MANAGED "Compile using Mem managed" OFF)
option(USE_NCCL "Compile using NCCL" OFF)

option(USE_SENSEI "Compile using SENSEI" OFF)

option(DEBUG_MODE "Set debug mode" OFF)

# Folder with files configuring extra CMake options
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Set MPI Compilers if required
if(USE_CRAY)
    message("Compiling with CrayMpi compilers")
    set(ENV{CC} "cc")
    set(ENV{CXX} "CC")
    set(ENV{FC} "ftn")
else()
    set(ENV{CC} "mpicc")
    set(ENV{CXX} "mpicxx")
    set(ENV{FC} "mpif90")
endif()

# Set project name, version and languages
project(sod2d LANGUAGES C CXX Fortran)

# Configure MPI module
include(mpi)

# Configure compiler options
include(compilerOps)

# Configure HDF5 module
include(hdf5)

# Include SmartRedis
include(smartredis)

# Configure GEMPA
include(gempa)

# Configure JFORTRAN module
include(json-fortran)

# Configure SENSEI
include(sensei)

# Include external projects
#include(ExternalProject) (I think this does nothing)
add_subdirectory(external)

if(DEBUG_MODE)
    set(CMAKE_BUILD_TYPE Debug)
else()
    set(CMAKE_BUILD_TYPE Release)
endif()

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

add_subdirectory(src)

if(TOOL_MESHPART)
    #include(gempa)
    #add_subdirectory(external)
    add_subdirectory(tool_meshConversorPar)
endif()

if(TOOL_COPYRES)
    add_subdirectory(tool_copyResults)
endif()

if(TOOL_COMMPERF)
    add_subdirectory(tool_commsPerformance)
endif()

#enable folder unitt
#add_subdirectory(unitt)

# Enable testing
if(UNIT_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
