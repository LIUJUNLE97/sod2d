# Declare minimum version of CMake required
cmake_minimum_required(VERSION 3.15)

#----------------OPTIONS------------------------------------------------------------------

option(USE_PCPOWER "Compiling in PC-POWER" OFF)
option(USE_MN "Compiling in MareNostrum" OFF)
option(USE_PREDEFINED_COMPILERS "Use predefined compilers" ON)

#option(USE_GEMPA "Compile using GEMPA" ON)
option(USE_MPI "Compile using MPI" ON)
option(USE_HDF5 "Compile using HDF5" ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Configure compiler options
message("-- Selecting compiler Ops...")
if(USE_PCPOWER)
  include(powerpcCompilers)
elseif(USE_MN)
  include(marenostrumCompilers)
else()
  include(localCompilers)
endif()

# Set project name, version and languages
project(tool_commsPerformance LANGUAGES C CXX Fortran)

# Configure MPI module
include(mpi)

# Configure HDF5 module
include(hdf5)

# Include external projects
#include(ExternalProject)
# Enable GEMPA if requested
#include(gempa)
#add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../external ${CMAKE_CURRENT_BINARY_DIR}/external/ )

# Set a default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#----------------COMPILE------------------------------------------------------------------
# Compile
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_constants.f90 
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_mpi.f90 
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_utils.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_custom_types.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_mpi_mesh.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_hdf5.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_comms.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_comms_boundaries.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_nvtx.f90
   ${CMAKE_CURRENT_LIST_DIR}/sources/mod_comms_performance.f90
   ${CMAKE_CURRENT_LIST_DIR}/sources/tool_comms_performance.f90) 
   
add_executable(${PROJECT_NAME} ${SRC_FILES})   # Build executable from listed sources
set_property(TARGET ${PROJECT_NAME} PROPERTY LINKER_LANGUAGE Fortran)

set_mpi()
set_gempa()
set_hdf5()

target_link_libraries(${PROJECT_NAME} lib_gempaInterface)

#-----------------------------------------------------------------------------------

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
	"Choose the installation directory; by default it installs in the SOD2D directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


#set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wsign-conversion -std=c++11 -Wold-style-cast")
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++11")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#Enable C++ 11
#set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_Fortran_FLAGS         "-cpp -lstdc++ -DNOACC")
#set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -lstdc++ -DNOACC")
#set(CMAKE_Fortran_FLAGS_DEBUG   "-cpp -lstdc++ -DNOACC")

#------OLD-----------------------
#project(testComms)
#file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90)
#include_directories(${lib_sod2d_INCLUDE_DIRS})
#add_definitions(-D__BUILD_VERSION__="${TEST_VERSION}")
#add_executable(${PROJECT_NAME} ${SRC_FILES})
#set_mpi()
#set_gempa()
#set_hdf5()
#target_link_libraries(${PROJECT_NAME} lib_sod2d)
#set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME} VERSION ${TEST_VERSION} SOVERSION ${TEST_VERSION_MAJOR})
#install(TARGETS ${PROJECT_NAME} DESTINATION bin)
