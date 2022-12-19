# Declare minimum version of CMake required
cmake_minimum_required(VERSION 3.15)

#----------------OPTIONS------------------------------------------------------------------

message("Verify that machine selection is identical to the one in the main CMakeLists.txt file!")
option(USE_PCPOWER "Compiling in PC-POWER" OFF)
option(USE_MN "Compiling in MareNostrum" OFF)
option(USE_PREDEFINED_COMPILERS "Use predefined compilers" ON)

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
project(tool_meshConversorPar LANGUAGES C CXX Fortran)

# Configure MPI module
include(mpi)

# Configure HDF5 module
include(hdf5)

# Set a default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

#----------------COMPILE------------------------------------------------------------------
# Compile
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_constants.f90 
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_mpi.f90 
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_utils.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_gmsh_indices.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_mpi_mesh.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_comms.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_mesh_boundaries.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_comms_boundaries.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_hdf5.f90
   ${CMAKE_CURRENT_LIST_DIR}/src/mod_meshConversorTool.f90
   ${CMAKE_CURRENT_LIST_DIR}/src/tool_meshConversorPar.f90) 

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

#set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} "-DGEMPAINTERFACE")
#string(REPLACE ";" " " CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")