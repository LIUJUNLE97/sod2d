# Declare minimum version of CMake required
cmake_minimum_required(VERSION 3.15)

project(tool_postProc)           # Project name
enable_language(Fortran) # Supported language Fortran

option(USE_PCPOWER "Compiling in PC-POWER" OFF)
option(USE_MN "Compiling in MareNostrum" OFF)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

#----------------COMPILER------------------------------------------------------------------
message("-- Selecting compiler Ops...")
if(USE_PCPOWER)
  message("-- Compiling in PCPOWER...")
  set(CMAKE_C_COMPILER       "/apps/OPENMPI/3.0.0/GCC/bin/mmpicc")
  set(CMAKE_CXX_COMPILER     "/apps/OPENMPI/3.0.0/GCC/bin/mpic++")
  set(CMAKE_Fortran_COMPILER "/apps/OPENMPI/3.0.0/GCC/bin/mpif90")
elseif(USE_MN)
  message("-- Compiling in MN...")
  set(CMAKE_C_COMPILER       "/apps/INTEL/2018.4.057/impi/2018.4.274/bin64/mpicc")
  set(CMAKE_CXX_COMPILER     "/apps/INTEL/2018.4.057/impi/2018.4.274/bin64/mpicxx")
  set(CMAKE_Fortran_COMPILER "/apps/INTEL/2018.4.057/impi/2018.4.274/bin64/mpif90")
else()
  message("-- Compiling in Local...")
  set(CMAKE_C_COMPILER       "/usr/bin/mpicc")
  set(CMAKE_CXX_COMPILER     "/usr/bin/mpic++")
  set(CMAKE_Fortran_COMPILER "/usr/bin/mpif90")
endif()

#----------------COMPILER------------------------------------------------------------------

# Compile
file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_constants.f90 
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_mpi.f90 
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_utils.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_mpi_mesh.f90
   ${CMAKE_CURRENT_LIST_DIR}/../src/lib_sod2d/sources/mod_hdf5.f90
   ${CMAKE_CURRENT_LIST_DIR}/src/*.f90) 

add_executable(${PROJECT_NAME} ${SRC_FILES})   # Build executable from listed sources
set_property(TARGET ${PROJECT_NAME} PROPERTY LINKER_LANGUAGE Fortran)

#----------------MPI------------------------------------------------------------------
message("-- Linking MPI...")
if(USE_PCPOWER)
  link_directories("/apps/OPENMPI/3.0.0/GCC/lib/")
  target_link_libraries(${PROJECT_NAME} mpi_mpifh mpi mpi_cxx)
elseif(USE_MN)

  link_directories("/apps/INTEL/2018.4.057/impi/2018.4.274/lib64")
  target_link_libraries(${PROJECT_NAME} mpifort mpi mpicxx)
else()
  link_directories("/usr/lib/x86_64-linux-gnu/")
  target_link_libraries(${PROJECT_NAME} mpif90 mpi mpicxx)
endif()
#-----------------------------------------------------------------------------------
#----------------HDF5----------------------------------
message("-- Linking HDF5...")

if(USE_PCPOWER)
  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 REQUIRED COMPONENTS Fortran C)
  include_directories(${HDF5_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} hdf5 hdf5_fortran)
elseif(USE_MN)
  find_package(HDF5 REQUIRED COMPONENTS Fortran C)
  include_directories(${HDF5_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} hdf5 hdf5_fortran)
else()
  include_directories("/home/jmuela/libs/hdf5-gcc/include/")
  add_library(hdf5 SHARED IMPORTED)
  set_property(TARGET hdf5 PROPERTY IMPORTED_LOCATION "/home/jmuela/libs/hdf5-gcc/lib/libhdf5.so")
  add_library(hdf5_fortran SHARED IMPORTED)
  set_property(TARGET hdf5_fortran PROPERTY IMPORTED_LOCATION "/home/jmuela/libs/hdf5-gcc/lib/libhdf5_fortran.so")
  #link_directories("/home/jmuela/libs/hdf5-gcc/lib/")
  target_link_libraries(${PROJECT_NAME} hdf5 hdf5_fortran)
endif()

#-----------------------------------------------------------------------------------
#----------------CGNS----------------------------------
message("-- Configuring CGNS...")
if(USE_PCPOWER)
  include_directories("/apps/CGNS/4.3.0/GCC/OPENMPI/include/")
  add_library(cgns SHARED IMPORTED)
  set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/apps/CGNS/4.3.0/GCC/OPENMPI/lib/libcgns.so")
  target_link_libraries(${PROJECT_NAME} cgns)
elseif(USE_MN)
  include_directories( "/apps/CGNS/4.3.0/INTEL/IMPI/include")
  add_library(cgns STATIC IMPORTED)
  set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/apps/CGNS/4.3.0/INTEL/IMPI/lib/libcgns.a")
  target_link_libraries(${PROJECT_NAME} cgns)
else()
  include_directories("/home/jmuela/libs/cgns/include/")
  add_library(cgns STATIC IMPORTED)
  set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/home/jmuela/libs/cgns/lib/libcgns.a")
  target_link_libraries(${PROJECT_NAME} cgns)
endif()
#-----------------------------------------------------------------------------------

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
	"Choose the installation directory; by default it installs in the SOD2D directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

set(CMAKE_Fortran_FLAGS         "-cpp -lstdc++ -DGEMPAINTERFACE -DNOACC")
set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -lstdc++ -DGEMPAINTERFACE -DNOACC")
set(CMAKE_Fortran_FLAGS_DEBUG   "-cpp -lstdc++ -DGEMPAINTERFACE -DNOACC")


