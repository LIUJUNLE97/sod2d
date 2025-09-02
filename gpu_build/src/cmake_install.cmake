# Install script for directory: /lscratch/junle/CODE_SOD2D/sod2d_gitlab/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/src/app_sod2d/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/config_files" TYPE FILE FILES
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/ABlFlowSolverIncomp.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/BLFlowSolverIncomp.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/BLFlowSolverIncompAFC.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/BLFlowSolverIncompDRL.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/BluffBody3DSolver.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/BluffBody3DSolverIncomp.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/BluffBodySolver.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/BluffBodySolverIncomp.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/BluffBodySolverIncompAFC.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/BluffBodySolverIncompDRL.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/ChannelFlowSolver.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/ChannelFlowSolverIncomp.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/MappedInletIncomp.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/SupersonicForwardStep.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/TGVCompSolver.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/TGVSolver.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/TGVSolverIncomp.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/ThermalBubbleSolver.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/ThermalChannelFlowSolver.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/TransientInletSolverIncomp.json"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/src/lib_mainBaseClass/config_files/WindFarmSolverIncomp.json"
    )
endif()

