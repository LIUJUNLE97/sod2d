# Install script for directory: /home/lucas/sod2d_github/unitt

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/lucas")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
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
  # Include the install script for each subdirectory.
  include("/home/lucas/sod2d_github/build/unitt/unitt_mass_convec/cmake_install.cmake")
  include("/home/lucas/sod2d_github/build/unitt/unitt_mom_convec/cmake_install.cmake")
  include("/home/lucas/sod2d_github/build/unitt/unitt_mom_diffu/cmake_install.cmake")
  include("/home/lucas/sod2d_github/build/unitt/unitt_ener_convec/cmake_install.cmake")
  include("/home/lucas/sod2d_github/build/unitt/unitt_generic_scalar_convec/cmake_install.cmake")
  include("/home/lucas/sod2d_github/build/unitt/unitt_3d_jacobian_inverse/cmake_install.cmake")
  include("/home/lucas/sod2d_github/build/unitt/unitt_char_length/cmake_install.cmake")

endif()

