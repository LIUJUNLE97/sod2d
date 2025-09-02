# Install script for directory: /lscratch/junle/CODE_SOD2D/sod2d_gitlab/external/json-fortran

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

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so.8.3.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so.8.3"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib" TYPE SHARED_LIBRARY FILES
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/lib/libjsonfortran.so.8.3.0"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/lib/libjsonfortran.so.8.3"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so.8.3.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so.8.3"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib" TYPE SHARED_LIBRARY FILES "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/lib/libjsonfortran.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib" TYPE STATIC_LIBRARY FILES "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/lib/libjsonfortran.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(GLOB_RECURSE MODULE_FILES "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/include/*.mod")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(GLOB_RECURSE SUBMOD_FILES "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/include/*.smod")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL ${MODULE_FILES} DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL ${SUBMOD_FILES} DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/cmake/jsonfortran-nvhpc-targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/cmake/jsonfortran-nvhpc-targets.cmake"
         "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/CMakeFiles/Export/a9c012c6aa457e57546ebbe30c688dad/jsonfortran-nvhpc-targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/cmake/jsonfortran-nvhpc-targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/cmake/jsonfortran-nvhpc-targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/cmake" TYPE FILE FILES "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/CMakeFiles/Export/a9c012c6aa457e57546ebbe30c688dad/jsonfortran-nvhpc-targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/cmake" TYPE FILE FILES "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/CMakeFiles/Export/a9c012c6aa457e57546ebbe30c688dad/jsonfortran-nvhpc-targets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/cmake" TYPE FILE FILES
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/pkg/jsonfortran-nvhpc-config.cmake"
    "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/jsonfortran-nvhpc-config-version.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/pkgconfig" TYPE FILE FILES "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/json-fortran/json-fortran.pc")
endif()

