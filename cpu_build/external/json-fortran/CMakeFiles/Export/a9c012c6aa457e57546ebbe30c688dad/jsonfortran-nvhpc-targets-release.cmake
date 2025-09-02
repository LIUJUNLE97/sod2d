#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "jsonfortran-nvhpc::jsonfortran" for configuration "Release"
set_property(TARGET jsonfortran-nvhpc::jsonfortran APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(jsonfortran-nvhpc::jsonfortran PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so.8.3.0"
  IMPORTED_SONAME_RELEASE "libjsonfortran.so.8.3"
  )

list(APPEND _cmake_import_check_targets jsonfortran-nvhpc::jsonfortran )
list(APPEND _cmake_import_check_files_for_jsonfortran-nvhpc::jsonfortran "${_IMPORT_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.so.8.3.0" )

# Import target "jsonfortran-nvhpc::jsonfortran-static" for configuration "Release"
set_property(TARGET jsonfortran-nvhpc::jsonfortran-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(jsonfortran-nvhpc::jsonfortran-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.a"
  )

list(APPEND _cmake_import_check_targets jsonfortran-nvhpc::jsonfortran-static )
list(APPEND _cmake_import_check_files_for_jsonfortran-nvhpc::jsonfortran-static "${_IMPORT_PREFIX}/jsonfortran-nvhpc-8.3.0/lib/libjsonfortran.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
