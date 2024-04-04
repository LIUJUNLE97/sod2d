message("-- Checking if GEMPA submodule has already been cloned...")
set(ENABLE_TESTS ON CACHE INTERNAL "")
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/gempa/CMakeLists.txt")
  message("  -- GEMPA submodule has not been cloned yet, cloning...")
  include(FetchContent)
  FetchContent_Declare(
    gempa
    GIT_REPOSITORY https://gitlab.com/rickbp/gempa.git
    GIT_TAG        master
  )
  FetchContent_GetProperties(gempa)
  if(NOT gempa_POPULATED)
    FetchContent_Populate(gempa)
    add_subdirectory(${gempa_SOURCE_DIR} ${gempa_BINARY_DIR})
  set(GEMPA_INCLUDE_DIR "${gempa_SOURCE_DIR}/include" CACHE INTERNAL "GEMPA: include directory" FORCE)
  endif()
else()
  message(" -- Adding GEMPA to compilation list...")
  add_subdirectory(gempa)
  set(GEMPA_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gempa/include" CACHE INTERNAL "GEMPA: include directory" FORCE)
endif()

message("-- Checking if JF submodule has already been cloned...")
set(SKIP_DOC_GEN TRUE CACHE INTERNAL "")
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/json-fortran/CMakeLists.txt")
  message("  -- JF submodule has not been cloned yet, cloning...")
  include(FetchContent)
  FetchContent_Declare(
    json-fortran
    GIT_REPOSITORY https://github.com/jacobwilliams/json-fortran.git
    GIT_TAG        master
  )
  FetchContent_GetProperties(json-fortran)
  if(NOT json-fortran_POPULATED)
    FetchContent_Populate(json-fortran)
    add_subdirectory(${json-fortran_SOURCE_DIR} ${json-fortran_BINARY_DIR})
  endif()
else()
  message(" -- Adding JF to compilation list...")
  add_subdirectory(json-fortran)
endif()