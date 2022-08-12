if (USE_GEMPA)
  message("-- Check is GEMPA submodule has already been cloned...")
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
endif()