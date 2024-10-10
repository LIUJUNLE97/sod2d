#libs
if (SHARED_LIBS)
    add_subdirectory(lib_sod2d)
    add_subdirectory(lib_sod2d_incomp)
    add_subdirectory(lib_species)
    add_subdirectory(lib_mainBaseClass)
endif()

#apps
add_subdirectory(app_sod2d)

# Find and install JSON files
file(GLOB_RECURSE CONFIG_FILES  ${CMAKE_CURRENT_LIST_DIR}/lib_mainBaseClass/config_files/*.json)
install(FILES   ${CONFIG_FILES} DESTINATION share/config_files)
