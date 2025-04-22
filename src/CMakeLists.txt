# Configure CONSTANTS
include(constants)

#libs
if (SHARED_LIBS)
    add_subdirectory(lib_sod2d)
    add_subdirectory(lib_meshElasticity)
    add_subdirectory(lib_sod2d_incomp)
    add_subdirectory(lib_species)
    add_subdirectory(lib_mainBaseClass)
endif()

#apps
add_subdirectory(app_sod2d)

# Find and install JSON files
file(GLOB_RECURSE CONFIG_FILES  ${CMAKE_CURRENT_LIST_DIR}/lib_mainBaseClass/config_files/*.json)
add_custom_target(sod2d-config ALL COMMAND cp -r ${CMAKE_CURRENT_LIST_DIR}/lib_mainBaseClass/config_files ${CMAKE_CURRENT_BINARY_DIR})
install(FILES ${CONFIG_FILES} DESTINATION share/config_files)
