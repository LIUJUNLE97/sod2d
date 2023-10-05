#libs
if (SHARED_LIBS)
    add_subdirectory(lib_sod2d)
    add_subdirectory(lib_sod2d_incomp)
    add_subdirectory(lib_mainBaseClass)
endif()

#apps
add_subdirectory(app_sod2d)