# Find and install bash files
file(GLOB_RECURSE BASH_FILES  ${CMAKE_CURRENT_LIST_DIR}/*.sh)
install(FILES   ${BASH_FILES} DESTINATION share/utils)

# Define custom targets for the python tools within utils
# so that they are given some love by cmake and managed during installation

# ANSA2SOD2D
add_custom_target(ansa2sod2d ALL 
    COMMAND cp ${CMAKE_CURRENT_LIST_DIR}/ansa2sod2d/ansa2sod2d.py ${CMAKE_CURRENT_BINARY_DIR}/ansa2sod2d.py
    COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/ansa2sod2d.py
)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ansa2sod2d.py RENAME ansa2sod2d DESTINATION bin)

# GMSH2SOD2D
add_custom_target(gmsh2sod2d ALL 
    COMMAND cp ${CMAKE_CURRENT_LIST_DIR}/gmsh2sod2d/gmsh2sod2d.py ${CMAKE_CURRENT_BINARY_DIR}/gmsh2sod2d.py
    COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/gmsh2sod2d.py
)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/gmsh2sod2d.py RENAME gmsh2sod2d DESTINATION bin)

# INTERPOLATE
add_custom_target(interpolate ALL 
    COMMAND cp ${CMAKE_CURRENT_LIST_DIR}/interpolate/interpolate.py ${CMAKE_CURRENT_BINARY_DIR}/interpolate.py
    COMMAND chmod +x ${CMAKE_CURRENT_BINARY_DIR}/interpolate.py
)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/interpolate.py RENAME interpolate DESTINATION bin)