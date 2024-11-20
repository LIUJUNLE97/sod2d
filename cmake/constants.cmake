message("-- SOD2D constants configuration:")

# Give some user info
message("  -- Using porder: " ${USE_PORDER})
message("  -- Using real precision: " ${USE_RP})
message("  -- Using vtk precision: " ${USE_RP_VTK})
message("  -- Using average precision: " ${USE_RP_AVG})

# Add definitions for porder and precision
add_definitions(-D__PORDER__=${USE_PORDER})
add_definitions(-D__RP__=${USE_RP})
add_definitions(-D__RP_VTK__=${USE_RP_VTK})
add_definitions(-D__RP_AVG__=${USE_RP_AVG})