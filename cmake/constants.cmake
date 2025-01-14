message("-- <${PROJECT_NAME}> constants configuration:")

# First checks
#------------------------------------------------------------------------------------------
if(USE_PORDER LESS 2)
    message(FATAL_ERROR "USE_PORDER must be 2 or higher. Current value: ${USE_PORDER}")
endif()
if(NOT (USE_RP STREQUAL "4" OR USE_RP STREQUAL "8"))
    message(FATAL_ERROR "USE_RP must be either 4 or 8. Current value: ${USE_RP}")
endif()
if(NOT (USE_RP_VTK STREQUAL "4" OR USE_RP_VTK STREQUAL "8"))
    message(FATAL_ERROR "USE_RP_VTK must be either 4 or 8. Current value: ${USE_RP_VTK}")
endif()
if(NOT (USE_RP_AVG STREQUAL "4" OR USE_RP_AVG STREQUAL "8"))
    message(FATAL_ERROR "USE_RP_AVG must be either 4 or 8. Current value: ${USE_RP_AVG}")
endif()
#------------------------------------------------------------------------------------------

# Give some user info
#------------------------------------------------------------------------------------------
message("  -- Using porder: " ${USE_PORDER})
message("  -- Using real precision: " ${USE_RP})
message("  -- Using vtk precision: " ${USE_RP_VTK})
message("  -- Using average precision: " ${USE_RP_AVG})

# Add definitions for porder and precision
#------------------------------------------------------------------------------------------
add_definitions(-D__PORDER__=${USE_PORDER})
add_definitions(-D__RP__=${USE_RP})
add_definitions(-D__RP_VTK__=${USE_RP_VTK})
add_definitions(-D__RP_AVG__=${USE_RP_AVG})