##------------------------------------------------------------------------------------------
## CMakeFile to compile everything as a whole, not using the libs structure, 
## since it is very problematic when enabling acc and cuda flags due to nvlink issues
##
##------------------------------------------------------------------------------------------

project(sod2d)

file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/../lib_sod2d/sources/*.f90 
   ${CMAKE_CURRENT_LIST_DIR}/../lib_mainBaseClass/sources/*.f90 
   ${CMAKE_CURRENT_LIST_DIR}/sources/*.f90) 
add_executable(${PROJECT_NAME} ${SRC_FILES})   # Build executable from listed sources
set_property(TARGET ${PROJECT_NAME} PROPERTY LINKER_LANGUAGE Fortran)

#add_definitions(-D__BUILD_VERSION__="${TEST_VERSION}")
#add_executable(${PROJECT_NAME} ${SRC_FILES})

set_mpi()
set_hdf5()

#set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME} VERSION ${TEST_VERSION} SOVERSION ${TEST_VERSION_MAJOR})

#install(TARGETS ${PROJECT_NAME} DESTINATION bin)