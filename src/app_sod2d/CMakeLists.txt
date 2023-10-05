##------------------------------------------------------------------------------------------
## CMakeFile to compile everything as a whole, not using the libs structure, 
## since it is very problematic when enabling acc and cuda flags due to nvlink issues
##
##------------------------------------------------------------------------------------------

project(sod2d)

if(SHARED_LIBS)
   file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90)
else()
   file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/../lib_sod2d/sources/*.f90 
                       ${CMAKE_CURRENT_LIST_DIR}/../lib_sod2d_incomp/sources/*.f90 
                       ${CMAKE_CURRENT_LIST_DIR}/../lib_mainBaseClass/sources/*.f90 
                       ${CMAKE_CURRENT_LIST_DIR}/sources/*.f90)
endif()

add_executable(${PROJECT_NAME} ${SRC_FILES})   # Build executable from listed sources
set_property(TARGET ${PROJECT_NAME} PROPERTY LINKER_LANGUAGE Fortran)

if(SHARED_LIBS)
   target_link_libraries(${PROJECT_NAME} lib_mainBaseClass)
else()
   set_mpi()
   set_hdf5()
endif()

install(TARGETS ${PROJECT_NAME} DESTINATION bin)