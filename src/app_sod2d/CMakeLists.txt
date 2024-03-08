project(sod2d)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

if(SHARED_LIBS)
   file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90)
else()
   file(GLOB SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/../lib_sod2d/sources/*.f90 
                       ${CMAKE_CURRENT_LIST_DIR}/../lib_sod2d_incomp/sources/*.f90 
                       ${CMAKE_CURRENT_LIST_DIR}/../lib_mainBaseClass/sources/*.f90 
                       ${CMAKE_CURRENT_LIST_DIR}/sources/*.f90)
endif()
include_directories(${CMAKE_BINARY_DIR}/modules)

set(CPP_SRC_FILES 
${CMAKE_CURRENT_LIST_DIR}/../lib_sod2d/sources/senseiCInterface.h
${CMAKE_CURRENT_LIST_DIR}/../lib_sod2d/sources/senseiCInterface.cpp
${CMAKE_CURRENT_LIST_DIR}/../lib_sod2d/sources/senseiDataAdaptor.h
${CMAKE_CURRENT_LIST_DIR}/../lib_sod2d/sources/senseiDataAdaptor.cpp)

add_executable(${PROJECT_NAME} ${SRC_FILES} ${CPP_SRC_FILES})   # Build executable from listed sources
set_property(TARGET ${PROJECT_NAME} PROPERTY LINKER_LANGUAGE Fortran)

if(SHARED_LIBS)
   target_link_libraries(${PROJECT_NAME} lib_mainBaseClass)
else()
   set_mpi()
   set_hdf5()
   set_jsonfortran()
endif()

find_package(SENSEI REQUIRED)

find_file(SENSEI_CONFIG_FILE DOC "config for SENSEI's configuarable analysis adaptor" REQUIRED)
target_compile_definitions(${PROJECT_NAME} PRIVATE SENSEI_CONFIG_FILE="${SENSEI_CONFIG_FILE}")
target_link_libraries(${PROJECT_NAME} sensei)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION include)