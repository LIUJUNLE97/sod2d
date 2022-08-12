project(testComms)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90)

include_directories(${lib_sod2d_INCLUDE_DIRS})

add_definitions(-D__BUILD_VERSION__="${TEST_VERSION}")

add_executable(${PROJECT_NAME} ${SRC_FILES})

set_mpi()
set_gempa()
set_hdf5()

target_link_libraries(${PROJECT_NAME} lib_sod2d)

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME} VERSION ${TEST_VERSION} SOVERSION ${TEST_VERSION_MAJOR})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)