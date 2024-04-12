project(lib_senseiInterface)
# borrowing ideas from https://github.com/LucasGasparinoBSC/CmakeTemplate_Fortran/blob/develop/src/lib1/

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.cpp)
set(HEADER_DIR ${CMAKE_CURRENT_LIST_DIR}/header)
file(GLOB HEADER_FILES ${HEADER_DIR}/*.h)

include_directories(${HEADER_DIR})

add_library(${PROJECT_NAME} ${LIBRARY_TYPE} ${SRC_FILES} ${HEADER_FILES})

##set_mpi()
set_sensei()

set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${TEST_VERSION} SOVERSION ${TEST_VERSION_MAJOR})
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(${PROJECT_NAME}_INCLUDE_DIRS ${HEADER_DIR} CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(FILES ${HEADER_FILES} DESTINATION include)

