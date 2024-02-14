function(set_json-fortran)
  include_directories(${GEMPA_INCLUDE_DIR})
  target_link_libraries(${PROJECT_NAME} jsonfortran)
endfunction()