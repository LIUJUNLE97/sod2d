function(set_gempa)
  include_directories(${GEMPA_INCLUDE_DIR})
  target_link_libraries(${PROJECT_NAME} gempa)
endfunction()
