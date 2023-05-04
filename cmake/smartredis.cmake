if (USE_SMARTREDIS)
	message("-- Configuring SmartRedis...")
	if(DEFINED ENV{SMARTREDIS_DIR})
	  message("-- SMARTREDIS_DIR properly set in environmental var:" $ENV{SMARTREDIS_DIR})
	  set(SMARTREDIS_DIR $ENV{SMARTREDIS_DIR})
	elseif(SMARTREDIS_DIR)
	  message("-- SMARTREDIS_DIR properly set in cmake:" ${SMARTREDIS_DIR})
	else()
	  message("-- POSSIBLE ERROR! Environmental variables SMARTREDIS_DIR not properly set!")
	  message("-- If Cmake does not find SMARTREDIS package you have two options:")
	  message("--   1. Define environmental variable SMARTREDIS_DIR")
	  message("--   2. Define SMARTREDIS_DIR in cmake/smartredis.cmake")
	  message("-- If Cmake finds SMARTREDIS everything is Ok!")
	endif()
else()
	message("  -- Not using SmartRedis...")
endif()

function(set_smartredis)
	if (USE_SMARTREDIS)
		find_library(SMARTREDIS_LIB smartredis
			PATHS ${SMARTREDIS_DIR}/lib
			NO_DEFAULT_PATH REQUIRED
		)
		if(NOT SMARTREDIS_LIB)
			message(FATAL_ERROR "smartredis library not found in" ${SMARTREDIS_DIR}/lib)
		else()
			message("-- smartredis library found in" ${SMARTREDIS_DIR}/lib)
		endif()

		find_library(SMARTREDIS_FORTRAN_LIB smartredis-fortran
            PATHS ${SMARTREDIS_DIR}/lib
            NO_DEFAULT_PATH REQUIRED
		)
		if(NOT SMARTREDIS_FORTRAN_LIB)
			message(FATAL_ERROR "smartredis-fortran library not found in" ${SMARTREDIS_DIR}/lib)
		else()
			message("-- smartredis-fortran library found in" ${SMARTREDIS_DIR}/lib)
  		endif()

		include_directories(${SMARTREDIS_DIR}/include)
		target_link_libraries(${PROJECT_NAME} ${SMARTREDIS_LIB} ${SMARTREDIS_FORTRAN_LIB})
		add_definitions(-DSMARTREDIS=1)
	else()
		message("  -- Not using SmartRedis...")
		add_definitions(-DSMARTREDIS=0)
	endif()
endfunction()
