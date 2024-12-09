function(set_git_version)
	# Obtain the git hash and store it to a variable
	execute_process(
		COMMAND git log -1 --pretty=format:%H 
		COMMAND cut -c1-7
		OUTPUT_VARIABLE GITHASH
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
    message("-- By <" $ENV{USER} "> at git commit <" ${GITHASH} ">")
#    set(GITHASH "${GITHASH}" PARENT_SCOPE)
	# Add the githash definition
#	add_definitions(-DGITHASH="${GITHASH}")
	# Add username definition
#	add_definitions(-DUSER="$ENV{USER}")
endfunction()