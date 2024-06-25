# Use NVIDIA-SMI tool to automatically find the
# compute capability of the GPU
function(set_cc)
	# Obtain the git hash and store it to a variable
	execute_process(
		COMMAND nvidia-smi --query-gpu=compute_cap --format=csv
		COMMAND tail +2
		COMMAND tr -d .
		OUTPUT_VARIABLE GPU_CC
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	message("GPU compute capability detected: cc" ${GPU_CC})
	set(GPU_CC ${GPU_CC} PARENT_SCOPE)
endfunction()