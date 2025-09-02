# CMake generated Testfile for 
# Source directory: /scratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_small_parallel
# Build directory: /scratch/junle/CODE_SOD2D/sod2d_gitlab/cpu_build_barce_p6/external/gempa/tests/test_small_parallel
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_small_parallel_TEST "/scratch/muma/libs/nvhpc/Linux_x86_64/24.5/comm_libs/mpi/bin/mpiexec" "-n" "3" "./test_small_parallel")
set_tests_properties(test_small_parallel_TEST PROPERTIES  WORKING_DIRECTORY "/scratch/junle/CODE_SOD2D/sod2d_gitlab/cpu_build_barce_p6/external/gempa/tests/test_small_parallel" _BACKTRACE_TRIPLES "/scratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_small_parallel/CMakeLists.txt;37;add_test;/scratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_small_parallel/CMakeLists.txt;0;")
add_test(test_small_parallel_DIFF "/usr/bin/cmake" "-E" "compare_files" "partition.txt" "/scratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_small_parallel/output/partition.txt")
set_tests_properties(test_small_parallel_DIFF PROPERTIES  DEPENDS "test_small_parallel_TEST" WORKING_DIRECTORY "/scratch/junle/CODE_SOD2D/sod2d_gitlab/cpu_build_barce_p6/external/gempa/tests/test_small_parallel" _BACKTRACE_TRIPLES "/scratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_small_parallel/CMakeLists.txt;41;add_test;/scratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_small_parallel/CMakeLists.txt;0;")
