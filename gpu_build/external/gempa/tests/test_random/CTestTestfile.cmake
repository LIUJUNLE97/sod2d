# CMake generated Testfile for 
# Source directory: /lscratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_random
# Build directory: /lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/gempa/tests/test_random
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_random_TEST "/scratch/muma/libs/nvhpc/Linux_x86_64/24.5/comm_libs/mpi/bin/mpiexec" "-n" "3" "./test_random")
set_tests_properties(test_random_TEST PROPERTIES  WORKING_DIRECTORY "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/gpu_build/external/gempa/tests/test_random" _BACKTRACE_TRIPLES "/lscratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_random/CMakeLists.txt;37;add_test;/lscratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_random/CMakeLists.txt;0;")
