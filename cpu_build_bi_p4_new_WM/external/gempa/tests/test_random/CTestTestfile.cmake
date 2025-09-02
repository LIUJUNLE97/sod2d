# CMake generated Testfile for 
# Source directory: /scratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_random
# Build directory: /scratch/junle/CODE_SOD2D/sod2d_gitlab/cpu_build_bi_p4_new_WM/external/gempa/tests/test_random
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_random_TEST "/scratch/muma/libs/nvhpc/Linux_x86_64/24.5/comm_libs/mpi/bin/mpiexec" "-n" "3" "./test_random")
set_tests_properties(test_random_TEST PROPERTIES  WORKING_DIRECTORY "/scratch/junle/CODE_SOD2D/sod2d_gitlab/cpu_build_bi_p4_new_WM/external/gempa/tests/test_random" _BACKTRACE_TRIPLES "/scratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_random/CMakeLists.txt;37;add_test;/scratch/junle/CODE_SOD2D/sod2d_gitlab/external/gempa/tests/test_random/CMakeLists.txt;0;")
