# Compressible TGV tests

Set of small TGV cases for testing basic functionality of SOD2D, using the compressible configuration. References obtained using the meshes and configurations contained in the folders (*.hdf and *.json files), using the following hardware:

- Intel(R) Core(TM) i9-8950HK CPU @ 4.80GHz
- NVIDIA(R) Quadro P4200

NOTE: in case a bug is fixed in the code, or some alteration is made that drastically changes the behaviour of the code, the reference results may not be valid anymore. In this case, the reference results should be updated, otherwise a merge will not be possible!

## GPU tests

Uses a larger mesh to account for the GPUs higher performance. Both the IMEX and LS-RK 5 stage configurations are tested under the flow conditions described on the *.json file. Given the hardware limitations, only a single GPU is used to process the case. To obtain the reference results, the code was compiled with NVIDIA's NVHPC-SDK version 24.11, in HPCX mode.

## CPU tests

Similar idea to the GPU tests, except on a smaller mesh. However, both tests are run using 4 CPU cores, so tests account for parallelism in the code. The compiler set used for these configurations was:

- GNU v14.2.1
- OpenMPI v5.0.5