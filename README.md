# SOD2D

**SOD2D**: **S**pectral high-**O**rder co**D**e **2** solve partial **D**ifferential equations

## Description

Before starting: we have a Wiki, please take a look before proceeding :D

This code implements a numerical solution for the equations governing compressible and incompressible fluid flow in three dimensions. The code is based on the spectral element method (SEM) and is designed to be used for scale-resolving simulations (LES and DNS).

It is written in Fortran, and uses MPI and OpenACC to provide parallelism at both coarse and fine-grained levels. The mesh is partitioned using the GEMPA library (https://gitlab.com/rickbp/gempa), which is included as a submodule in this repository. It also uses HDF5 for I/O, which must be installed appropriately according to the desired platform (see DEPENDENCIES section).

## Dependencies

The code depends on the following libraries:

* HDF5
* MPI
* OpenACC (optional)
* GEMPA (included as a submodule)
* CMake (building)

The HDF5 library must be properly compiled with the same compiler that will be used for the code. For example, if the code is compiled with the PGI compiler, then HDF5 must be compiled with the PGI compiler.

## Building

The code is built using CMake. The following commands will build the code in the `build` directory:

```bash
mkdir build
cd build
cmake ..
make
```

Selecting a non-default compiler requires setting the following options to CMake:

```bash
cmake -DCMAKE_Fortran_COMPILER=[FC] -DCMAKE_C_COMPILER=[CC] -DCMAKE_CXX_COMPILER=[CXX] ..
```

Building with the PGI compiler requires that the NVIDIA HPC SDK library be properly installed, including the MPI library packaged in the SDK.

At a minimum, some form of OpenMPI must be installed in the system. The code will not build without MPI.

## Running

The code is run using the `sod2d` executable, found in the `build` directory. The executable takes no arguments. Configuration is done using the implemented classes, found inside the `src/classes` directory. The `src/main.f90` file contains the main program, which is responsible for initializing the classes and running the simulation.

## Acknowledgements

The research leading to this software has received funding from the European High-Performance Computing Joint Undertaking (JU) under grant agreement No 956104. The JU receives support from the European Union’s Horizon 2020 research and innovation programme and Spain, France, Germany.
O. Lehmkuhl work is financed by a Ramón y Cajal postdoctoral contract by the Ministerio de Economía y Competitividad, Secretaría de Estado de Investigación, Desarrollo e Innovación, Spain (RYC2018-025949-I).


