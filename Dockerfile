### Basic OS image that allows SOD2D to be built and tested

## Import the base image
FROM ubuntu:22.04

## Update and install basic system packages
RUN apt-get -y update && apt-get -y upgrade
RUN apt-get install -y build-essential git cmake gfortran ninja-build wget

## Set the working dir to home
WORKDIR /home

## Create directory for installing additional libraries
RUN mkdir -p apps/libraries
RUN mkdir -p apps/compilers

## Download and install OpenMPI-4.1.4
WORKDIR /home/apps/compilers
RUN mkdir -p openMPI/4.1.4
WORKDIR /home/apps/compilers/openMPI/4.1.4
RUN wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.4.tar.gz
RUN tar -xvzf openmpi-4.1.4.tar.gz
WORKDIR /home/apps/compilers/openMPI/4.1.4/openmpi-4.1.4
RUN ./configure --prefix=/apps/compilers/openMPI/4.1.4/gnu
RUN make -j 12 && make install
# Set the PATH and LD_LIBRARY_PATH environment variables
ENV PATH="/apps/compilers/openMPI/4.1.4/gnu/bin:${PATH}"
ENV LD_LIBRARY_PATH="/apps/compilers/openMPI/4.1.4/gnu/lib:${LD_LIBRARY_PATH}"

## Download and install HDF5-1.12.0 enabling parallel and fortran support
WORKDIR /home/apps/libraries
RUN mkdir -p hdf5/1.12.0
WORKDIR /home/apps/libraries/hdf5/1.12.0
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz
RUN tar -xvzf hdf5-1.12.0.tar.gz
WORKDIR /home/apps/libraries/hdf5/1.12.0/hdf5-1.12.0
RUN CPP=cpp CC=mpicc CXX=mpicxx FC=mpif90 ./configure --prefix=/apps/libraries/hdf5/1.12.0/gnu --enable-threadsafe --enable-cxx --enable-fortran --enable-unsupported --enable-parallel
RUN make -j 12 && make install
# Set the PATH and LD_LIBRARY_PATH environment variables
ENV PATH="/apps/libraries/hdf5/1.12.0/gnu/bin:${PATH}"
ENV LD_LIBRARY_PATH="/apps/libraries/hdf5/1.12.0/gnu/lib:${LD_LIBRARY_PATH}"
ENV HDF5_ROOT="/apps/libraries/hdf5/1.12.0/gnu"
ENV HDF5_DIR="/apps/libraries/hdf5/1.12.0/gnu"