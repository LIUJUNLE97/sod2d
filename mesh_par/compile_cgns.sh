#!/bin/bash

#----------------------------------------------------------------------------------
PROGRAM=read_and_create_mesh
PROGRAM2=load_cgnsmesh
GEMPA_INT=gempa_interface
COMPILER_F=mpif90
FLAGS_F="-cpp -lstdc++"
COMPILER_CPP=mpic++
FLAGS_CPP="-O3"
MODS_PATH=""
#----------------------------------------------------------------------------------
#LIBS=$LIBS" -L/home/jmuela/libs/gempa/lib -lgempa"
#INCLUDES=$INCLUDES" -I/home/jmuela/libs/gempa/include"

LIBS="-L/usr/lib -llapack -L/usr/lib -lblas -L$HOME/libs/cgns/lib -lcgns -L$HOME/libs/hdf5/lib -lhdf5"
INCLUDES="-I$HOME/libs/cgns/include/ -I$HOME/libs/hdf5/include/ -I/usr/local/include"

LIBS=$LIBS" -L/home/jmuela/libs/gempa/lib -lgempa"
INCLUDES=$INCLUDES" -I/home/jmuela/libs/gempa/include"

echo 'Compiling Gempa interface...'
$COMPILER_CPP -c $GEMPA_INT'.cpp' $FLAGS_CPP $LIBS $INCLUDES
echo 'GEMPA Interface Compiled!'
#----------------------------------------------------------------------------------

echo 'Compiling mods...'

LIBS="-L/usr/lib -llapack -L/usr/lib -lblas -L$HOME/libs/cgns/lib -lcgns -L$HOME/libs/hdf5/lib -lhdf5"
INCLUDES="-I$HOME/libs/cgns/include/ -I$HOME/libs/hdf5/include/ -I/usr/local/include"

LIBS=$LIBS" -L/home/jmuela/libs/gempa/lib -lgempa"
INCLUDES=$INCLUDES" -I/home/jmuela/libs/gempa/include"

MOD1="mod_mpi"
MOD2="mod_utils"
MOD3="mod_mpi_mesh"
MOD4="mod_comms"
MOD5="mod_cgns_mesh"

$COMPILER_F $INCLUDES -c $MODS_PATH$MOD1'.f90' -cpp -o $MODS_PATH$MOD1'.o' $LIBS #-J $MODS_PATH
$COMPILER_F $INCLUDES -c $MODS_PATH$MOD2'.f90' -cpp -o $MODS_PATH$MOD2'.o' $LIBS #-J $MODS_PATH
$COMPILER_F $INCLUDES -c $MODS_PATH$MOD3'.f90' -cpp -o $MODS_PATH$MOD3'.o' $LIBS #-J $MODS_PATH
$COMPILER_F $INCLUDES -c $MODS_PATH$MOD4'.f90' -cpp -o $MODS_PATH$MOD4'.o' $LIBS #-J $MODS_PATH
$COMPILER_F $INCLUDES -c $MODS_PATH$MOD5'.f90' -cpp -o $MODS_PATH$MOD5'.o' $LIBS #-J $MODS_PATH

MODS="$MODS_PATH$MOD1.o $MODS_PATH$MOD2.o $MODS_PATH$MOD3.o $MODS_PATH$MOD4.o $MODS_PATH$MOD5.o"
MODS="$MODS $GEMPA_INT.o"
#----------------------------------------------------------------------------------


echo 'Compiling '$PROGRAM
echo ' with mods: '$MODS

$COMPILER_F $MODS $INCLUDES $PROGRAM'.f90' $FLAGS_F -o $PROGRAM $LIBS

echo 'Compiling '$PROGRAM2

$COMPILER_F $MODS $INCLUDES $PROGRAM2'.f90' $FLAGS_F -o $PROGRAM2 $LIBS

         
echo 'Compilation done!'
