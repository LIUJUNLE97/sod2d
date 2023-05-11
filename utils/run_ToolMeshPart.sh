#!/bin/bash

## Small script to rrun the mesh_part tool of Jordi.
## Based on the arguments, it will create a "part.dat"
## file with the data required for the tool to run:

## - gmsh_filePath ""
## - gmsh_fileName "cube_per40"
## - mesh_h5_filePath ""
## - mesh_h5_fileName "cube_per40"
## - num_partitions 1

## Usage: ./run_ToolMeshPart.sh -N [numCPUs]
## -N [numCPUs] : Number of CPUs to run the tool.

## Check that the -N flag is passed
if [ $1 != "-N" ]; then
    echo "Usage: ./run_ToolMeshPart.sh -N [numCPUs]"
    echo "  -N [numCPUs] : Number of CPUs to run the tool."
    exit 1
fi
## Check that the number of CPUs is passed
if [ -z $2 ]; then
    echo "Usage: ./run_ToolMeshPart.sh -N [numCPUs]"
    echo "  -N [numCPUs] : Number of CPUs to run the tool."
    exit 1
elif [ $2 -lt 1 ]; then
    echo "Error: Number of CPUs must be greater than 0."
    exit 1
fi

## If a part.dat file does not exist, create it
if [ ! -f part.dat ]; then
    touch part.dat
    ## Ask the user for the path of the input *.h5 mesh
    echo "Please enter the path of the *.h5 mesh file:"
    read gmsh_filePath

    ## Ask for the name of the input *.h5 mesh
    echo "Please enter the name of the *.h5 mesh file:"
    read gmsh_fileName

    ## Ask for the path of the output *.h5 mesh (partitioned)
    echo "Please enter the path of the output *.h5 mesh file:"
    read mesh_h5_filePath

    ## Ask for the name of the output *.h5 mesh (partitioned)
    echo "Please enter the name of the output *.h5 mesh file:"
    read mesh_h5_fileName

    ## Ask for the number of partitions
    echo "Please enter the number of desired mesh partitions:"
    read num_partitions

    ## Check number of partitions
    ## Entry is not a number
    if ! [[ $num_partitions =~ ^[0-9]+$ ]]; then
        echo "ERROR: Number of partitions must be an integer."
        exit 1
    ## Entry is smaller than 1
    elif [ $num_partitions -lt 1 ]; then
        echo "ERROR: Number of partitions must be greater than 0."
        exit 1
    ## Entry is smaller than number of requested CPUs
    elif [ $num_partitions -lt $2 ]; then
        echo "WARNING: Too many CPUs requested! Some resources will be idle."
    fi

    ## Create the part.dat file
    echo "gmsh_filePath \"$gmsh_filePath\"" > part.dat
    echo "gmsh_fileName \"$gmsh_fileName\"" >> part.dat
    echo "mesh_h5_filePath \"$mesh_h5_filePath\"" >> part.dat
    echo "mesh_h5_fileName \"$mesh_h5_fileName\"" >> part.dat
    echo "num_partitions $num_partitions" >> part.dat
fi

## Request the path of the tool_meshConversorPar executable
echo "Please enter the path of the tool_meshConversorPar executable:"
read exec_path

## Check that the tool_meshConversorPar is on the requested path
if [ ! -f $exec_path/tool_meshConversorPar ]; then
    echo "ERROR: tool_meshConversorPar not found on the requested path."
    exit 1
fi

## Run the tool
echo "Running the tool..."
mpirun -np $2 $exec_path/tool_meshConversorPar part.dat
