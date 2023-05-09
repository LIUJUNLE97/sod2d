#!/bin/bash

# Use $PMI_RANK for MPICH and $SLURM_PROCID with srun.
if [ $OMPI_COMM_WORLD_RANK -eq 0 ]; then
  nsys profile -e NSYS_MPI_STORE_TEAMS_PER_RANK=1 -t cuda,openacc,nvtx,mpi --cuda-memory-usage=true -o tgv_40x40 "$@"
else
  "$@"
fi
