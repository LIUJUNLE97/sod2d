#!/bin/bash
if [[ $OMPI_COMM_WORLD_RANK == 0 ]]; then
   ncu --set full -k regex:"elem_diffu_full_diffusion_ijk_539|elem_convec_full_convec_ijk_934" --launch-skip 4 --launch-count 2 -o kernels_${OMPI_COMM_WORLD_RANK} -f --target-processes all "$@"
else
   "$@"
fi
