#!/bin/bash
module load openmpi/1.10.2-gnu-pmi2
srun ./heat_2d_loc_blk_MPI
# execute a 10000 x 10000 point 2d-heat transfer problem
# for 5 iternations and suppress its output
