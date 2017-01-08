#!/bin/bash

#------ pjsub option --------#

#PJM -L "rscgrp=debug"
#PJM -L "node=2x2"
#PJM --mpi "proc=4"
#PJM -L "elapse=10:00"
#PJM -e out/e%J
#PJM -o out/o%J
#PJM -s
#PJM --spath out/i%J

#------ program execution ------#
export OMP_NUM_THREADS=16

mkdir -p out
matrix_dir= # Set properly
size=400
matrix=$matrix_dir/ELSES_MATRIX_VCNT${size}std_A.mtx

mpiexec ./bin/eigen -s scalapack_all -c 10 -p 1,10 $matrix
