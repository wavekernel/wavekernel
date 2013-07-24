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

matrix_dir=/work/k0099/k009914/matrix
size=400
matrix=$matrix_dir/ELSES_MATRIX_VCNT${size}std_A.mtx

#mpiexec ./bin/eigen -s scalapack_select -n 200 -c 5 $matrix
mpiexec ./bin/eigen -s scalapack_all $matrix
