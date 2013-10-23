#!/bin/bash
#BSUB -J eigen_test
#BSUB -q F4
#BSUB -n 4
#BSUB -o stdout
#BSUB -e stderr

# Usage: bsub -W 5 < job_general_isspB.sh
export OMP_NUM_THREADS=8

matrix_dir= # Set properly
size=900
matrixA=$matrix_dir/ELSES_MATRIX_VCNT${size}_A.mtx
matrixB=$matrix_dir/ELSES_MATRIX_VCNT${size}_B.mtx

hybrid -np 2 "./bin/eigen -s general_scalapack_all -c 10 -p 1,10 $matrixA $matrixB"
