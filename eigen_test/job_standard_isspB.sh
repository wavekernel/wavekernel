#!/bin/bash
#BSUB -J eigen_test
#BSUB -q F4
#BSUB -n 4

# Usage: bsub -W 1 < job.sh
export OMP_NUM_THREADS=8

matrix_dir=/work/k0099/k009914/matrix
size=400
matrix=$matrix_dir/ELSES_MATRIX_VCNT${size}std_A.mtx

hybrid -np 2 "./bin/eigen -s scalapack_all $matrix > out/o_eigen"
