#!/bin/bash
#BSUB -J eigen_test
#BSUB -q F4
#BSUB -n 4
#BSUB -o stdout
#BSUB -e stderr

# Usage: bsub -W 5 < job_standard_isspB.sh
export OMP_NUM_THREADS=8

matrix_dir= # Set properly
size=400
matrix=$matrix_dir/ELSES_MATRIX_VCNT${size}std_A.mtx

hybrid -np 2 "./bin/eigen -s scalapack_all -c 10 -p 1,10 $matrix"
