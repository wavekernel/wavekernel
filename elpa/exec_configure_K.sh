#!/bin/sh

. /work/system/Env_base

./configure --prefix=$HOME MPIFC=mpifrt SCALAPACK_LDFLAGS="-SCALAPACK -SSL2BLAMP" SCALAPACK_FCFLAGS="-SCALAPACK -SSL2BLAMP" FCFLAGS="-Kfast" --with-only-real-generic-kernel --with-only-complex-generic-kernel --build=sparcv8
