#!/bin/sh
ELSES_SRC=./src

if [ "$1" = "true" ]; then
    echo "switch to true"
    cp -p ${ELSES_SRC}/code_for_mpi/elses-lib-mpi-wrapper.f90 ${ELSES_SRC}/elses-lib-mpi-wrapper-compile.f90
    cp -p ${ELSES_SRC}/code_optional/elses-qm-wavepacket-main-true.f90 ${ELSES_SRC}/elses-qm-wavepacket-main.f90
    cp -p ${ELSES_SRC}/code_optional/elses-la-eigen-solver-main-true.f90 ${ELSES_SRC}/elses-la-eigen-solver-main.f90
elif [ "$1" = "mpi" ]; then
    echo "switch to mpi"
    cp -p ${ELSES_SRC}/code_for_mpi/elses-lib-mpi-wrapper.f90 ${ELSES_SRC}/elses-lib-mpi-wrapper-compile.f90
    cp -p ${ELSES_SRC}/code_optional/elses-qm-wavepacket-main-dummy.f90 ${ELSES_SRC}/elses-qm-wavepacket-main.f90
    cp -p ${ELSES_SRC}/code_optional/elses-la-eigen-solver-main-dummy.f90 ${ELSES_SRC}/elses-la-eigen-solver-main.f90
elif [ "$1" = "dummy" ]; then
    echo "switch to dummy"
    cp -p ${ELSES_SRC}/code_for_mpi/elses-lib-mpi-wrapper-dummy.f90 ${ELSES_SRC}/elses-lib-mpi-wrapper-compile.f90
    cp -p ${ELSES_SRC}/code_optional/elses-qm-wavepacket-main-dummy.f90 ${ELSES_SRC}/elses-qm-wavepacket-main.f90
    cp -p ${ELSES_SRC}/code_optional/elses-la-eigen-solver-main-dummy.f90 ${ELSES_SRC}/elses-la-eigen-solver-main.f90
else
    echo "command: true / mpi / dummy"
fi
