FC = mpif90
FFLAGS = -O3
#FFLAGS = -O3 -mcmodel=medium -i-dynamic -openmp
#FFLAGS = -g -assume nobuffered_io -check all -traceback -ftrapuv
LDFLAGS = $(FFLAGS)
LIBS = -L /usr/lib/openmpi/lib -lscalapack -llapack -lrefblas -lgomp -lmpi
#LIBS = -L /opt/intel/Compiler/11.1/069/mkl/lib/em64t -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#LIBS  = /opt/acml4.4.0/ifort64_mp/lib/libacml_mp.a /opt/acml4.4.0/ifort64_mp/lib/libacml_mv.a  -lguide -lpthread
######################################################
# Setting for MacOSX with gfortran
#FC = gfortran
#FFLAGS = -Wall -pedantic -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
#LDFLAGS = -fopenmp
#LIBS = -framework vecLib
######################################################
#
TARGET = eigen

.PHONY: all
all: $(TARGET)

F90 ?= $(FC)

.SUFFIXES: .f90

.f90.o:
	$(F90) -c $(FFLAGS) $< -o $@

OBJS-TEST = eigen-test-main.o

include lib_misc/Makefile.in
include lib_eigen/Makefile.in

eigen : $(OBJS-TEST) $(OBJS-LIB-EIGEN) $(OBJ-LIB-MISC)
	$(F90) $(LDFLAGS) $(OBJS-TEST) $(OBJS-LIB-EIGEN) $(OBJ-LIB-MISC) $(LIBS) -o $@
	cp $@ ../bin/

.PHONY: clean dep
clean:
	@rm -f *.exe *.o *.mod *~ ../bin/* lib_eigen/*.o lib_eigen/*.mod lib_eigen/*~ lib_misc/*.o lib_misc/*.mod $(TARGET)

dep:
	makedepf90 -nosrc *.f90 lib_eigen/*.f90 lib_misc/*.f90 > Makefile.dep

include Makefile.dep
