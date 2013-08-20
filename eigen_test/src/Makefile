include ../Makefile.in

#FFLAGS = -O3 -mcmodel=medium -i-dynamic -openmp
#FFLAGS = -g -assume nobuffered_io -check all -traceback -ftrapuv
LDFLAGS = $(FFLAGS)
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

.SUFFIXES: .f .f90

.f.o:
	$(F90) -c $(FFLAGS) $(LIBS) $< -o $@

.f90.o:
	$(F90) -c $(FFLAGS) $(LIBS) $< -o $@

OBJS-TEST = main.o

include util/Makefile_util.in
include solver/Makefile_solver.in

eigen : $(OBJS-TEST) $(OBJS-LIB-EIGEN) $(OBJ-LIB-MISC)
	$(F90) $(LDFLAGS) $(OBJS-TEST) $(OBJS-LIB-EIGEN) $(OBJ-LIB-MISC) $(LIBS) -o $@
	mkdir -p ../bin
	mv $@ ../bin/

.PHONY: clean dep
clean:
	@rm -f *.exe *.o *.mod *~ ../bin/* */*.o */*.mod */*~ $(TARGET)

dep:
	makedepf90 -nosrc *.f90 util/*.f90 solver/*.f90 > Makefile.dep

include Makefile.dep
