include ../Makefile.in

LDFLAGS = $(FFLAGS)

TARGET = eigen

.PHONY: all clean dep

all: $(TARGET)

.SUFFIXES: .f .f90

OBJS_MAIN = main.o

OBJS_SOLVER = solver/solver_main.o solver/solver_lapack.o solver/solver_scalapack_all.o solver/solver_scalapack_select.o

OBJS_EIGENEXA = solver/solver_eigenexa.o
ifdef Use_EigenExa
	SRCS_EIGENEXA = solver/solver_eigenexa.f90
else
	SRCS_EIGENEXA = solver/solver_eigenexa_dummy.f90
endif

OBJS_UTIL_F77 = util/eigentest_pdlaprnt.o util/mmio.o

OBJS_UTIL = util/time.o util/descriptor_parameters.o util/command_argument.o util/matrix_io.o util/distribute_matrix.o util/processes.o util/eigenpairs_types.o

OBJS = $(OBJS_MAIN) $(OBJS_SOLVER) $(OBJS_EIGENEXA) $(OBJS_UTIL) $(OBJS_UTIL_F77)

$(OBJS_UTIL_F77): %.o: %.f
	$(FC) -c $(FFLAGS) $(LIBS) $< -o $@

$(OBJS_MAIN) $(OBJS_SOLVER) $(OBJS_UTIL): %.o: %.f90
	$(FC) -c $(FFLAGS) $(LIBS) $< -o $@

$(OBJS_EIGENEXA): $(SRCS_EIGENEXA)
	$(FC) -c $(FFLAGS) $(LIBS) $< -o $@

$(TARGET): $(OBJS)
	$(FC) $(LDFLAGS) $^ $(LIBS) -o $@
	@mkdir -p ../bin
	cp $(TARGET) ../bin/$(TARGET)

.PHONY: clean dep
clean:
	@rm -f *.exe *.o *.mod *~ ../bin/* */*.o */*.mod */*~ $(TARGET)

dep:
	makedepf90 -nosrc *.f90 util/*.f90 solver/*.f90 > Makefile.dep

include Makefile.dep
