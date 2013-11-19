include ../Makefile.in

LDFLAGS = $(FFLAGS)

TARGET = eigen

.PHONY: all clean dep

all: $(TARGET)

OBJS_MAIN = main.o

OBJS_SOLVER = solver/solver_main.o solver/solver_lapack.o solver/solver_scalapack_all.o solver/solver_scalapack_select.o

OBJS_EIGENEXA = solver/solver_eigenexa.o
ifdef Use_EigenExa
	SRCS_EIGENEXA = solver/solver_eigenexa.f90
	FFLAGS_EIGENEXA = $(FFLAGS)
else
	SRCS_EIGENEXA = solver/solver_eigenexa_dummy.f90
	FFLAGS_EIGENEXA = $(FFLAGS) -Wno-unused-dummy-argument
endif

OBJS_UTIL_F77 = util/eigentest_pdlaprnt.o util/mmio.o

OBJS_UTIL = util/time.o util/descriptor_parameters.o util/command_argument.o util/matrix_io.o util/distribute_matrix.o util/processes.o util/eigenpairs_types.o

OBJS = $(OBJS_MAIN) $(OBJS_SOLVER) $(OBJS_EIGENEXA) $(OBJS_UTIL) $(OBJS_UTIL_F77)

$(OBJS_UTIL_F77): %.o: %.f
	$(FC) -c $(FFLAGS) $(LIBS) $< -o $@

$(OBJS_MAIN) $(OBJS_SOLVER) $(OBJS_UTIL): %.o: %.f90
	$(FC) -c $(FFLAGS) $(LIBS) $< -o $@

$(OBJS_EIGENEXA): $(SRCS_EIGENEXA)
	$(FC) -c $(FFLAGS_EIGENEXA) $(LIBS) $< -o $@

$(TARGET): $(OBJS)
	$(FC) $(LDFLAGS) $^ $(LIBS) -o $@
	@mkdir -p ../bin
	cp $(TARGET) ../bin/$(TARGET)

clean:
	@rm -f *.exe *.o *.mod *~ ../bin/* */*.o */*.mod */*~ $(TARGET)

dep:
	find . -name \*.f90 -and ! -name \*_dummy.f90 | xargs makedepf90 -nosrc > Makefile.dep

include Makefile.dep
