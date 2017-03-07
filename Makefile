-include ./Makefile.inc

.PHONY: all main lib clean dep

MAIN = bin/wavepacket
LIB = libwavekernel.a
OBJS_MAIN = src/main.o
OBJS_LIB = \
	src/main_aux.o \
	src/initialization.o \
	src/matrix_generation.o \
	src/conversion.o \
	src/state.o \
	src/charge.o \
	src/atom.o \
	src/output.o \
	src/time_evolution.o \
	src/linear_algebra.o \
	src/matrix_io.o \
	src/stats.o \
	src/util.o \
	src/setting.o \
	src/event_logger.o \
	src/fson.o \
	src/distribute_matrix.o \
	src/processes.o \
	src/descriptor_parameters.o \
	src/global_variables.o
OBJS_LIB_F77 = src/mmio.o src/eigentest_pdlaprnt.o

all: main lib

main: $(MAIN)

lib: $(LIB)

$(MAIN): $(OBJS_MAIN) lib
	@mkdir -p bin
	$(FC) $(LDFLAGS) -o $@ $< $(LIB) $(LIBS)

$(OBJS_MAIN): %.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJS_LIB): %.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJS_LIB_F77): %.o: %.f
	$(FC) -c $(FFLAGS_NO_WARN) $< -o $@

$(LIB): $(OBJS_LIB) $(OBJS_LIB_F77)
	$(AR) r $@ $^

clean:
	@rm -f $(MAIN) $(LIB) src/*.o *.mod

dep:
	find . -name \*.f90 | xargs makedepf90 -nosrc > src/Makefile.dep

include src/Makefile.dep
