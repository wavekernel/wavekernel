-include ../../../Makefile.inc
include ../Makefile.in

.PHONY: all bin lib clean dep

BIN = ../bin/wavepacket
LIB = libwavepacket.a
OBJS_LIB = main_aux.o initialization.o matrix_generation.o conversion.o state.o charge.o atom.o output.o time_evolution.o linear_algebra.o matrix_io.o stats.o util.o setting.o event_logger.o fson.o distribute_matrix.o processes.o descriptor_parameters.o global_variables.o
OBJS_LIB_F77 = mmio.o eigentest_pdlaprnt.o

all: bin

bin: $(BIN)

lib: $(LIB)

$(BIN): main.f90 lib
	$(FC) $(FFLAGS) -o $@ $< $(LIB) $(LIBS)

$(OBJS_LIB_F77): %.o: %.f
	$(FC) -c $(FFLAGS_NO_WARN) $< -o $@

$(OBJS_LIB): %.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

$(LIB): $(OBJS_LIB) $(OBJS_LIB_F77)
	mkdir -p `dirname $(BIN)`
	$(AR) r $@ $^

clean:
	rm -f $(BIN) $(LIB) *.o *.mod

dep:
	find . -name \*.f90 | xargs makedepf90 -nosrc > Makefile.dep

include Makefile.dep
