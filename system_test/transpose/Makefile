include ../../eigen_test/Makefile.in

.PHONY: all clean

BIN = transpose

all: $(BIN)

$(BIN): main.f90 processes.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

processes.o: ../../eigen_test/src/util/processes.f90
	$(FC) $(FFLAGS) -c -o $@ $^ $(LIBS)

clean:
	rm -f $(BIN) *.o
