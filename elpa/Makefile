include ./Makefile.in

.PHONY: all clean

all: elpa_generalized

elpa_generalized: elpa_generalized.f90
	$(FC) $(FFLAGS) $< -o $@ $(LIBS)

clean:
	rm -rf *.o elpa_generalized
