include ./Makefile.in

.PHONY: all clean

all: test_real2

test_real2: test_real2.f90
	$(FC) $(FFLAGS) $< -o $@ -I. -L. -lelpa $(LIBS)

clean:
	rm -rf *.o test_real2
