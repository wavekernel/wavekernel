include ../Makefile.in

LDFLAGS = $(FFLAGS)

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

$(TARGET): $(OBJS-TEST) $(OBJS-LIB-EIGEN) $(OBJ-LIB-MISC)
	$(F90) $(LDFLAGS) $^ $(LIBS) -o $@
	@mkdir -p ../bin
	cp $(TARGET) ../bin/$(TARGET)

.PHONY: clean dep
clean:
	@rm -f *.exe *.o *.mod *~ ../bin/* */*.o */*.mod */*~ $(TARGET)

dep:
	makedepf90 -nosrc *.f90 util/*.f90 solver/*.f90 > Makefile.dep

include Makefile.dep
