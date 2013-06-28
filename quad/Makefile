include ./Makefile.in

TARGET = pi-quad-prec

.PHONY: all
all: $(TARGET)

F90 ?= $(FC)

.SUFFIXES: .f90

.f90.o:
	$(F90) -c $(FFLAGS) $(LIBS) $< -o $@

OBJ = pi-quad-prec.o

$(TARGET) : $(OBJ)
	$(F90) $(LDFLAGS) $(OBJ) $(LIBS) -o $@

.PHONY: clean
clean:
	@rm -f $(TARGET) *~ *.o
