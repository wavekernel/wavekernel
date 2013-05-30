include ./Makefile.in

LDFLAGS = $(FFLAGS)

TARGET = test_blacs_overhead.out

OBJS = test_blacs_overhead.o

.PHONY: all
all: $(TARGET)

.SUFFIXES: .f90

.f90.o:
	$(FC) -c $(FFLAGS) $< -o $@

$(TARGET) : $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

.PHONY: clean dep
clean:
	@rm -f *.o *~ $(TARGET)
