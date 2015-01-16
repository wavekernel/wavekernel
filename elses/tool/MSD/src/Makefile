#================================================================
# ELSES version 0.04
# Copyright (C) ELSES. 2007-2014 all rights reserved
#================================================================

include ../../../Makefile.inc

TARGET = msd

F90 ?= $(FC)

OBJS = calc_msd.o

.PHONY: all clean

.SUFFIXES: .f90

.f90.o:
	$(F90) $(FFLAGS) -c $<

all: $(TARGET)

clean:
	@rm -f *.o *~ $(TARGET)
	@rm -f ../../../bin/$(TARGET)

msd: $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
	@cp $@ ../../../bin/
