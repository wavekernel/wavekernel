#================================================================
# ELSES version 0.04
# Copyright (C) ELSES. 2007-2014 all rights reserved
#================================================================
#
# Makefile for ELSES
#

include ../../../Makefile.inc

TARGET = elses-generate-cubefile

.PHONY: all


all: $(TARGET)

F90 ?= $(FC)
F77 ?= $(FC)

.SUFFIXES: .f90 .f

.f90.o:
	$(F90) -c $(FFLAGS) $< -o $@

.f.o:
	$(F77) -c $(FFLAGS) $< -o $@

OBJS = elses-generate-cubefile.o elses-lib-math-cubic-harmonics.o elses-lib-sto-radial-func.o

elses-generate-cubefile: $(OBJS)
	$(F90) $(LDFLAGS) $(OBJS) $(LIBS) -o $@
	cp $@ ../../../bin/

.PHONY: clean dep etags
clean:
	@rm -f *.o *.mod *~ $(TARGET)
	@rm -f ../bin/$(TARGET)

dep:
	makedepf90 -nosrc *.f90 > Makefile.dep

etags:
	etags *.{f90}

include Makefile.dep
