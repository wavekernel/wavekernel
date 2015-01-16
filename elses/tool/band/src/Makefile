#================================================================
# ELSES version 0.04
# Copyright (C) ELSES. 2007-2014 all rights reserved
#================================================================
#
# Makefile for Tools for ELSES
#

TOPDIR= ../../..
include $(TOPDIR)/Makefile.inc
TARGETS= band dos
BINDIR=$(TOPDIR)/bin

F90 ?= $(FC)
OBJS = plot_data.o band.o diag.o

.PHONY: all clean
.SUFFIXES: .f90

.f90.o:
	$(F90) $(FFLAGS) -c $<

all: $(TARGETS)

clean:
	@rm -f *.o *.mod $(TARGETS)

band : $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
	@cp $@ $(BINDIR)/$@

dos : dos.o
	$(F90) $(LDFLAGS) -o $@ dos.o $(LIBS)
	@cp $@ $(BINDIR)/$@

band.o: diag.o



