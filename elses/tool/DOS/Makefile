#================================================================
# ELSES version 0.05
# Copyright (C) ELSES. 2007-2014 all rights reserved
#================================================================
#
# Makefile for ELSES
#

TOPDIR = ../..
include $(TOPDIR)/Makefile.inc

TARGET = elses-dos

.PHONY: all

all: $(TARGET)

FFLAGS +=  $(FFLAGS-Nrl) $(FFLAGS-XML)
LIBS += $(LIBS-XML)

#for debugging
#FFLAGS += -traceback

F90 ?= $(FC)
F77 ?= $(FC)

.SUFFIXES: .f90 .f

.f90.o:
	$(F90) -c $(FFLAGS) $< -o $@

.f.o:
	$(F77) -c $(FFLAGS) $< -o $@

OBJS = $(shell ls $(TOPDIR)/src/elses-mod-file-io.o $(TOPDIR)/src/elses-lib-phys-const.o $(TOPDIR)/src/elses-lib-math-func.o $(TOPDIR)/src/elses-io-ctrl-dos-plot.o |grep -v 00-main|grep -v selftest)
OBJS += elses-dos.o

elses-dos: $(OBJS)
	$(F90) $(LDFLAGS) $(OBJS) $(LIBS) -o $@
	cp $@ $(BINDIR)

.PHONY: clean dep etags
clean:
	@rm -f *.o *.mod *~ $(TARGET)
	@rm -f ../bin/$(TARGET)

dep:
	makedepf90 -nosrc *.f90 ../../src/*.{f,f90} > Makefile.dep

etags:
	etags *.{f90} ../../src/*.{f,f90}

include Makefile.dep
