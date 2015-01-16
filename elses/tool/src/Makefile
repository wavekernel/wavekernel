#================================================================
# ELSES version 0.05
# Copyright (C) ELSES. 2007-2014 all rights reserved
#================================================================
#
# Makefile for Tools for ELSES
#

TOPDIR=../..
include $(TOPDIR)/Makefile.inc

TARGETS_PHONY =

TARGETS = elses-xml-generate elses-xml-generate-vasp \
 elses-xml-generate-nrl-param elses-rdf elses-msd

.PHONY: all clean $(TARGETS_PHONY)

all: $(TARGETS) $(TARGETS_PHONY)

FFLAGS += $(FFLAGS-XML)

LIBS = $(OBJS-XML) $(LIBS-XML)

F90 ?= $(FC)
F77 ?= $(FC)

.SUFFIXES: .f90 .f

.f90.o:
	$(F90) -c $(FFLAGS) $< -o $@

.f.o:
	$(F77) -c $(FFLAGS) $< -o $@

all: $(TARGETS)
	@cp $(TARGETS) $(BINDIR)

elses-xml-generate: elses-xml-generate.o
	$(F90) $(LDFLAGS) elses-xml-generate.o $(LIBS) -o $@

elses-xml-generate-vasp: elses-xml-generate-vasp.o
	$(F90) $(LDFLAGS) elses-xml-generate-vasp.o $(LIBS) -o $@

elses-xml-generate-nrl-param: elses-xml-generate-nrl-param.o
	$(F90) $(LDFLAGS) elses-xml-generate-nrl-param.o $(LIBS) -o $@

elses-rdf: elses-rdf.o
	$(F90) $(LDFLAGS) elses-rdf.o $(LIBS) -o $@

elses-msd: elses-msd.o
	$(F90) $(LDFLAGS) elses-msd.o $(LIBS) -o $@

clean:
	@rm -f *.o $(TARGETS)
