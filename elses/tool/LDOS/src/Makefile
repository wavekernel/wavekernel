#================================================================
# ELSES version 0.04
# Copyright (C) ELSES. 2007-2014 all rights reserved
#================================================================
#
# Makefile for Tools for ELSES
#

TOPDIR= ../../..
include $(TOPDIR)/Makefile.inc
TARGETS= Green_R Green_C Green_SSW Green_diag
BINDIR=$(TOPDIR)/bin

F90 ?= $(FC)

.f.o:
	$(F90) $(FFLAGS) -c $<

.PHONY: all clean

all: $(TARGETS)

clean:
	@rm -f *.o *.mod $(TARGETS)


Green_R    : Green_R.o subroutines.o
	$(F90) $(LDFLAGS) -o $@ Green_R.o subroutines.o $(LIBS)
	@cp $@ $(BINDIR)/$@

Green_C    : Green_C.o subroutines.o
	$(F90) $(LDFLAGS) -o $@ Green_C.o subroutines.o $(LIBS)
	@cp $@ $(BINDIR)/$@

Green_SSW  : Green_SSW.o subroutines.o
	$(F90) $(LDFLAGS) -o $@ Green_SSW.o subroutines.o $(LIBS)
	@cp $@ $(BINDIR)/$@

Green_diag : Green_diag.o subroutines.o diag.o
	$(F90) $(LDFLAGS) -o $@ Green_diag.o subroutines.o diag.o $(LIBS)
	@cp $@ $(BINDIR)/$@

Green_Lanc : Green_Lanc.o subroutines.o
	$(F90) $(LDFLAGS) -o $@ Green_Lanc.o subroutines.o $(LIBS)
	@cp $@ $(BINDIR)/$@

Green_R.o : subroutines.o
Green_C.o : subroutines.o
Green_SSW.o : subroutines.o
Green_diag.o : subroutines.o diag.o
Green_Lanc.o : subroutines.o
