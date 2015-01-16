#================================================================
# ELSES version 0.04
# Copyright (C) ELSES. 2007-2014 all rights reserved
#================================================================
TOPDIR=../../..
include $(TOPDIR)/Makefile.inc

TARGET_PHONY = 
TARGET = xml-check

OBJS = elses-xml-check.o

.PHONY: all clean $(TARGET_PHONY)

FFLAGS += $(FFLAGS-NRL) $(FFLAGS-XML)
LIBS += $(OBJS-XML) $(LIBS-XML)

F90 ?= $(FC)
.SUFFIXES: .f90 .f

.f90.o:
	$(F90) -c $(FFLAGS) $< -o $@

.f.o:
	$(F77) -c $(FFLAGS) $< -o $@

all: $(TARGET)

clean:
	@rm -f *.o *.mod $(TARGET)

xml-check : $(OBJS)
	$(F90) $(LDFLAGS) $(OBJS) $(LIBS) -o $@
	@cp $@ $(BINDIR)
