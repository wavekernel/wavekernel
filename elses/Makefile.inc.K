#================================================================
# ELSES version 0.06
# Copyright (C) ELSES. 2007-2015 all rights reserved
#================================================================

TOPDIR ?= .

# Common settings
MOD_OPT=-I
#XML_SELFTEST = './elses-xml-selftest'
# For cross compiler systems,
XML_SELFTEST = echo './elses-xml-selftest'
%.o:%.mod

XML-DIR = xmlf90-1.2g-elses

BINDIR=$(TOPDIR)/bin
SRCDIR=$(TOPDIR)/src

FFLAGS-NRL =$(MOD_OPT)$(TOPDIR)/src/Nrl
FFLAGS-XML =\
 $(MOD_OPT)$(TOPDIR)/$(XML-DIR)/strings\
 $(MOD_OPT)$(TOPDIR)/$(XML-DIR)/sax\
 $(MOD_OPT)$(TOPDIR)/$(XML-DIR)/dom\
 $(MOD_OPT)$(SRCDIR)
LIBS-XML = -L$(TOPDIR)/$(XML-DIR)/macros/lib -lflib
OBJS-XML = $(SRCDIR)/elses-xml-config.o $(SRCDIR)/elses-config.o\
 $(SRCDIR)/elses-xml-misc.o $(SRCDIR)/elses-xml-sax-counter.o\
 $(SRCDIR)/elses-xml-sax-handler.o $(SRCDIR)/elses-xml-sax-parser.o

# Dummy setting for LIBS-EXT 
LIBS-EXT = 

## Fujitsu FX10
FC = mpifrtpx
FFLAGS = -Kfast,openmp
LDFLAGS = $(FFLAGS)
LIBS = -SSL2

