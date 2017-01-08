#================================================================
# ELSES version 0.05
# Copyright (C) ELSES. 2007-2014 all rights reserved
#================================================================
#
# Makefile for Tools for ELSES
#

TOPDIR=..
include $(TOPDIR)/Makefile.inc

TARGET_PHONY =\
 tool-LDOS tool-Perl tool-band tool-eigenstates\
 tool-DOS tool-MSD tool-src tool-xml-check 

TARGET =

.PHONY: all clean $(TARGET_PHONY)

all: $(TARGET) $(TARGET_PHONY)

clean:
	@rm -f *.exe src/*.o src/*.mod src/*~ *.mod $(TARGET)
	@$(MAKE) -C MSD/src clean
	@$(MAKE) -C band/src clean
	@$(MAKE) -C plot_eigenstates/src clean
	@$(MAKE) -C LDOS/src clean
	@$(MAKE) -C DOS -f Makefile clean
	@$(MAKE) -C src clean
	@$(MAKE) -C xml-check/src clean


tool-src:
	make -C src all

tool-LDOS:
	@$(MAKE) -C LDOS/src all

tool-DOS:
	@$(MAKE) -C DOS all

tool-band:
	@$(MAKE) -C band/src  all

tool-Perl:
	@cp Perl/*.pl ../bin

tool-eigenstates:
	@$(MAKE) -C plot_eigenstates/src all

tool-MSD:
	@$(MAKE) -C MSD/src all

tool-xml-check: tool-src
	@$(MAKE) -C xml-check/src -f Makefile all

#tool-Perl:
#	( pwd=`pwd`; for i in $${pwd}/tool/Perl/*.pl ; do fn=`basename $$i`;\
#	if [ ! -L bin/$$fn ] ; then ln -s $$i bin/$$fn ; fi ; done )
