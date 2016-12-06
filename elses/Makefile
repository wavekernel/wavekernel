#================================================================
# ELSES version 0.06
# Copyright (C) ELSES. 2007-2015 all rights reserved
#================================================================
#
# Makefile for ELSES
#

VERSION=0.06.03

include Makefile.inc

PackageContents=xml elses tool ext-lib

.PHONY: all clean clean_elses clean_tool clean_extlib release revision $(PackageContents)

all: $(PackageContents)

ext-lib:
	@$(MAKE) -C ext-lib all

xml:
	@$(MAKE) -C $(XML-DIR) -f Makefile-elses all

elses: ext-lib xml
	@$(MAKE) -C src all

# Type "make relase" to make a release package.
#
# The directory named "sample" is excluded from the package.
# To append "sample" directory to the package,
# insert the string "sample,"
# just after the string "$(CUR_DIR)/{bin,src,"
# at the last line of following block.
CUR_DIR=elses-v$(VERSION)
TODAY=`date +%Y%m%d`
PACKAGE_NAME=elses-v$(VERSION)-$(TODAY)-rel.tgz
release: clean
	@touch 00_elses_build_time.txt
	@rm -f bin/*
	@diff -q Makefile.inc Makefile.inc_release; \
	if [ $$? -ne 0 ] ; then \
      echo "Makefile.inc is modified.";\
      echo "If it is what you expected, renew Makefile.inc_release, too.";\
      echo "If not so, cp Makefile.inc_release Makefile.inc, by your own hand."; \
      exit 1; fi
	@$(MAKE) -C src ctree
	@(grep 'version_info='  src/elses-00-version-info.f90 | grep -q $(VERSION)) ;\
	  if [ $$?  -ne 0 ]; \
      then echo "Inconsistent VERSION. Check VERSION in Makefile."; exit 1 ; fi
	@if [ `basename \`pwd\` ` != "$(CUR_DIR)" ] ; \
       then echo "Inconsistent VERSION. Check the name of the current directory." ;\
       exit 1 ; fi
	@cd ../ ;\
      tar czf $(PACKAGE_NAME) \
      --exclude '*.o' --exclude '*.mod' --exclude '*,v' --exclude '*~' \
      --exclude '*old' --exclude 'doc/msc/src' --exclude 'tool/TB-PAM/doc/doc_work' \
      $(CUR_DIR)/{bin,ext-lib,src,tool,xml*,Makefile*,sample,doc,*.patch,*.ReadMe,*.txt}
	@echo "Release package $(PACKAGE_NAME) is created in the parent directory."

revision: clean
	@rm -f bin/* tool/LDOS/src/Green_{C,R,SSW,diag}
	@cp Makefile.inc Makefile.inc.old
	@cp Makefile.inc_release Makefile.inc
	@(grep 'version_info='  src/elses-00-version-info.f90 | grep -q $(VERSION)) ;\
	  if [ $$?  -ne 0 ]; \
      then echo "Inconsistent VERSION. Check VERSION in Makefile."; exit 1 ; fi
	@if [ `basename \`pwd\` ` != "$(CUR_DIR)" ] ; \
       then echo "Inconsistent VERSION. Check the name of the current directory." ;\
       exit 1 ; fi
	@cd ../ ;\
      tar czf $(PACKAGE_NAME) \
      --exclude '*.o' --exclude '*.mod' --exclude '*,v' --exclude '*~' \
      --exclude '*old' --exclude 'doc/msc/src' \
      $(CUR_DIR)/{bin,ext-lib,src,tool,xml*,Makefile*,doc,*.patch,*.ReadMe,*.txt}
	@echo "Release package $(PACKAGE_NAME) is created in the parent directory."

rc: release
	@echo "Release package $(PACKAGE_NAME) is renamed to $(PACKAGE_NAME:%-rel.tgz=%-rc.tgz)"
	@cd ..; mv $(PACKAGE_NAME) $(PACKAGE_NAME:%-rel.tgz=%-rc.tgz)

tool: xml
	@$(MAKE) -C tool all

clean: clean_elses clean_tool clean_extlib

clean_elses:
	@$(MAKE) -C $(XML-DIR) -f Makefile-elses clean
	@$(MAKE) -C src clean

clean_tool:
	@$(MAKE) -C tool clean

clean_extlib:
	@$(MAKE) -C ext-lib clean
