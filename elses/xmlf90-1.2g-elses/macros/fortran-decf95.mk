#
# Fortran macros for DEC's f95
#
.SUFFIXES: .SUFFIXES .f .f90 .f95 .F .F90 .F95 .fpp
#
LIBRARY=libflib.a
#
LIB_STD=$(FLIB_ROOT)/lib/
MOD_STD=$(FLIB_ROOT)/modules/
INC_STD=$(FLIB_ROOT)/include/
BIN_STD=$(FLIB_ROOT)/bin/
#
FC=f95
CFLAGS= 
FFLAGS= -fast -tune host   $(CFLAGS)
FFLAGS_DEBUG= -g 
LDFLAGS=      $(CFLAGS)
#
INC_PREFIX=-I
MOD_PREFIX=-I
LIB_PREFIX=-L
#
MOD_EXT=mod
MOD_SEARCH_STD= $(MOD_PREFIX)$(MOD_STD) 
MOD_SEARCH= $(MOD_SEARCH_STD) $(MOD_SEARCH_OTHER)
#INC_SEARCH= $(INC_PREFIX)$(INC_STD)
#
#
AR=ar
RANLIB=ranlib
#
CPP=/bin/cpp -P
COCO=$(BIN_STD)coco -I $(INC_STD)
DEFS=
#
.F.o:
	$(CPP) $(DEFS) $< > $*.f
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $*.f
	@rm -f $*.f
.f.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
.F90.o:
	$(CPP) $(DEFS) $< > $*.f90
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $*.f90
	@rm -f $*.f90
.f90.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
.F95.o:
	$(CPP) $(DEFS) $< > $*.f90
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $*.f90
	@rm -f $*.f90
.f95.o:
	cp $< $*.f90
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $*.f90
	@rm -f $*.f90
#
.fpp.f90:
	$(COCO) $*
#









