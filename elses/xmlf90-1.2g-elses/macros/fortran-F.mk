#
# Fortran macros for F on Linux.
# The F compiler is based on NAG's f95.
#
# This macro file defines __F__ and __NAG__ for conditional compilation
#
#
.SUFFIXES: .SUFFIXES .f .f90 .f95 .F .F90 .F95 .fpp
#
LIB_STD=$(FLIB_ROOT)/lib/
MOD_STD=$(FLIB_ROOT)/modules/
INC_STD=$(FLIB_ROOT)/include/
BIN_STD=$(FLIB_ROOT)/bin/
#
FC=F
F77=g77
CFLAGS= 
FFLAGS= -O  $(CFLAGS)
F77_FLAGS= -O -f-no-second-underscore
FFLAGS_DEBUG= -g  
FFLAGS_CHECK= -g  --chk s,u,e,a
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
AR=ar
RANLIB=ranlib
#
CPP=/bin/cpp -P
DEFS= -D__NAG__ -D__F__
#
# Experimental : the following deactivates an implicit rule
# which breaks havoc with the operation of this makefile
# It works at least with GNU make
%.o : %.mod
#
.F.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)  $(DEFS) $<
.f.o:
	$(F77) -c $(INC_SEARCH) $(F77_FFLAGS)   $<
.F90.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(DEFS) $<
.f90.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
.F95.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(DEFS) $<
.f95.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
#









