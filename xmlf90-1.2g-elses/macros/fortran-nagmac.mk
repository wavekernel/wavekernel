#
# Fortran macros for NAG's f95 on a Mac
# Special steps taken to work around case-insensitivity of the file system
#
# This macro file defines the symbol __NAG__ for conditional compilation
#
.SUFFIXES: .SUFFIXES .f .f90 .f95 .F .F90 .F95 
#
LIB_STD=$(FLIB_ROOT)/lib/
MOD_STD=$(FLIB_ROOT)/modules/
INC_STD=$(FLIB_ROOT)/include/
BIN_STD=$(FLIB_ROOT)/bin/
#
FC=f95 -u
F77=f95 -u
CFLAGS= 
# Could use also -dcfuns -dusty
FFLAGS= -O2  $(CFLAGS)       
F77_FLAGS= -O 
FFLAGS_DEBUG= -g  
FFLAGS_CHECK= -g -mtrace=all -C
FFLAGS_PROFILE= -pg
LDFLAGS=      $(CFLAGS)
#
INC_PREFIX=-I
MOD_PREFIX=-I
LIB_PREFIX=-L
#
MOD_EXT=.mod
MOD_SEARCH_STD= $(MOD_PREFIX)$(MOD_STD) 
MOD_SEARCH= $(MOD_SEARCH_STD) $(MOD_SEARCH_OTHER)
#INC_SEARCH= $(INC_PREFIX)$(INC_STD)
#
#
AR=ar
RANLIB=ranlib
#
CPP=/usr/local/lib/NAGWare/fpp -P
#
COCO=$(BIN_STD)coco -I $(INC_STD)
#DEFS=-D__F__
DEFS=-D__NAG__
#
# Experimental : the following deactivates an implicit rule
# which breaks havoc with the operation of this makefile
# It works at least with GNU make
%.o : %.mod
#
.F.o:
	$(CPP) -fixed $(DEFS) $*.F > aux_$*.f
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)  $(DEFS) aux_$*.f
	@rm -f aux_$*.f
	@mv aux_$*.o $*.o
.f.o:
	$(FC) -c $(FFLAGS)   $<
.F90.o:
	$(CPP) -free $(DEFS) $*.F90 > aux_$*.f90
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(DEFS) aux_$*.f90
	@rm -f aux_$*.f90
	@mv aux_$*.o $*.o
.f90.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
.F95.o:
	$(CPP) -free $(DEFS) $*.F95 > aux_$*.f95
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS) $(DEFS) aux_$*.f95
	@rm -f aux_$*.f95
	@mv aux_$*.o $*.o
.f95.o:
	$(FC) -c $(MOD_SEARCH) $(INC_SEARCH) $(FFLAGS)   $<
#










