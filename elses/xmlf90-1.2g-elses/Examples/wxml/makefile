#
# Makefile for WXML examples
#
default: all
all:  simple pseudo
#
#---------------------------
MK=$(FLIB_ROOT)/fortran.mk
include $(MK)
#---------------------------
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
LIBS=$(LIB_PREFIX)$(LIB_STD) -lflib
#
simple: simple.o
	$(FC) $(LDFLAGS) -o simple simple.o  $(LIBS)
#
pseudo: m_pseudo_utils.o pseudo.o
	$(FC) $(LDFLAGS) -o pseudo m_pseudo_utils.o pseudo.o  $(LIBS)
#
clean: 
	rm -f simple pseudo  *.o *.$(MOD_EXT)
#










