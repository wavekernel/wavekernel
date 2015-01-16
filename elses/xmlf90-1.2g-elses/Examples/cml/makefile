#
# Makefile for CML examples
#
default: all
all: example
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
example: example.o
	$(FC) $(LDFLAGS) -o example example.o  $(LIBS)
#
clean: 
	rm -f example *.o *.$(MOD_EXT)
#










