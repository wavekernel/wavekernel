#
# Makefile for Xpath examples
#
default: all
all: relative pseudo text
#
#---------------------------
MK=$(FLIB_ROOT)/fortran.mk
include $(MK)
#---------------------------
#
# Uncomment the following line for debugging support
#
FFLAGS=$(FFLAGS_DEBUG)
#
LIBS=$(LIB_PREFIX)$(LIB_STD) -lflib
#
relative:  relative.o
	$(FC) $(LDFLAGS) -o relative relative.o  $(LIBS)
pseudo:  m_pseudo_types.o pseudo.o 
	$(FC) $(LDFLAGS) -o pseudo m_pseudo_types.o pseudo.o  $(LIBS)
text:  text.o
	$(FC) $(LDFLAGS) -o text text.o  $(LIBS)
#
clean: 
	rm -f relative pseudo text *.o *.$(MOD_EXT)
#










