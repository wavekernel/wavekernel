#
# Makefile for Pseudo XML processing
#
default: pseudo
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
OBJS= m_pseudo_types.o m_pseudo.o pseudo.o
     
pseudo: $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)
#
clean: 
	rm -f *.o pseudo *.$(MOD_EXT)
#










