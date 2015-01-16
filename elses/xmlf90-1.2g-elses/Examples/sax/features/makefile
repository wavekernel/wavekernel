#
# Makefile for example of XML processing with all kinds of handlers
# and features.
#
default: example  xmlcheck
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
OBJS= m_handlers.o example.o
     
example:  $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS)  $(LIBS)
xmlcheck: xmlcheck.o
	$(FC) $(LDFLAGS) -o $@ xmlcheck.o  $(LIBS)
#
clean: 
	rm -f *.o example xmlcheck *.$(MOD_EXT)
#










