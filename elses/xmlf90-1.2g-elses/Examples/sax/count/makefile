#
# Makefile for count (count elements in an xml document)
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
OBJS= m_count.o count.o 
# 
count: $(OBJS)
	$(FC) $(LDFLAGS) -o $@  $(OBJS)  $(LIBS)
#
clean: 
	rm -f *.o count *.$(MOD_EXT)
#










