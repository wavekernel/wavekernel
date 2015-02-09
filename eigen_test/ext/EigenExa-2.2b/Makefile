INSTALL_DIR = $(HOME)/lib
#INSTALL_DIR = $(HOME)/lib64
#INSTALL_DIR =

#
# Inclusion of architecture-depended rules
# select an appropriate make_inc file
#
#include make_incs/make_inc.gcc
#include make_incs/make_inc.Intel
#include make_incs/make_inc.Intel.shared
#include make_incs/make_inc.BX900
#include make_incs/make_inc.K_FX10
#include make_incs/make_inc.K_FX10.shared
#include make_incs/make_inc.BlueGeneQ
include make_incs/make_inc


#
# Definition of source objects
#
OBJS_TRD  = \
	eigen_trd.o \
	eigen_trd_t1.o eigen_trd_t2.o eigen_trd_t4.o \
	eigen_trd_t5.o eigen_trd_t5x.o \
	eigen_trd_t6_3.o eigen_trd_t7.o eigen_trd_t8.o \
	eigen_prd.o \
	eigen_prd_t1.o eigen_prd_t2.o eigen_prd_t4x.o \
	eigen_prd_t5.o \
	eigen_prd_t6_3.o eigen_prd_t7.o eigen_prd_t8.o
OBJS_DCX  = \
	dcx.o dc2.o bisect.o bisect2.o \
	mx_pdstedc.o mx_pdlaedz.o mx_pdlaed0.o mx_pdlaed1.o \
	mx_pdlaed2.o mx_pdlaed3.o \
	my_pdsxedc.o my_pdlaedz.o my_pdlaed0.o my_pdlaed1.o \
	my_pdlaed2.o my_pdlaed3.o my_pdlasrt.o lapack_eigen.o \
	dlaed3.o dlaed6.o dlaed6_init.o \
	dc_redist1.o dc_redist2.o dlaset.o dlacpy.o
OBJS_TBK  = \
	trbakwy4.o trbakwy4_body.o
OBJS_ABC  = \
	eigen_sx.o \
	eigen_s.o \
	eigen_scaling.o \
	comm.o matrix_adjustment.o \
	CSTAB.o miscC.o
OBJS_MOD  = \
	eigen_libs.o eigen_blacs.o eigen_devel.o \
	eigen_house.o eigen_dc.o eigen_trbak.o \
	dlaed6_common.o

OBJS      = $(OBJS_TRD) $(OBJS_DCX) $(OBJS_TBK) $(OBJS_ABC) $(OBJS_MOD) \
		$(EXTRAOBJS)

OBJS_TEST  = main2.o mat_set.o ev_test.o w_test.o
OBJS_TEST1 = samples/main_comm_create.o mat_set.o ev_test.o
OBJS_TEST2 = samples/main_comm_create_null.o mat_set.o ev_test.o
OBJS_TEST3 = samples/main_cartesian.o mat_set.o ev_test.o

MODEIGEN_PUB  = \
	eigen_blacs.mod eigen_libs.mod
MODEIGEN_SYS  = \
	eigen_devel.mod \
	eigen_house.mod eigen_dc.mod eigen_trbak.mod \
	dlaed6_common.mod

MODEIGEN  = $(MODEIGEN_PUB) $(MODEIGEN_SYS)

ifdef SHARED
LIBEIGEN  = libEigenExa.so
else
LIBEIGEN  = libEigenExa.a
endif
LIBEIGEN2 = libEigenExa.a libEigenExa.so
TESTEIGEN = eigenexa_benchmark a.out comm_create.exe comm_null.exe cart.exe

LD_EIGEN  = -L./ -lEigenExa

# MACROS
AT_SYNC   = \
	$(MACRO_D_PREFIX)USE_SYNC=0 \
	$(MACRO_D_PREFIX)_DEBUG_=0 \
	$(MACRO_D_PREFIX)_DEBUG2_=0
AT_BCAST_OVERLAP = \
	$(MACRO_D_PREFIX)AT_BCAST_OVERLAP=1 \
	$(MACRO_D_PREFIX)OVERLAP_DECISION_TYPE=3
DEBUGFLAG = \
	$(MACRO_D_PREFIX)TIMER_PRINT=1
#DEBUGFLAG = \
#	$(MACRO_D_PREFIX)TIMER_PRINT=0

BASIC_CCFLAG = $(CCFLAG) $(INCFLAG) $(MODFLAG) $(CPPFLAG) $(DEBUGFLAG)
NOOPT_CCFLAG = $(CCFLAG) $(INCFLAG) $(MODFLAG) $(CPPFLAG) $(DEBUGFLAG)
BASIC_LDFLAG = $(LD_EIGEN) $(LDFLAG) $(LD_EIGEN) $(LDFLAG)


#
# Start Make rule
#
all : $(LIBEIGEN)
test : $(TESTEIGEN)


#
# Suffix rules
#
.o.mod:
	@if [ ! -e $(@:.mod=.F) ]; then \
		exit 1; \
	else \
		if [ ! -e $@ ]; then \
			\rm $<; make $<; \
		fi \
	fi
.F.o:
	@if [ -e $(<:.F=.mod) ]; then \
		\rm $(<:.F=.mod); \
	fi
	$(FF) -o $@ -c $< $(BASIC_CCFLAG)
	@if [ -e $(<:.F=.mod) ]; then \
		touch $(<:.F=.mod); \
	fi
.f.o:
	$(FF) -o $@ -c $< $(CCFLAG) $(INCFLAG) $(MODFLAG) $(DEBUGFLAG)
.c.o:
	$(CC) -o $@ -c $< $(CCFLAG) $(INCFLAG)


#
# libEigenExa
#
$(LIBEIGEN) : $(OBJS)
ifdef SHARED
	$(CC) $(SHARED_LDFLAG) -o $(LIBEIGEN) $(OBJS) $(SHARED_LIBEIGEN_LDFLAG)
else
	$(AR) r $(LIBEIGEN) $(OBJS)
endif

#
# benchmark and test codes
#
eigenexa_benchmark : $(OBJS_TEST) $(LIBEIGEN)
	$(FF) -o $@ $(OBJS_TEST) $(BASIC_LDFLAG)
a.out : eigenexa_benchmark
	@if [ -e $@ ]; then \
		\rm $@; \
	fi
	@ln -s $< $@
comm_create.exe : $(OBJS_TEST1) $(LIBEIGEN)
	$(FF) -o $@ $(OBJS_TEST1) $(BASIC_LDFLAG)
comm_null.exe : $(OBJS_TEST2) $(LIBEIGEN)
	$(FF) -o $@ $(OBJS_TEST2) $(BASIC_LDFLAG)
cart.exe : $(OBJS_TEST3) $(LIBEIGEN)
	$(FF) -o $@ $(OBJS_TEST3) $(BASIC_LDFLAG)


#
# special rures
#
eigen_libs.o : eigen_libs.F
	$(FF) -o $@ -c $< $(NOOPT_CCFLAG)
eigen_devel.o : eigen_devel.F
	$(FF) -o $@ -c $< $(NOOPT_CCFLAG)
comm.o : comm.F
	$(FF) -o $@ -c $< $(NOOPT_CCFLAG)
eigen_trd.o : eigen_trd.F
	$(FF) -o $@ -c $< $(BASIC_CCFLAG) $(AT_SYNC)
eigen_prd.o : eigen_prd.F
	$(FF) -o $@ -c $< $(BASIC_CCFLAG) $(AT_SYNC)
eigen_trd_t2.o : eigen_trd_t2.F
	$(FF) -o $@ -c $< $(BASIC_CCFLAG) $(AT_SYNC)
eigen_prd_t2.o : eigen_prd_t2.F
	$(FF) -o $@ -c $< $(BASIC_CCFLAG) $(AT_SYNC)
trbakwy4.o : trbakwy4.F
	$(FF) -o $@ -c $< $(BASIC_CCFLAG) $(AT_BCAST_OVERLAP)
trbakwy4_body.o : trbakwy4_body.F
	$(FF) -o $@ -c $< $(BASIC_CCFLAG) $(AT_BCAST_OVERLAP)
mat_set.o : mat_set.f
	$(FF) -o $@ -c $< $(BASIC_CCFLAG) $(TEST_FLAG)
ev_test.o : ev_test.f
	$(FF) -o $@ -c $< $(BASIC_CCFLAG) $(TEST_FLAG)


#
# dependencies among .o and .mod
#
$(OBJS_TRD) $(OBJS_DCX) $(OBJS_TBK) $(OBJS_ABC) : $(MODEIGEN)
$(OBJS_TEST) : $(MODEIGEN)

#
# dependencies among .mod's
#
eigen_libs.mod : eigen_blacs.mod eigen_devel.mod
eigen_libs.o : eigen_blacs.o eigen_devel.o


#
# installation
#
install : all
	@if [ '$(INSTALL_DIR)X' == 'X' ]; then \
		echo 'Please setup INSTALL_DIR in Makefile'; \
	else \
		if [ ! -d $(INSTALL_DIR) ]; then \
			mkdir $(INSTALL_DIR); \
		fi; \
		for f in $(LIBEIGEN) $(MODEIGEN_PUB) ; do \
			echo 'install ./'$$f' $(INSTALL_DIR)/' ; \
			install ./$$f $(INSTALL_DIR)/; \
		done; \
	fi

#
# clean
#
clean :
	@echo 'rm *.o samples/*.o *.s *.lst core.* *.mod '$(LIBEIGEN2)' '$(TESTEIGEN)''
	@for f in \
		  *.o samples/*.o *.s *.lst core.* *.mod  $(LIBEIGEN2)   $(TESTEIGEN); do \
		if [ -e $$f ]; then \
		if [ -f $$f ]; then \
			\rm $$f; \
		fi; \
		fi; \
		if [ -h $$f ]; then \
			\rm $$f; \
		fi; \
	done
