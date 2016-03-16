#================================================================
# ELSES version 0.05
# Copyright (C) ELSES. 2007-2014 all rights reserved
#================================================================
#
# Makefile for ELSES
#

TOPDIR=..

include $(TOPDIR)/Makefile.inc

#TARGET = elses-xml-selftest elses
TARGET = elses

.PHONY: all
all: $(TARGET)

F90 ?= $(FC)
F77 ?= $(FC)

FFLAGS += $(FFLAGS-NRL) $(FFLAGS-XML) -I$(TOPDIR)/ext-lib

.SUFFIXES: .f90 .f

.f90.o:
	$(F90) -c $(FFLAGS) $< -o $@

.f.o:
	$(F77) -c $(FFLAGS) $< -o $@


OBJS-XML = \
	elses-lib-dst-info.o  elses-lib-wall-clock-time.o \
	elses-xml-config.o elses-xml-misc.o elses-options.o elses-config.o elses-lib-mpi-wrapper-compile.o \
	elses-lib-quad-precision.o \
	elses-xml-sax-counter.o elses-xml-sax-handler.o  elses-xml-sax-parser.o elses-xml-compat-check.o

OBJS-NRL = \
	Nrl/nrlvar.o Nrl/nrladapt.o Nrl/nrlsystem.o Nrl/nrlpair.o \
	Nrl/sktable.o Nrl/skderiv.o Nrl/nrlcut.o Nrl/nrlhop.o \
	Nrl/nrlonsite.o Nrl/nrlmatrix.o Nrl/nrlio.o Nrl/nrlldos.o \
	Nrl/nrlrep.o

OBJS-ELSES = \
	elses-lib-phys-const.o elses-xml-get-value-with-unit.o  elses-lib-element-name.o \
	elses-lib-db-ion-ene.o elses-lib-db-polarization.o elses-lib-db-covalent-radius.o \
	elses-mod-01.o elses-mod-02.o elses-mod-03.o elses-mod-04.o \
	elses-lib-math-func.o elses-lib-sort.o \
	elses-mod-file-io.o elses-mod-nrl.o elses-md-motion.o elses-md-booking.o elses-md-dst-cell.o elses-md-cell-change.o elses-md-constraint.o elses-md-motion-verlet.o \
	elses-qm-geno-rest-kernal.o \
	elses-qm-domain.o elses-qm-engine.o elses-qm-engine-nrl.o elses-qm-geno-output.o  elses-qm-engine-csc.o elses-qm-engine-csc-tb0.o \
	elses-qm-geno.o elses-ini-load-common.o elses-ini-load-geno.o elses-io-ctrl-dos-plot.o elses-io-ctrl-plot-eigenstates.o \
	elses-qm-vdW-parameters.o elses-ini-dst.o \
	elses-qm-geno-Huckel-atom-parameters.o elses-xml-check-element-info.o \
	elses-qm-geno-Huckel.o elses-qm-geno-csc-basic.o elses-qm-geno-csc.o  elses-qm-solver-eig-geno.o \
	elses-00-main.o elses-md-main.o elses-md-init.o elses-md-velocity-routines.o \
	elses-qm-main.o elses-kr-main.o \
	elses-kr-sub.o elses-kr-setkry.o  elses-hami-01.o \
	elses-lib-random-num.o elses-lib-omp.o \
	elses-xml-02.o elses-xml-03.o elses-mod-eig-leg.o \
	elses-qm-nrl.o elses-eig-leg.o elses-eig-leg2.o \
	elses-hami-gaas-molteni.o \
	elses-ana-plot-hami.o \
	elses-vis-main.o elses-vis-rasmol.o elses-ana-icohp.o elses-md-save-struct.o elses-md-output.o \
	elses-md-output-atom-energy.o elses-md-save-restart-xml.o \
	elses-vis-xyz-extended.o elses-00-version-info.o \
	elses-qm-solver-geno-main.o \
	elses-qm-solver-gscocg.o elses-qm-projection.o elses-la-matvec-routines.o \
	elses-la-matvec-io.o elses-la-gscocg-cg-s-mat.o \
	elses-la-gscocg-main.o elses-la-gscocg-tru-res.o \
	elses-la-gkrylov-hist.o elses-qm-solver-krgl.o elses-la-krgl-main.o elses-la-gkrylov-output.o \
	elses-band-calc.o elses-qm-output-eigenstates.o elses-qm-output-atom-charge.o elses-qm-output-basis-info.o \
	elses-qm-output-icohp.o elses-qm-output-wfn-charge.o \
	elses-la-gkrylov-main.o elses-qm-solver-gkrylov.o elses-md-dst.o \
	elses-qm-solver-gkrylov-dst.o elses-qm-engine-dst.o elses-qm-dst-global-weight.o \
	elses-qm-dst-global-dens-mat.o elses-md-dst-booking.o elses-qm-dst-domain.o \
	elses-qm-geno-dst.o elses-qm-dst-global-ham-mat.o elses-la-matvec-dst.o elses-qm-geno-csc-dst.o \
	elses-qm-dst-projection-cell.o elses-qm-population.o elses-io-dst-write-log.o elses-io-dst-system-load.o \
	elses-qm-dst-micro-mat.o elses-md-get-distance.o elses-qm-charge-from-weight-dst.o \
	elses-la-matvec-routines-dstm.o elses-qm-partial-trace-dstm.o  elses-la-gkrylov-main-dstm.o \
	elses-qm-output-matrices.o elses-qm-output-levels.o elses-la-eigen-solver.o \
	elses-qm-wavepacket-main.o \
	elses-md-make-bk-list-w-cell.o elses-md-dst-get-atom-range.o elses-qm-geno-rest-dstm.o \
	elses-md-dst-update-cell-booking.o elses-qm-dstm-ini-setting.o \
	elses-qm-alloc-charge.o elses-qm-projection-list-dstm.o elses-io-read-xyz-file.o \
	elses-lib-stop-signal.o elses-qm-projection-suggest.o elses-la-eigen-tridiag.o \
	elses-qm-cohp-dstm.o elses-qm-cohp-dstm-plot.o elses-io-update-snapshot.o elses-io-pair-list-for-cohp.o \
	elses-qm-set-s-as-i.o elses-qm-output-participation.o elses-qm-output-rest-energy-pair.o elses-eig-standard.o \
	elses-lib-report-mpi-time.o elses-qm-population-ortho.o \
	elses-qm-output-atom-charge-ortho.o elses-xml-sax-data-sync.o elses-qm-flexible-cutoff.o \
	elses-qm-dstm-get-nna-distance.o elses-md-alloc.o elses-ini-set-str-from-vatom.o elses-ini-set-str-from-matom.o \
	elses-xml-sax-data-sync-matom.o elses-qm-chk-allocated.o elses-md-virial-pressure.o elses-la-eigen-solver-main.o \
	elses-gid-initial.o elses-gid-basic-routines.o elses-gid-constraint.o elses-lib-read-w-comment.o elses-lib-read-w-split.o \
	elses-ext-atom-info.o elses-ext-matrix-data.o elses-ext-test-quantum-dynamics.o \
	elses-la-matvec-dstm-crs.o elses-la-matvec-dstm-dens.o elses-la-matvec-dstm-bcrs.o elses-la-matvec-dstm-sss.o elses-qm-dstm-reordering.o

#elses-xml-selftest: elses-xml-selftest.o $(OBJS-XML)
#	$(F90)  $(LDFLAGS) elses-xml-selftest.o $(OBJS-XML) $(LIBS-XML) -o $@
#	$(XML_SELFTEST) selftest/config.xml || exit

elses: $(OBJS-XML) $(OBJS-NRL) $(OBJS-ELSES)
	$(F90) $(LDFLAGS) $(OBJS-XML) $(OBJS-NRL) $(OBJS-ELSES) $(LIBS-XML) $(LIBS-EXT) $(LIBS) -o $@
	cp $@ $(BINDIR)

.PHONY: clean dep etags ctree
clean:
	@rm -f *.exe *.o *.mod *~ Nrl/*.o Nrl/*.mod Nrl/*~ $(TARGET)
	@rm -f $(BINDIR)/elses*

dep:
	makedepf90 -nosrc *.{f,f90} Nrl/*.{f,f90} > Makefile.dep

etags:
	etags *.{f,f90} Nrl/*.{f,f90}

ctree:
	@rm -f HTML/*
	@xmlfiles=`find $(TOPDIR)/$(XML-DIR) ! -path '*/Examples/*' ! -path '*/Tutorial/*' \( -name '*.f' -o -name '*.f90' \) -print`;\
     fortspec -d HTML *.{f,f90} Nrl/*.{f,f90} $$xmlfiles

include Makefile.dep
