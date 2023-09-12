# Makefile for the WEST software

include ../make.inc

default:
	@echo "Welcome to WEST!"
	@echo
	@echo "to install WEST, type at the shell prompt:"
	@echo '  make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"'
	@echo "  make [-j] target"
	@echo
	@echo "where target identifies one or multiple packages"
	@echo "  wstat        calculation of static dielectric response using PDEP"
	@echo "  wfreq        calculation of dynamical dielectric response and GW self-energy"
	@echo "  westpp       postprocessing programs"
	@echo "  wbse         calculation of BSE"
	@echo '  all          same as "make wstat wfreq westpp wbse"'
	@echo
	@echo "where target is one of the following operations:"
	@echo "  doc          build WEST documentation"
	@echo
	@echo "  clean        remove executables and objects"
	@echo '  veryclean    remove files produced by "make conf" as well'

conf: \
update_make_inc \
include_make_inc \
update_version \
depend

update_make_inc:
	@[ "${PYT}" ] || ( echo ">> PYT is not set. Ex: make conf PYT=python3"; exit 1 )
	@echo "Welcome to WEST!"
	@echo
	@echo "version : `${PYT} ./Pytools/read_json.py VERSION.json version`"
	@echo "url : `${PYT} ./Pytools/read_json.py VERSION.json url`"
	@echo "license : `${PYT} ./Pytools/read_json.py VERSION.json license`"
	@echo " " > west_make.inc
	@echo WESTDIR=`pwd` >> west_make.inc
	@echo PYT=${PYT} >> west_make.inc
	@echo PYT_LDFLAGS=${PYT_LDFLAGS} >> west_make.inc
	@echo
	@echo "Generated file: west_make.inc"
	@cat west_make.inc
	@echo

include_make_inc:
	@[ -f "west_make.inc" ] || ( echo '>> Cannot find west_make.inc. Run "make conf" first.'; exit 1 )
	$(eval include ./west_make.inc)

update_version:
	if test -d Modules ; then \
	( cd Modules ; ./update_west_version ${WESTDIR} `${PYT} ../Pytools/read_json.py ../VERSION.json version`; ) ; fi
	@echo "Generated file: west_version.f90"
	@echo

depend:
	if test -x west_makedeps ; then ./west_makedeps ; fi
	@echo

report_build_vars:
	@echo
	@echo "##############"
	@echo "# Build vars #"
	@echo "##############"
	@echo
	@[ "${MPIF90}" ] || ( echo ">> MPIF90 is not set."; exit 1 )
	@[ "${CC}" ] || ( echo ">> CC is not set."; exit 1 )
	@echo "# WEST_VERSION_NUMBER : `${PYT} ./Pytools/read_json.py VERSION.json version`"
	@echo "# WESTDIR : ${WESTDIR}"
	@echo "# FDFLAGS : ${FDFLAGS}"
	@echo "# IFLAGS : ${IFLAGS}"
	@echo "# MOD_FLAG : ${MOD_FLAG}"
	@echo "# MPIF90 : ${MPIF90}"
	@echo "# CC : ${CC}"
	@echo "# F77 : ${F77}"
	@echo "# CPP : ${CPP}"
	@echo "# CPPFLAGS : ${CPPFLAGS}"
	@echo "# CFLAGS : ${CFLAGS}"
	@echo "# F90FLAGS : ${F90FLAGS}"
	@echo "# FFLAGS : ${FFLAGS}"
	@echo "# FFLAGS_NOOPT : ${FFLAGS_NOOPT}"
	@echo "# LD : ${LD}"
	@echo "# LDFLAGS : ${LDFLAGS}"
	@echo "# LD_LIBS : ${LD_LIBS}"
	@echo "# BLAS_LIBS : ${BLAS_LIBS}"
	@echo "# LAPACK_LIBS : ${LAPACK_LIBS}"
	@echo "# SCALAPACK_LIBS : ${SCALAPACK_LIBS}"
	@echo "# FFT_LIBS : ${FFT_LIBS}"
	@echo "# MPI_LIBS : ${MPI_LIBS}"
	@echo "# MASS_LIBS : ${MASS_LIBS}"
	@echo "# AR : ${AR}"
	@echo "# ARFLAGS : ${ARFLAGS}"
	@echo "# RANLIB : ${RANLIB}"
	@echo "# FLIB_TARGETS : ${FLIB_TARGETS}"
	@echo "# WGET : ${WGET}"
	@echo "# PYT : ${PYT}"
	@echo "# PYT_LDFLAGS : ${PYT_LDFLAGS}"
	@echo


pytools: \
include_make_inc \
report_build_vars \
pytools_do

wstat: \
pytools \
libraries_do \
modules_do \
tools_do \
fft_kernel_do \
coulomb_kernel_do \
para_kernel_do \
hamiltonian_kernel_do \
dfpt_kernel_do \
io_kernel_do \
wstat_do

wfreq: \
pytools \
wstat \
wfreq_do

westpp: \
pytools \
wstat \
wbse \
westpp_do

wbse: \
pytools \
wstat \
wbse_do

all: \
pytools \
wstat \
wfreq \
wbse \
westpp

doc: include_make_inc
	if test -d doc ; then \
	( cd doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) html; \
	else $(MAKE) $(MFLAGS) html ; fi ) ; fi
	@echo "Open the file: ${WESTDIR}/doc/_build/html/index.html"

pytools_do:
	if test -d Pytools ; then \
	( cd Pytools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

libraries_do:
	if test -d Libraries ; then \
	( cd Libraries ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

modules_do: update_version libraries_do
	if test -d Modules ; then \
	( cd Modules ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

tools_do: modules_do libraries_do
	if test -d Tools ; then \
	( cd Tools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

fft_kernel_do: modules_do
	if test -d FFT_kernel ; then \
	( cd FFT_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

coulomb_kernel_do: tools_do modules_do
	if test -d Coulomb_kernel ; then \
	( cd Coulomb_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

para_kernel_do: tools_do
	if test -d Para_kernel ; then \
	( cd Para_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

hamiltonian_kernel_do: modules_do
	if test -d Hamiltonian_kernel ; then \
	( cd Hamiltonian_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

dfpt_kernel_do: para_kernel_do fft_kernel_do tools_do modules_do
	if test -d DFPT_kernel ; then \
	( cd DFPT_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

io_kernel_do: para_kernel_do tools_do fft_kernel_do modules_do libraries_do
	if test -d IO_kernel ; then \
	( cd IO_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

wstat_do: io_kernel_do dfpt_kernel_do para_kernel_do coulomb_kernel_do fft_kernel_do tools_do modules_do libraries_do
	if test -d Wstat ; then \
	( cd Wstat ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; \
	else $(MAKE) $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; fi ) ; fi

wfreq_do: io_kernel_do dfpt_kernel_do para_kernel_do coulomb_kernel_do fft_kernel_do tools_do modules_do libraries_do
	if test -d Wfreq ; then \
	( cd Wfreq ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; \
	else $(MAKE) $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; fi ) ; fi

wbse_do: wstat_do io_kernel_do dfpt_kernel_do para_kernel_do coulomb_kernel_do fft_kernel_do tools_do modules_do libraries_do
	if test -d ../LR_Modules; then \
	( cd ../LR_Modules ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; \
	else $(MAKE) $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; fi ) ; fi
	cd ../West
	if test -d Wbse ; then \
	( cd Wbse ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; \
	else $(MAKE) $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; fi ) ; fi

westpp_do: wbse_do io_kernel_do para_kernel_do coulomb_kernel_do fft_kernel_do tools_do modules_do libraries_do
	if test -d Westpp ; then \
	( cd Westpp ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; \
	else $(MAKE) $(MFLAGS) all PYT_LDFLAGS="${PYT_LDFLAGS}"; fi ) ; fi

clean: \
pytools_undo \
libraries_undo \
modules_undo \
tools_undo \
fft_kernel_undo \
coulomb_kernel_undo \
para_kernel_undo \
hamiltonian_kernel_undo \
dfpt_kernel_undo \
io_kernel_undo \
wstat_undo \
wfreq_undo \
wbse_undo \
westpp_undo

pytools_undo:
	if test -d Pytools ; then \
	( cd Pytools ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

libraries_undo:
	if test -d Libraries ; then \
	( cd Libraries ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

modules_undo:
	if test -d Modules ; then \
	( cd Modules ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

tools_undo:
	if test -d Tools ; then \
	( cd Tools ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

fft_kernel_undo:
	if test -d FFT_kernel ; then \
	( cd FFT_kernel ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

coulomb_kernel_undo:
	if test -d Coulomb_kernel ; then \
	( cd Coulomb_kernel ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

para_kernel_undo:
	if test -d Para_kernel ; then \
	( cd Para_kernel ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

hamiltonian_kernel_undo:
	if test -d Hamiltonian_kernel ; then \
	( cd Hamiltonian_kernel ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

dfpt_kernel_undo:
	if test -d DFPT_kernel ; then \
	( cd DFPT_kernel ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

io_kernel_undo:
	if test -d IO_kernel ; then \
	( cd IO_kernel ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

wstat_undo:
	if test -d Wstat ; then \
	( cd Wstat ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

wfreq_undo:
	if test -d Wfreq ; then \
	( cd Wfreq ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

wbse_undo:
	if test -d Wbse ; then \
	( cd Wbse ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

westpp_undo:
	if test -d Westpp ; then \
	( cd Westpp ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

unconf:
	[ -f "west_make.inc" ] && ( rm west_make.inc )

undepend:
	find . -name make.depend -delete

veryclean: clean \
unconf \
undepend
