# Makefile for the WEST software 

default: all

pytools: \
pytools_do

wstat: \
pytools_do \
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
wstat \
wfreq_do

westpp: \
wstat \
wfreq \
westpp_do

all: \
pytools \
wstat \
wfreq \
westpp 

pytools_do:
	if test -d Pytools ; then \
	( cd Pytools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

libraries_do:
	if test -d Libraries ; then \
	( cd Libraries ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

modules_do:
	if test -d Modules ; then \
	( cd Modules ; sh update_west_version; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

tools_do:
	if test -d Tools ; then \
	( cd Tools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

fft_kernel_do:
	if test -d FFT_kernel ; then \
	( cd FFT_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

coulomb_kernel_do:
	if test -d Coulomb_kernel ; then \
	( cd Coulomb_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

para_kernel_do:
	if test -d Para_kernel ; then \
	( cd Para_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

hamiltonian_kernel_do:
	if test -d Hamiltonian_kernel ; then \
	( cd Hamiltonian_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

dfpt_kernel_do:
	if test -d DFPT_kernel ; then \
	( cd DFPT_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

io_kernel_do:
	if test -d IO_kernel ; then \
	( cd IO_kernel ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

wstat_do:
	if test -d Wstat ; then \
	( cd Wstat ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

wfreq_do:
	if test -d Wfreq ; then \
	( cd Wfreq ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

westpp_do:
	if test -d Westpp ; then \
	( cd Westpp ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

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

westpp_undo:
	if test -d Westpp ; then \
	( cd Westpp ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

distclean: clean 
