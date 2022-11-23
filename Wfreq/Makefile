# Makefile for Wfreq

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../Coulomb_kernel \
           $(MOD_FLAG)../FFT_kernel \
           $(MOD_FLAG)../IO_kernel \
           $(MOD_FLAG)../Libraries/Base64 \
           $(MOD_FLAG)../Libraries/Json \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG)../Para_kernel \
           $(MOD_FLAG)../Tools \
           $(MOD_FLAG).

WFREQ_OBJS = \
calc_corr.o \
calc_exx2.o \
calc_vxc.o \
chi_invert.o \
diago_lanczos.o \
get_brak_hyper_parallel.o \
set_freqlists.o \
solve_deflated_lanczos_w_full_ortho.o \
solve_eri.o \
solve_gfreq.o \
solve_h1e.o \
solve_hf.o \
solve_qp.o \
solve_wfreq.o \
wfreq.o \
wfreq_memory_report.o \
wfreq_readin.o \
wfreq_restart.o \
wfreq_setup.o

PWOBJS = ../../PW/src/libpw.a

QEMODS = ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a $(BASEMODS)

WESTLIBS = \
../IO_kernel/lib_io_kernel.a \
../DFPT_kernel/lib_dfpt_kernel.a \
../Hamiltonian_kernel/lib_hamiltonian_kernel.a \
../Para_kernel/lib_para_kernel.a \
../Coulomb_kernel/lib_coulomb_kernel.a \
../FFT_kernel/lib_fft_kernel.a \
../Tools/libtools.a \
../Modules/libmodules.a \
../Libraries/Forpy/libforpy.a \
../Libraries/Json/libjson.a \
../Libraries/Base64/libbase64.a

TLDEPS = bindir mods pwlibs pw

all : title tldeps wfreq.x

wfreq.x : $(WFREQ_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(QELIBS) ${PYT_LDFLAGS}
	- ( cd ../../bin ; ln -fs ../West/Wfreq/wfreq.x . )

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/wfreq.x

title :
	@echo
	@echo "##############"
	@echo "### Wfreq ####"
	@echo "##############"
	@echo

include make.depend
