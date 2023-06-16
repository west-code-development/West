# Makefile for Wbse

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../../LR_Modules \
           $(MOD_FLAG)../Coulomb_kernel \
           $(MOD_FLAG)../DFPT_kernel \
           $(MOD_FLAG)../FFT_kernel \
           $(MOD_FLAG)../IO_kernel \
           $(MOD_FLAG)../Libraries/Base64 \
           $(MOD_FLAG)../Libraries/Forpy \
           $(MOD_FLAG)../Libraries/Json \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG)../Para_kernel \
           $(MOD_FLAG)../Tools \
           $(MOD_FLAG)../Wstat \
           $(MOD_FLAG).

WBSE_OBJS = \
bse_kernel.o \
build_rhs_zvector_eq.o \
calc_tau.o \
check_ovl_wfc.o \
hybrid_kernel.o \
lanczos_restart.o \
td_liouville_oper.o \
wbse_calc_dens.o \
wbse_calc_forces.o \
wbse_davidson_diago.o \
wbse_dot.o \
wbse_dv.o \
wbse_init_readin.o \
wbse_init_restart.o \
wbse_init_setup.o \
wbse_lanczos_diago.o \
wbse_localization.o\
wbse_memory_report.o \
wbse_readin.o \
wbse_setup.o \
wbse_solve_e_psi.o \
wbse_tools.o

PWOBJS = ../../PW/src/libpw.a

QEMODS = ../../LR_Modules/liblrmod.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a $(BASEMODS)

WESTLIBS = \
../Wstat/libwstat.a \
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

all : title tldeps wbse.x wbse_init.x

wbse.x : wbse.o $(WBSE_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(QELIBS) $(PYT_LDFLAGS)
	- ( cd ../../bin ; ln -fs ../West/Wbse/wbse.x . )

wbse_init.x : wbse_init.o $(WBSE_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(QELIBS) $(PYT_LDFLAGS)
	- ( cd ../../bin ; ln -fs ../West/Wbse/wbse_init.x . )

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/wbse.x ../../bin/wbse_init.x

title :
	@echo
	@echo "#############"
	@echo "### Wbse ####"
	@echo "#############"
	@echo

include make.depend
