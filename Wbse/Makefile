# Makefile for Wbse

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../Coulomb_kernel \
           $(MOD_FLAG)../DFPT_kernel \
           $(MOD_FLAG)../FFT_kernel \
           $(MOD_FLAG)../IO_kernel \
           $(MOD_FLAG)../Libraries/Base64 \
           $(MOD_FLAG)../Libraries/Forpy \
           $(MOD_FLAG)../Libraries/iotk \
           $(MOD_FLAG)../Libraries/Json \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG)../Para_kernel \
           $(MOD_FLAG)../Tools \
           $(MOD_FLAG).

WBSE_OBJS = \
bse_hqp_psi.o \
bse_localization.o\
bse_module.o \
bse_setup.o \
bse_tools.o \
check_ovl_wfc.o \
lanczos_db.o \
lanczos_restart.o \
plep_db.o \
plep_io.o \
td_liouville_oper.o \
wbse_bse_kernel.o \
wbse_calc_dens.o \
wbse_davidson_diago.o \
wbse_davidson_diago_square.o\
wbse_dot.o \
wbse_init_pdep.o \
wbse_init_qboxcoupling.o \
wbse_init_restart.o \
wbse_lanczos_diago.o \
wbse_lanczos_diago.o \
wbse_memory_report.o \
wbse_readin.o \
wbse_restart.o  \
wbse_setup.o \
wbse_solve_e_psi.o \
wbse_tools.o \
west_dv_of_drho.o \
west_dv_setup.o

PWOBJS = ../../PW/src/libpw.a

QEMODS = ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a $(BASEMODS)

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

all : title tldeps libwbse.a wbse.x wbse_init.x

libwbse.a : $(WBSE_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

wbse.x : wbse.o libwbse.a $(PWOBJS) $(WBSE_OBJS) $(QEMODS) $(WESTLIBS)
	$(LD) $(LDFLAGS) -o $@ wbse.o libwbse.a \
	        $(WESTLIBS) $(PWOBJS) $(QEMODS) $(QELIBS) $(PYT_LDFLAGS)
	- ( cd ../../bin ; ln -fs ../West/Wbse/wbse.x . )

wbse_init.x : wbse_init.o libwbse.a $(PWOBJS) $(WBSE_OBJS) $(QEMODS) $(WESTLIBS)
	$(LD) $(LDFLAGS) -o $@ wbse_init.o libwbse.a \
		$(WESTLIBS) $(PWOBJS) $(QEMODS) $(QELIBS) $(PYT_LDFLAGS)
	- ( cd ../../bin ; ln -fs ../West/Wbse/wbse_init.x . )

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/wbse.x ../../bin/wbse_init.x

title :
	@echo
	@echo "##############"
	@echo "### Wbse ####"
	@echo "##############"
	@echo

include make.depend
