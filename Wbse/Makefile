# Makefile for Wbse

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../Libraries/Json \
          $(MOD_FLAG)../Libraries/Base64 \
          $(MOD_FLAG)../Libraries/Forpy \
          $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../Tools \
          $(MOD_FLAG)../FFT_kernel \
          $(MOD_FLAG)../Coulomb_kernel \
          $(MOD_FLAG)../Para_kernel \
          $(MOD_FLAG)../Hamiltonian_kernel \
          $(MOD_FLAG)../DFPT_kernel \
          $(MOD_FLAG)../IO_kernel \
          $(MOD_FLAG)../Wstat \
          $(MOD_FLAG).


#  wbse_fetch_namelist.o \
#  wbsecom.o \
IFLAGS=
WBSE_OBJS = \
bse_module.o \
wbse_setup.o \
bse_setup.o \
wbse_readin.o \
plep_io.o \
plep_db.o \
wbse_restart.o  \
wbse_memory_report.o \
wbse_tools.o \
wbse_davidson_diago.o \
wbse_calc_dens.o \
td_liouville_oper.o \
west_dv_of_drho.o \
west_dv_setup.o \
wbse_bse_kernel.o \
bse_tools.o \
bse_hqp_psi.o \
wbse_solve_e_psi.o \
wbse_dot.o \
lanzcos_db.o \
lanzcos_restart.o \
wbse_lanzcos_diago.o \
wbse_lanzcos_diago.o \
check_ovl_wfc.o \
bse_localization.o\
wbse_init_restart.o \
wbse_init_pdep.o \
wbse_davidson_diago_square.o\
wbse_init_qboxcoupling.o \
wbse_init_fock_energy.o 

PWOBJS = ../../PW/src/libpw.a 
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../LR_Modules/liblrmod.a ../Libraries/Json/libjson.a ../Libraries/Base64/libbase64.a ../Libraries/Forpy/libforpy.a  #../../iotk/src/libiotk.a
WESTLIBS = \
../Wstat/libwstat.a \
../IO_kernel/lib_io_kernel.a \
../DFPT_kernel/lib_dfpt_kernel.a \
../Tools/libtools.a \
../Hamiltonian_kernel/lib_hamiltonian_kernel.a \
../Para_kernel/lib_para_kernel.a \
../Coulomb_kernel/lib_coulomb_kernel.a \
../FFT_kernel/lib_fft_kernel.a \
../Modules/libmodules.a 


TLDEPS= bindir mods libs pw

all : title tldeps libwbse.a wbse.x wbse_init.x

libwbse.a : $(WBSE_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

wbse.x : libwbse.a wbse.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ wbse.o libwbse.a \
	        $(WESTLIBS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(PYT_LDFLAGS)
	- ( cd ../../bin ; ln -fs ../West/Wbse/$@ . )

wbse_init.x : wbse_init.o libwbse.a $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ wbse_init.o libwbse.a \
		$(WESTLIBS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(PYT_LDFLAGS)
	- ( cd ../../bin ; ln -fs ../West/Wbse/wbse_init.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/wbse.x ../../bin/wbse_init.x

title :
	@echo "              "
	@echo "##############"
	@echo "### Wbse ####"
	@echo "##############"
	@echo "              "

include make.depend
