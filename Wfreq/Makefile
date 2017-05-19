# Makefile for Wfreq

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Libraries/Json \
          $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../Tools \
          $(MOD_FLAG)../FFT_kernel \
          $(MOD_FLAG)../Coulomb_kernel \
          $(MOD_FLAG)../Para_kernel \
          $(MOD_FLAG)../Hamiltonian_kernel \
          $(MOD_FLAG)../DFPT_kernel \
          $(MOD_FLAG)../IO_kernel \
          $(MOD_FLAG).
IFLAGS=


WFREQ_OBJS = \
wfreq_readin.o \
wfreq_setup.o \
io_eigenfreq.o \
diago_lanczos.o \
solve_deflated_lanczos_w_full_ortho.o \
get_brak_hyper_parallel.o \
solve_wfreq.o \
set_freqlists.o \
chi_invert.o \
solve_gfreq.o \
calc_vxc.o \
calc_exx2.o \
solve_hf.o \
calc_corr.o \
solve_qp.o \
wfreq.o

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../Libraries/Json/libjson.a
WESTLIBS = \
../IO_kernel/lib_io_kernel.a \
../DFPT_kernel/lib_dfpt_kernel.a \
../Hamiltonian_kernel/lib_hamiltonian_kernel.a \
../Para_kernel/lib_para_kernel.a \
../Coulomb_kernel/lib_coulomb_kernel.a \
../FFT_kernel/lib_fft_kernel.a \
../Tools/libtools.a \
../Modules/libmodules.a


TLDEPS= bindir mods libs pw

all : title tldeps wfreq.x

wfreq.x : $(PWOBJS) $(WFREQ_OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o wfreq.x \
		$(WFREQ_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../West/Wfreq/wfreq.x . )
tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/wfreq.x

title :
	@echo "              "
	@echo "##############"
	@echo "### Wfreq ####"
	@echo "##############"
	@echo "              "

include make.depend
