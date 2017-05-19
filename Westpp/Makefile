# Makefile for Westpp

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
          $(MOD_FLAG)../IO_kernel \
          $(MOD_FLAG)../DFPT_kernel \
          $(MOD_FLAG).
IFLAGS=


WESTPP_OBJS = \
westpp_readin.o \
westpp_setup.o \
write_wfc_1d_r.o \
write_wfc_spav.o \
dump_r.o \
do_rho.o \
do_wfc2.o \
do_eigenpot2.o \
do_sxx.o \
westpp.o

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../Libraries/Json/libjson.a
WESTLIBS = \
../DFPT_kernel/lib_dfpt_kernel.a \
../IO_kernel/lib_io_kernel.a \
../Hamiltonian_kernel/lib_hamiltonian_kernel.a \
../Para_kernel/lib_para_kernel.a \
../Coulomb_kernel/lib_coulomb_kernel.a \
../FFT_kernel/lib_fft_kernel.a \
../Tools/libtools.a \
../Modules/libmodules.a


TLDEPS= bindir mods libs pw

all : title tldeps westpp.x

westpp.x : $(PWOBJS) $(WESTPP_OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o westpp.x \
		$(WESTPP_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../West/Westpp/westpp.x . )
tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/westpp.x

title :
	@echo "              "
	@echo "##############"
	@echo "### Westpp ###"
	@echo "##############"
	@echo "              "

include make.depend
