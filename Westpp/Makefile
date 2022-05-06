# Makefile for Westpp

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../Coulomb_kernel \
           $(MOD_FLAG)../FFT_kernel \
           $(MOD_FLAG)../IO_kernel \
           $(MOD_FLAG)../Libraries/Json \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG)../Para_kernel \
           $(MOD_FLAG)../Tools \
           $(MOD_FLAG).

WESTPP_OBJS = \
do_eigenpot2.o \
do_rho.o \
do_sxx.o \
do_wfc2.o \
dump_r.o \
westpp.o \
westpp_readin.o \
westpp_setup.o \
write_wfc_1d_r.o \
write_wfc_spav.o

PWOBJS = ../../PW/src/libpw.a

QEMODS = \
../../Modules/libqemod.a \
../../KS_Solvers/libks_solvers.a \
../../FFTXlib/libqefft.a \
../../upflib/libupf.a \
../../LAXlib/libqela.a \
../../UtilXlib/libutil.a \
../../XClib/xc_lib.a

WESTLIBS = \
../DFPT_kernel/lib_dfpt_kernel.a \
../IO_kernel/lib_io_kernel.a \
../Hamiltonian_kernel/lib_hamiltonian_kernel.a \
../Para_kernel/lib_para_kernel.a \
../Coulomb_kernel/lib_coulomb_kernel.a \
../FFT_kernel/lib_fft_kernel.a \
../Tools/libtools.a \
../Modules/libmodules.a \
../Libraries/Forpy/libforpy.a \
../Libraries/Json/libjson.a \
../Libraries/Base64/libbase64.a

TLDEPS = bindir mods libs pw

all : title tldeps westpp.x

westpp.x : $(PWOBJS) $(WESTPP_OBJS) $(QEMODS) $(WESTLIBS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o westpp.x \
		$(WESTPP_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS) $(LIBS) ${PYT_LDFLAGS}
	- ( cd ../../bin ; ln -fs ../West/Westpp/westpp.x . )

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/westpp.x

title :
	@echo
	@echo "##############"
	@echo "### Westpp ###"
	@echo "##############"
	@echo

include make.depend
