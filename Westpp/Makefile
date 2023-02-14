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
do_dip.o \
do_eigenpot2.o \
do_exc.o \
do_rho.o \
do_sxx.o \
do_wann.o \
do_wfc2.o \
do_loc.o \
do_resp.o \
dump_r.o \
westpp.o \
westpp_readin.o \
westpp_setup.o \
write_wfc_1d_r.o \
write_wfc_spav.o

PWOBJS = ../../PW/src/libpw.a

QEMODS = ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a $(BASEMODS)

WESTLIBS = \
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

TLDEPS = bindir mods pwlibs pw

all : title tldeps westpp.x

westpp.x : $(WESTPP_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(QELIBS) $(PYT_LDFLAGS)
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
