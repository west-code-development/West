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
           $(MOD_FLAG)../Libraries/Json \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG)../Para_kernel \
           $(MOD_FLAG)../Tools \
           $(MOD_FLAG).

WBSEPP_OBJS = \
wbsepp_decompose_eig_contributions.o \
wbsepp_ads_spectrum.o \
wbsepp_plot_exc.o \
wbsepp_plot_charged_density_res_exc.o \
wbsepp.o

PWOBJS = ../../PW/src/libpw.a

QEMODS = ../../KS_Solvers/libks_solvers.a $(BASEMODS)

WESTLIBS = \
../Wstat/libwstat.a \
../Wbse/libwbse.a \
../Westpp/libwestpp.a \
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

all : title tldeps wbsepp.x

wbsepp.x : $(PWOBJS) $(WBSEPP_OBJS) $(QEMODS) $(WESTLIBS)
	$(LD) $(LDFLAGS) -o wbsepp.x \
	        $(WBSEPP_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS) $(QELIBS) $(PYT_LDFLAGS)
	- ( cd ../../bin ; ln -fs ../West/Wbsepp/wbsepp.x . )

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/wbsepp.x

title :
	@echo
	@echo "##############"
	@echo "### Wbsepp ####"
	@echo "##############"
	@echo

include make.depend
