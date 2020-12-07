# Makefile for Wbse

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
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
          $(MOD_FLAG)../Wbse \
          $(MOD_FLAG)../Westpp \
          $(MOD_FLAG).

#wbsepp_fetch_namelist.o \
#wbseppcom.o \

IFLAGS=

WBSE_OBJS = \
wbsepp_decompose_eig_contributions.o \
wbsepp_ads_spectrum.o \
wbsepp_plot_exc.o \
wbsepp_plot_charged_density_res_exc.o \
wbsepp.o

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a  ../Libraries/Json/libjson.a ../Libraries/Base64/libbase64.a ../Libraries/Forpy/libforpy.a
WESTLIBS = \
../Wstat/libwstat.a \
../Wbse/libwbse.a \
../Westpp/libwestpp.a \
../IO_kernel/lib_io_kernel.a \
../DFPT_kernel/lib_dfpt_kernel.a \
../Tools/libtools.a \
../Hamiltonian_kernel/lib_hamiltonian_kernel.a \
../Para_kernel/lib_para_kernel.a \
../Coulomb_kernel/lib_coulomb_kernel.a \
../FFT_kernel/lib_fft_kernel.a \
../Modules/libmodules.a

TLDEPS= bindir mods libs pw

all : title tldeps wbsepp.x

wbsepp.x : $(PWOBJS) $(WBSE_OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o wbsepp.x \
	        $(WBSE_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(PYT_LDFLAGS)
	- ( cd ../../bin ; ln -fs ../West/Wbsepp/wbsepp.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/wbsepp_spectrum.x

title :
	@echo "              "
	@echo "##############"
	@echo "### Wbsepp ####"
	@echo "##############"
	@echo "              "

include make.depend
