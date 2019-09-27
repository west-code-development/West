# Makefile for Wstat

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Libraries/Forpy \
          $(MOD_FLAG)../Libraries/Json \
          $(MOD_FLAG)../Libraries/Base64 \
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


WSTAT_OBJS = \
wstat_readin.o \
wstat_setup.o \
wstat_restart.o \
wstat_memory_report.o \
wstat_tools.o \
davidson_diago.o \
wstat.o

PWOBJS = ../../PW/src/libpw.a 
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../Libraries/Forpy/libforpy.a ../Libraries/Json/libjson.a ../Libraries/Base64/libbase64.a 
LIBS += `python3-config --ldflags`
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

all : title tldeps wstat.x

wstat.x : $(PWOBJS) $(WSTAT_OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o wstat.x \
		$(WSTAT_OBJS) $(WESTLIBS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../West/Wstat/wstat.x . )
tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/wstat.x

title :
	@echo "              "
	@echo "##############"
	@echo "### Wstat ####"
	@echo "##############"
	@echo "              "

include make.depend
