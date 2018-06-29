# Makefile for FFT_kernel

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG).
IFLAGS=


FFT_KERNEL_OBJS = \
fft_at_gamma.o \
fft_at_k.o

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

TLDEPS= bindir mods libs pw

all : title tldeps lib_fft_kernel.a

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_fft_kernel.a : $(FFT_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f lib_fft_kernel.a

title :
	@echo "                   "
	@echo "###################"
	@echo "### FFT_kernel ####"
	@echo "###################"
	@echo "                   "

include make.depend
