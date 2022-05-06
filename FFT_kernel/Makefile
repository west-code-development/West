# Makefile for FFT_kernel

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG).

FFT_KERNEL_OBJS = \
fft_at_gamma.o \
fft_at_k.o \
fft_interpolation.o

TLDEPS = bindir mods pwlibs pw

all : title tldeps lib_fft_kernel.a

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_fft_kernel.a : $(FFT_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f lib_fft_kernel.a

title :
	@echo
	@echo "###################"
	@echo "### FFT_kernel ####"
	@echo "###################"
	@echo

include make.depend
