# Makefile for Coulomb_kernel

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../Tools \
          $(MOD_FLAG).
IFLAGS=


COULOMB_KERNEL_OBJS = \
store_sqvc.o 

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
WESTLIBS = ../Tools/libtools.a

TLDEPS= bindir mods libs pw

all : tldeps lib_coulomb_kernel.a

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_coulomb_kernel.a : $(COULOMB_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f lib_coulomb_kernel.a

title :
	@echo "                       "
	@echo "#######################"
	@echo "### Coulomb_kernel ####"
	@echo "#######################"
	@echo "                       "

include make.depend
