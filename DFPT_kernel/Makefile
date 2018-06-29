# Makefile for DFPT_kernel

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../Tools \
          $(MOD_FLAG)../FFT_kernel \
          $(MOD_FLAG)../Coulomb_kernel \
          $(MOD_FLAG).
IFLAGS=


DFPT_KERNEL_OBJS = \
apply_sternheimerop_to_m_wfcs.o \
precondition_m_wfcts.o \
linsolve_sternheimer_m_wfcts.o \
dfpt_module.o  

PWOBJS = ../../PW/src/libpw.a 
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

TLDEPS= bindir mods libs pw

all : title tldeps lib_dfpt_kernel.a

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_dfpt_kernel.a : $(DFPT_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f lib_dfpt_kernel.a

title :
	@echo "                    "
	@echo "####################"
	@echo "### DFPT_kernel ####"
	@echo "####################"
	@echo "                    "

include make.depend
