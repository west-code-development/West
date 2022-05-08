# Makefile for Coulomb_kernel

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG)../Tools \
           $(MOD_FLAG).

COULOMB_KERNEL_OBJS = \
class_coulomb.o \
types_coulomb.o

TLDEPS = bindir mods libs pw

all : tldeps lib_coulomb_kernel.a

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_coulomb_kernel.a : $(COULOMB_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f lib_coulomb_kernel.a

title :
	@echo
	@echo "#######################"
	@echo "### Coulomb_kernel ####"
	@echo "#######################"
	@echo

include make.depend
