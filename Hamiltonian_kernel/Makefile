# Makefile for Hamiltonian_kernel

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG).

HAMILTONIAN_KERNEL_OBJS = \
apply_alpha_pc_to_m_wfcs.o \
apply_alpha_pv_to_m_wfcs.o \
commut_Hx_psi.o \
glbrak.o \
k_psi.o

TLDEPS = bindir mods libs pw

all : title tldeps lib_hamiltonian_kernel.a

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_hamiltonian_kernel.a : $(HAMILTONIAN_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f lib_hamiltonian_kernel.a

title :
	@echo "                           "
	@echo "###########################"
	@echo "### Hamiltonian_kernel ####"
	@echo "###########################"
	@echo "                           "

include make.depend
