# Makefile for Para_kernel

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../Tools \
           $(MOD_FLAG).

PARA_KERNEL_OBJS = \
class_idistribute.o \
distribution_center.o \
parallel_distributed_diago.o \
west_mp.o

TLDEPS = bindir mods pwlibs pw

all : title tldeps lib_para_kernel.a

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_para_kernel.a : $(PARA_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f lib_para_kernel.a

title :
	@echo
	@echo "####################"
	@echo "### Para_kernel ####"
	@echo "####################"
	@echo

include make.depend
