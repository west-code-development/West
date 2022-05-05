# Makefile for IO_kernel

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../Libraries/Base64 \
           $(MOD_FLAG)../Libraries/Forpy \
           $(MOD_FLAG)../Libraries/Json \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG)../Para_kernel \
           $(MOD_FLAG)../Tools \
           $(MOD_FLAG).

IO_KERNEL_OBJS = \
cubefile.o \
function3d.o \
mod_west_io.o \
mpiio.o \
pdep_db.o \
pdep_io.o \
wfreq_db.o \
wfreq_io.o

TLDEPS = bindir mods libs pw

all : title tldeps lib_io_kernel.a

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_io_kernel.a : $(IO_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f lib_io_kernel.a

title :
	@echo "                  "
	@echo "##################"
	@echo "### IO_kernel ####"
	@echo "##################"
	@echo "                  "

include make.depend
