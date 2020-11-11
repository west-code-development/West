# Makefile for IO_kernel

include ../../make.inc

%.o : %.C
	$(CXX) -c $<


# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Libraries/Forpy \
          $(MOD_FLAG)../Libraries/Json \
          $(MOD_FLAG)../Libraries/Base64 \
          $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../Tools \
          $(MOD_FLAG)../FFT_kernel \
          $(MOD_FLAG)../Para_kernel \
          $(MOD_FLAG).
IFLAGS=


IO_KERNEL_OBJS = \
mpiio.o \
cubefile.o \
function3d.o \
pdep_io.o \
pdep_db.o \
wfreq_db.o \
mod_west_io.o \
wfreq_io.o \
wfreq_restart.o \
cpp_functions.o \
cpp_wrappers.o \
ft_interpolation.o \
qbox_interface.o  

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../Libraries/Forpy/libforpy.a ../Libraries/Json/libjson.a ../Libraries/Base64/libbase64.a

TLDEPS= bindir mods libs pw

all : title tldeps lib_io_kernel.a

tldeps:
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
