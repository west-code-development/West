# Makefile for Modules

include ../VERSION
include ../west_make.inc
include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Libraries/Forpy \
          $(MOD_FLAG).
IFLAGS=


MODULES_OBJS = \
west_version.o \
westcom.o \

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../Libraries/Forpy/libforpy.a

TLDEPS= bindir mods libs pw

all: title tldeps update libmodules.a

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

update:	
	./update_west_version ${WESTDIR} ${WEST_VERSION_NUMBER}

libmodules.a : $(MODULES_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L 
	- /bin/rm -f libmodules.a
	- /bin/rm -f west_version.f90

title :
	@echo "                "
	@echo "################"
	@echo "### Modules ####"
	@echo "################"
	@echo "                "

include make.depend
