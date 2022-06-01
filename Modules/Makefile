# Makefile for Modules

include ../west_make.inc
include ../../make.inc

# location of needed modules
MODFLAGS = $(MOD_FLAG)../../Modules \
           $(MOD_FLAG)../../FFTXlib \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../../UtilXlib \
           $(MOD_FLAG)../../upflib \
           $(MOD_FLAG)../Libraries/Forpy \
           $(MOD_FLAG).

MODULES_OBJS = \
west_version.o \
westcom.o

TLDEPS = bindir mods pwlibs pw

all : title tldeps update libmodules.a

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

update :
	./update_west_version ${WESTDIR} `${PYT} ../Pytools/read_json.py ../VERSION.json version`

libmodules.a : $(MODULES_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f libmodules.a
	- /bin/rm -f west_version.f90

title :
	@echo
	@echo "################"
	@echo "### Modules ####"
	@echo "################"
	@echo

include make.depend
