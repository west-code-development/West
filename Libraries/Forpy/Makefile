# Makefile for Forpy

include ../../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG).
IFLAGS=

FORPY_OBJS = \
forpy_mod.o

#TLDEPS= bindir mods libs pw

#all : title tldeps libjson.a
all : title libforpy.a

#tldeps:
#	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libforpy.a : $(FORPY_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f libforpy.a

title :
	@echo "              "
	@echo "##############"
	@echo "### Forpy ####"
	@echo "##############"
	@echo "              "

