# Makefile for Forpy

include ../../../make.inc

# location of needed modules
MODFLAGS = $(MOD_FLAG).

FORPY_OBJS = \
forpy_mod.o

all : title libforpy.a

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
