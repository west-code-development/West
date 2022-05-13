# Makefile for Base64

include ../../../make.inc

# location of needed modules
MODFLAGS = $(MOD_FLAG).

BASE64_OBJS = \
cbase64.o \
base64module.o

all : title libbase64.a

libbase64.a : $(BASE64_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f libbase64.a

title :
	@echo
	@echo "##############"
	@echo "### Base64 ###"
	@echo "##############"
	@echo
