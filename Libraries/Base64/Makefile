# Makefile for Base64

include ../../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG).
IFLAGS=

BASE64_OBJS = \
Base64.o \
Base64Transcoder.o \
Base64Module.o

.SUFFIXES : .o .c .f .f90 .C 
.C.o:
	$(CC) $(CFLAGS) -c $< -lsdtc++ 

#all : title tldeps libbase64.a
all : title libbase64.a

#tldeps:
#	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libbase64.a : $(BASE64_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f libbase64.a

title :
	@echo "              "
	@echo "##############"
	@echo "### Base64 ###"
	@echo "##############"
	@echo "              "

