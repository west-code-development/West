# Makefile for Json-test

include ../../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../Json \
          $(MOD_FLAG).
IFLAGS=


JSONTEST_OBJS = \
test.o 

JSON_OBJS = ../Json/libjson.a

TLDEPS= bindir mods libs pw

#all : title tldeps wstat.x
all : title tldeps wstat.x

wstat.x : $(JSON_OBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o test.x \
		$(JSONTEST_OBJS) $(JSON_OBJS) $(LIBOBJS) $(LIBS)
#	- ( cd ../../bin ; ln -fs ../West/Wstat/wstat.x . )
#tldeps:
#	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
#	- /bin/rm -f ../../bin/wstat.x

title :
	@echo "              "
	@echo "##############"
	@echo "# JSON test ##"
	@echo "##############"
	@echo "              "

include make.depend
