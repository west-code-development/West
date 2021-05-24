# Makefile for Json-test

include ../../../make.inc

# location of needed modules
MODFLAGS = $(MOD_FLAG)../Json \
          $(MOD_FLAG).

JSONTEST_OBJS = \
test.o

JSON_OBJS = ../Json/libjson.a

all : title test.x

test.x : $(JSONTEST_OBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o test.x \
		$(JSONTEST_OBJS) $(JSON_OBJS) $(LIBOBJS) $(LIBS)
	- ./test.x
	- python check.py

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L *.json

title :
	@echo "              "
	@echo "##############"
	@echo "# JSON test ##"
	@echo "##############"
	@echo "              "

include make.depend
