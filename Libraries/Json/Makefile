# Makefile for Json

include ../../../make.inc

# location of needed modules
MODFLAGS = $(MOD_FLAG).

JSON_OBJS = \
json_kinds.o \
json_parameters.o \
json_string_utilities.o \
json_value_module.o \
json_file_module.o \
json_module.o

all : title libjson.a

libjson.a : $(JSON_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f libjson.a

title :
	@echo
	@echo "##############"
	@echo "### Json  ####"
	@echo "##############"
	@echo

include make.depend
