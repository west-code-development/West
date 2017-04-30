# Makefile for Json

include ../../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG).
IFLAGS=

JSON_OBJS = \
json_kinds.o \
json_parameters.o \
json_string_utilities.o \
json_value_module.o \
json_file_module.o \
json_module.o  

#TLDEPS= bindir mods libs pw

#all : title tldeps libjson.a
all : title libjson.a copy

#tldeps:
#	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libjson.a : $(JSON_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

copy : 
	- cp libjson.a ../

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f libjson.a

title :
	@echo "              "
	@echo "##############"
	@echo "### Json  ####"
	@echo "##############"
	@echo "              "

include make.depend
