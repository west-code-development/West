# Makefile for Forpy

include ../../../make.inc

# location of needed modules
MODFLAGS = $(MOD_FLAG).

IOTK_OBJS = \
iotk_attr+CHARACTER1_0.o \
iotk_attr+COMPLEX1_0.o \
iotk_attr+COMPLEX1_3.o \
iotk_attr+COMPLEX1_6.o \
iotk_attr+COMPLEX2_0.o \
iotk_attr+COMPLEX2_3.o \
iotk_attr+COMPLEX2_6.o \
iotk_attr+COMPLEX3_0.o \
iotk_attr+COMPLEX3_3.o \
iotk_attr+COMPLEX3_6.o \
iotk_attr+COMPLEX4_0.o \
iotk_attr+COMPLEX4_3.o \
iotk_attr+COMPLEX4_6.o \
iotk_attr.o iotk_attr+INTEGER1_0.o \
iotk_attr+INTEGER1_3.o \
iotk_attr+INTEGER1_6.o \
iotk_attr+INTEGER2_0.o \
iotk_attr+INTEGER2_3.o \
iotk_attr+INTEGER2_6.o \
iotk_attr+INTEGER3_0.o \
iotk_attr+INTEGER3_3.o \
iotk_attr+INTEGER3_6.o \
iotk_attr+INTEGER4_0.o \
iotk_attr+INTEGER4_3.o \
iotk_attr+INTEGER4_6.o \
iotk_attr_interf.o \
iotk_attr+LOGICAL1_0.o \
iotk_attr+LOGICAL1_3.o \
iotk_attr+LOGICAL1_6.o \
iotk_attr+LOGICAL2_0.o \
iotk_attr+LOGICAL2_3.o \
iotk_attr+LOGICAL2_6.o \
iotk_attr+LOGICAL3_0.o \
iotk_attr+LOGICAL3_3.o \
iotk_attr+LOGICAL3_6.o \
iotk_attr+LOGICAL4_0.o \
iotk_attr+LOGICAL4_3.o \
iotk_attr+LOGICAL4_6.o \
iotk_attr+REAL1_0.o \
iotk_attr+REAL1_3.o \
iotk_attr+REAL1_6.o \
iotk_attr+REAL2_0.o \
iotk_attr+REAL2_3.o \
iotk_attr+REAL2_6.o \
iotk_attr+REAL3_0.o \
iotk_attr+REAL3_3.o \
iotk_attr+REAL3_6.o \
iotk_attr+REAL4_0.o \
iotk_attr+REAL4_3.o \
iotk_attr+REAL4_6.o \
iotk_base.o \
iotk_dat+CHARACTER1_0.o \
iotk_dat+CHARACTER1_3.o \
iotk_dat+CHARACTER1_6.o \
iotk_dat+COMPLEX1_0.o \
iotk_dat+COMPLEX1_3.o \
iotk_dat+COMPLEX1_6.o \
iotk_dat+COMPLEX2_0.o \
iotk_dat+COMPLEX2_3.o \
iotk_dat+COMPLEX2_6.o \
iotk_dat+COMPLEX3_0.o \
iotk_dat+COMPLEX3_3.o \
iotk_dat+COMPLEX3_6.o \
iotk_dat+COMPLEX4_0.o \
iotk_dat+COMPLEX4_3.o \
iotk_dat+COMPLEX4_6.o \
iotk_dat.o \
iotk_dat+INTEGER1_0.o \
iotk_dat+INTEGER1_3.o \
iotk_dat+INTEGER1_6.o \
iotk_dat+INTEGER2_0.o \
iotk_dat+INTEGER2_3.o \
iotk_dat+INTEGER2_6.o \
iotk_dat+INTEGER3_0.o \
iotk_dat+INTEGER3_3.o \
iotk_dat+INTEGER3_6.o \
iotk_dat+INTEGER4_0.o \
iotk_dat+INTEGER4_3.o \
iotk_dat+INTEGER4_6.o \
iotk_dat_interf.o \
iotk_dat+LOGICAL1_0.o \
iotk_dat+LOGICAL1_3.o \
iotk_dat+LOGICAL1_6.o \
iotk_dat+LOGICAL2_0.o \
iotk_dat+LOGICAL2_3.o \
iotk_dat+LOGICAL2_6.o \
iotk_dat+LOGICAL3_0.o \
iotk_dat+LOGICAL3_3.o \
iotk_dat+LOGICAL3_6.o \
iotk_dat+LOGICAL4_0.o \
iotk_dat+LOGICAL4_3.o \
iotk_dat+LOGICAL4_6.o \
iotk_dat+REAL1_0.o \
iotk_dat+REAL1_3.o \
iotk_dat+REAL1_6.o \
iotk_dat+REAL2_0.o \
iotk_dat+REAL2_3.o \
iotk_dat+REAL2_6.o \
iotk_dat+REAL3_0.o \
iotk_dat+REAL3_3.o \
iotk_dat+REAL3_6.o \
iotk_dat+REAL4_0.o \
iotk_dat+REAL4_3.o \
iotk_dat+REAL4_6.o \
iotk_error.o \
iotk_error_interf.o \
iotk_files.o \
iotk_files_interf.o \
iotk_fmt.o \
iotk_fmt_interf.o \
iotk_misc.o \
iotk_misc_interf.o \
iotk_module.o \
iotk_scan.o \
iotk_scan_interf.o \
iotk_stream.o \
iotk_stream_interf.o \
iotk_str.o \
iotk_str_interf.o \
iotk_tool.o \
iotk_tool_interf.o \
iotk_unit.o \
iotk_unit_interf.o \
iotk_unit_list.o \
iotk_write.o \
iotk_write_interf.o \
iotk_xtox.o \
iotk_xtox_interf.o

all : title libiotk.a

libiotk.a : $(IOTK_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f libiotk.a

title :
	@echo
	@echo "#############"
	@echo "### iotk ####"
	@echo "#############"
	@echo

iotk_attr+CHARACTER1_0.o : iotk_attr_interf.o
iotk_attr+CHARACTER1_0.o : iotk_base.o
iotk_attr+CHARACTER1_0.o : iotk_error_interf.o
iotk_attr+CHARACTER1_0.o : iotk_misc_interf.o
iotk_attr+CHARACTER1_0.o : iotk_str_interf.o
iotk_attr+COMPLEX1_0.o : iotk_attr_interf.o
iotk_attr+COMPLEX1_0.o : iotk_base.o
iotk_attr+COMPLEX1_0.o : iotk_error_interf.o
iotk_attr+COMPLEX1_0.o : iotk_fmt_interf.o
iotk_attr+COMPLEX1_0.o : iotk_misc_interf.o
iotk_attr+COMPLEX1_0.o : iotk_str_interf.o
iotk_attr+COMPLEX1_0.o : iotk_xtox_interf.o
iotk_attr+COMPLEX1_3.o : iotk_attr_interf.o
iotk_attr+COMPLEX1_3.o : iotk_base.o
iotk_attr+COMPLEX1_3.o : iotk_error_interf.o
iotk_attr+COMPLEX1_3.o : iotk_misc_interf.o
iotk_attr+COMPLEX1_3.o : iotk_str_interf.o
iotk_attr+COMPLEX1_6.o : iotk_attr_interf.o
iotk_attr+COMPLEX1_6.o : iotk_base.o
iotk_attr+COMPLEX1_6.o : iotk_error_interf.o
iotk_attr+COMPLEX1_6.o : iotk_misc_interf.o
iotk_attr+COMPLEX1_6.o : iotk_str_interf.o
iotk_attr+COMPLEX2_0.o : iotk_attr_interf.o
iotk_attr+COMPLEX2_0.o : iotk_base.o
iotk_attr+COMPLEX2_0.o : iotk_error_interf.o
iotk_attr+COMPLEX2_0.o : iotk_fmt_interf.o
iotk_attr+COMPLEX2_0.o : iotk_misc_interf.o
iotk_attr+COMPLEX2_0.o : iotk_str_interf.o
iotk_attr+COMPLEX2_0.o : iotk_xtox_interf.o
iotk_attr+COMPLEX2_3.o : iotk_attr_interf.o
iotk_attr+COMPLEX2_3.o : iotk_base.o
iotk_attr+COMPLEX2_3.o : iotk_error_interf.o
iotk_attr+COMPLEX2_3.o : iotk_misc_interf.o
iotk_attr+COMPLEX2_3.o : iotk_str_interf.o
iotk_attr+COMPLEX2_6.o : iotk_attr_interf.o
iotk_attr+COMPLEX2_6.o : iotk_base.o
iotk_attr+COMPLEX2_6.o : iotk_error_interf.o
iotk_attr+COMPLEX2_6.o : iotk_misc_interf.o
iotk_attr+COMPLEX2_6.o : iotk_str_interf.o
iotk_attr+COMPLEX3_0.o : iotk_attr_interf.o
iotk_attr+COMPLEX3_0.o : iotk_base.o
iotk_attr+COMPLEX3_0.o : iotk_error_interf.o
iotk_attr+COMPLEX3_0.o : iotk_fmt_interf.o
iotk_attr+COMPLEX3_0.o : iotk_misc_interf.o
iotk_attr+COMPLEX3_0.o : iotk_str_interf.o
iotk_attr+COMPLEX3_0.o : iotk_xtox_interf.o
iotk_attr+COMPLEX3_3.o : iotk_attr_interf.o
iotk_attr+COMPLEX3_3.o : iotk_base.o
iotk_attr+COMPLEX3_3.o : iotk_error_interf.o
iotk_attr+COMPLEX3_3.o : iotk_misc_interf.o
iotk_attr+COMPLEX3_3.o : iotk_str_interf.o
iotk_attr+COMPLEX3_6.o : iotk_attr_interf.o
iotk_attr+COMPLEX3_6.o : iotk_base.o
iotk_attr+COMPLEX3_6.o : iotk_error_interf.o
iotk_attr+COMPLEX3_6.o : iotk_misc_interf.o
iotk_attr+COMPLEX3_6.o : iotk_str_interf.o
iotk_attr+COMPLEX4_0.o : iotk_attr_interf.o
iotk_attr+COMPLEX4_0.o : iotk_base.o
iotk_attr+COMPLEX4_0.o : iotk_error_interf.o
iotk_attr+COMPLEX4_0.o : iotk_fmt_interf.o
iotk_attr+COMPLEX4_0.o : iotk_misc_interf.o
iotk_attr+COMPLEX4_0.o : iotk_str_interf.o
iotk_attr+COMPLEX4_0.o : iotk_xtox_interf.o
iotk_attr+COMPLEX4_3.o : iotk_attr_interf.o
iotk_attr+COMPLEX4_3.o : iotk_base.o
iotk_attr+COMPLEX4_3.o : iotk_error_interf.o
iotk_attr+COMPLEX4_3.o : iotk_misc_interf.o
iotk_attr+COMPLEX4_3.o : iotk_str_interf.o
iotk_attr+COMPLEX4_6.o : iotk_attr_interf.o
iotk_attr+COMPLEX4_6.o : iotk_base.o
iotk_attr+COMPLEX4_6.o : iotk_error_interf.o
iotk_attr+COMPLEX4_6.o : iotk_misc_interf.o
iotk_attr+COMPLEX4_6.o : iotk_str_interf.o
iotk_attr+INTEGER1_0.o : iotk_attr_interf.o
iotk_attr+INTEGER1_0.o : iotk_base.o
iotk_attr+INTEGER1_0.o : iotk_error_interf.o
iotk_attr+INTEGER1_0.o : iotk_fmt_interf.o
iotk_attr+INTEGER1_0.o : iotk_misc_interf.o
iotk_attr+INTEGER1_0.o : iotk_str_interf.o
iotk_attr+INTEGER1_0.o : iotk_xtox_interf.o
iotk_attr+INTEGER1_3.o : iotk_attr_interf.o
iotk_attr+INTEGER1_3.o : iotk_base.o
iotk_attr+INTEGER1_3.o : iotk_error_interf.o
iotk_attr+INTEGER1_3.o : iotk_misc_interf.o
iotk_attr+INTEGER1_3.o : iotk_str_interf.o
iotk_attr+INTEGER1_6.o : iotk_attr_interf.o
iotk_attr+INTEGER1_6.o : iotk_base.o
iotk_attr+INTEGER1_6.o : iotk_error_interf.o
iotk_attr+INTEGER1_6.o : iotk_misc_interf.o
iotk_attr+INTEGER1_6.o : iotk_str_interf.o
iotk_attr+INTEGER2_0.o : iotk_attr_interf.o
iotk_attr+INTEGER2_0.o : iotk_base.o
iotk_attr+INTEGER2_0.o : iotk_error_interf.o
iotk_attr+INTEGER2_0.o : iotk_fmt_interf.o
iotk_attr+INTEGER2_0.o : iotk_misc_interf.o
iotk_attr+INTEGER2_0.o : iotk_str_interf.o
iotk_attr+INTEGER2_0.o : iotk_xtox_interf.o
iotk_attr+INTEGER2_3.o : iotk_attr_interf.o
iotk_attr+INTEGER2_3.o : iotk_base.o
iotk_attr+INTEGER2_3.o : iotk_error_interf.o
iotk_attr+INTEGER2_3.o : iotk_misc_interf.o
iotk_attr+INTEGER2_3.o : iotk_str_interf.o
iotk_attr+INTEGER2_6.o : iotk_attr_interf.o
iotk_attr+INTEGER2_6.o : iotk_base.o
iotk_attr+INTEGER2_6.o : iotk_error_interf.o
iotk_attr+INTEGER2_6.o : iotk_misc_interf.o
iotk_attr+INTEGER2_6.o : iotk_str_interf.o
iotk_attr+INTEGER3_0.o : iotk_attr_interf.o
iotk_attr+INTEGER3_0.o : iotk_base.o
iotk_attr+INTEGER3_0.o : iotk_error_interf.o
iotk_attr+INTEGER3_0.o : iotk_fmt_interf.o
iotk_attr+INTEGER3_0.o : iotk_misc_interf.o
iotk_attr+INTEGER3_0.o : iotk_str_interf.o
iotk_attr+INTEGER3_0.o : iotk_xtox_interf.o
iotk_attr+INTEGER3_3.o : iotk_attr_interf.o
iotk_attr+INTEGER3_3.o : iotk_base.o
iotk_attr+INTEGER3_3.o : iotk_error_interf.o
iotk_attr+INTEGER3_3.o : iotk_misc_interf.o
iotk_attr+INTEGER3_3.o : iotk_str_interf.o
iotk_attr+INTEGER3_6.o : iotk_attr_interf.o
iotk_attr+INTEGER3_6.o : iotk_base.o
iotk_attr+INTEGER3_6.o : iotk_error_interf.o
iotk_attr+INTEGER3_6.o : iotk_misc_interf.o
iotk_attr+INTEGER3_6.o : iotk_str_interf.o
iotk_attr+INTEGER4_0.o : iotk_attr_interf.o
iotk_attr+INTEGER4_0.o : iotk_base.o
iotk_attr+INTEGER4_0.o : iotk_error_interf.o
iotk_attr+INTEGER4_0.o : iotk_fmt_interf.o
iotk_attr+INTEGER4_0.o : iotk_misc_interf.o
iotk_attr+INTEGER4_0.o : iotk_str_interf.o
iotk_attr+INTEGER4_0.o : iotk_xtox_interf.o
iotk_attr+INTEGER4_3.o : iotk_attr_interf.o
iotk_attr+INTEGER4_3.o : iotk_base.o
iotk_attr+INTEGER4_3.o : iotk_error_interf.o
iotk_attr+INTEGER4_3.o : iotk_misc_interf.o
iotk_attr+INTEGER4_3.o : iotk_str_interf.o
iotk_attr+INTEGER4_6.o : iotk_attr_interf.o
iotk_attr+INTEGER4_6.o : iotk_base.o
iotk_attr+INTEGER4_6.o : iotk_error_interf.o
iotk_attr+INTEGER4_6.o : iotk_misc_interf.o
iotk_attr+INTEGER4_6.o : iotk_str_interf.o
iotk_attr_interf.o : iotk_base.o
iotk_attr+LOGICAL1_0.o : iotk_attr_interf.o
iotk_attr+LOGICAL1_0.o : iotk_base.o
iotk_attr+LOGICAL1_0.o : iotk_error_interf.o
iotk_attr+LOGICAL1_0.o : iotk_fmt_interf.o
iotk_attr+LOGICAL1_0.o : iotk_misc_interf.o
iotk_attr+LOGICAL1_0.o : iotk_str_interf.o
iotk_attr+LOGICAL1_0.o : iotk_xtox_interf.o
iotk_attr+LOGICAL1_3.o : iotk_attr_interf.o
iotk_attr+LOGICAL1_3.o : iotk_base.o
iotk_attr+LOGICAL1_3.o : iotk_error_interf.o
iotk_attr+LOGICAL1_3.o : iotk_misc_interf.o
iotk_attr+LOGICAL1_3.o : iotk_str_interf.o
iotk_attr+LOGICAL1_6.o : iotk_attr_interf.o
iotk_attr+LOGICAL1_6.o : iotk_base.o
iotk_attr+LOGICAL1_6.o : iotk_error_interf.o
iotk_attr+LOGICAL1_6.o : iotk_misc_interf.o
iotk_attr+LOGICAL1_6.o : iotk_str_interf.o
iotk_attr+LOGICAL2_0.o : iotk_attr_interf.o
iotk_attr+LOGICAL2_0.o : iotk_base.o
iotk_attr+LOGICAL2_0.o : iotk_error_interf.o
iotk_attr+LOGICAL2_0.o : iotk_fmt_interf.o
iotk_attr+LOGICAL2_0.o : iotk_misc_interf.o
iotk_attr+LOGICAL2_0.o : iotk_str_interf.o
iotk_attr+LOGICAL2_0.o : iotk_xtox_interf.o
iotk_attr+LOGICAL2_3.o : iotk_attr_interf.o
iotk_attr+LOGICAL2_3.o : iotk_base.o
iotk_attr+LOGICAL2_3.o : iotk_error_interf.o
iotk_attr+LOGICAL2_3.o : iotk_misc_interf.o
iotk_attr+LOGICAL2_3.o : iotk_str_interf.o
iotk_attr+LOGICAL2_6.o : iotk_attr_interf.o
iotk_attr+LOGICAL2_6.o : iotk_base.o
iotk_attr+LOGICAL2_6.o : iotk_error_interf.o
iotk_attr+LOGICAL2_6.o : iotk_misc_interf.o
iotk_attr+LOGICAL2_6.o : iotk_str_interf.o
iotk_attr+LOGICAL3_0.o : iotk_attr_interf.o
iotk_attr+LOGICAL3_0.o : iotk_base.o
iotk_attr+LOGICAL3_0.o : iotk_error_interf.o
iotk_attr+LOGICAL3_0.o : iotk_fmt_interf.o
iotk_attr+LOGICAL3_0.o : iotk_misc_interf.o
iotk_attr+LOGICAL3_0.o : iotk_str_interf.o
iotk_attr+LOGICAL3_0.o : iotk_xtox_interf.o
iotk_attr+LOGICAL3_3.o : iotk_attr_interf.o
iotk_attr+LOGICAL3_3.o : iotk_base.o
iotk_attr+LOGICAL3_3.o : iotk_error_interf.o
iotk_attr+LOGICAL3_3.o : iotk_misc_interf.o
iotk_attr+LOGICAL3_3.o : iotk_str_interf.o
iotk_attr+LOGICAL3_6.o : iotk_attr_interf.o
iotk_attr+LOGICAL3_6.o : iotk_base.o
iotk_attr+LOGICAL3_6.o : iotk_error_interf.o
iotk_attr+LOGICAL3_6.o : iotk_misc_interf.o
iotk_attr+LOGICAL3_6.o : iotk_str_interf.o
iotk_attr+LOGICAL4_0.o : iotk_attr_interf.o
iotk_attr+LOGICAL4_0.o : iotk_base.o
iotk_attr+LOGICAL4_0.o : iotk_error_interf.o
iotk_attr+LOGICAL4_0.o : iotk_fmt_interf.o
iotk_attr+LOGICAL4_0.o : iotk_misc_interf.o
iotk_attr+LOGICAL4_0.o : iotk_str_interf.o
iotk_attr+LOGICAL4_0.o : iotk_xtox_interf.o
iotk_attr+LOGICAL4_3.o : iotk_attr_interf.o
iotk_attr+LOGICAL4_3.o : iotk_base.o
iotk_attr+LOGICAL4_3.o : iotk_error_interf.o
iotk_attr+LOGICAL4_3.o : iotk_misc_interf.o
iotk_attr+LOGICAL4_3.o : iotk_str_interf.o
iotk_attr+LOGICAL4_6.o : iotk_attr_interf.o
iotk_attr+LOGICAL4_6.o : iotk_base.o
iotk_attr+LOGICAL4_6.o : iotk_error_interf.o
iotk_attr+LOGICAL4_6.o : iotk_misc_interf.o
iotk_attr+LOGICAL4_6.o : iotk_str_interf.o
iotk_attr+REAL1_0.o : iotk_attr_interf.o
iotk_attr+REAL1_0.o : iotk_base.o
iotk_attr+REAL1_0.o : iotk_error_interf.o
iotk_attr+REAL1_0.o : iotk_fmt_interf.o
iotk_attr+REAL1_0.o : iotk_misc_interf.o
iotk_attr+REAL1_0.o : iotk_str_interf.o
iotk_attr+REAL1_0.o : iotk_xtox_interf.o
iotk_attr+REAL1_3.o : iotk_attr_interf.o
iotk_attr+REAL1_3.o : iotk_base.o
iotk_attr+REAL1_3.o : iotk_error_interf.o
iotk_attr+REAL1_3.o : iotk_misc_interf.o
iotk_attr+REAL1_3.o : iotk_str_interf.o
iotk_attr+REAL1_6.o : iotk_attr_interf.o
iotk_attr+REAL1_6.o : iotk_base.o
iotk_attr+REAL1_6.o : iotk_error_interf.o
iotk_attr+REAL1_6.o : iotk_misc_interf.o
iotk_attr+REAL1_6.o : iotk_str_interf.o
iotk_attr+REAL2_0.o : iotk_attr_interf.o
iotk_attr+REAL2_0.o : iotk_base.o
iotk_attr+REAL2_0.o : iotk_error_interf.o
iotk_attr+REAL2_0.o : iotk_fmt_interf.o
iotk_attr+REAL2_0.o : iotk_misc_interf.o
iotk_attr+REAL2_0.o : iotk_str_interf.o
iotk_attr+REAL2_0.o : iotk_xtox_interf.o
iotk_attr+REAL2_3.o : iotk_attr_interf.o
iotk_attr+REAL2_3.o : iotk_base.o
iotk_attr+REAL2_3.o : iotk_error_interf.o
iotk_attr+REAL2_3.o : iotk_misc_interf.o
iotk_attr+REAL2_3.o : iotk_str_interf.o
iotk_attr+REAL2_6.o : iotk_attr_interf.o
iotk_attr+REAL2_6.o : iotk_base.o
iotk_attr+REAL2_6.o : iotk_error_interf.o
iotk_attr+REAL2_6.o : iotk_misc_interf.o
iotk_attr+REAL2_6.o : iotk_str_interf.o
iotk_attr+REAL3_0.o : iotk_attr_interf.o
iotk_attr+REAL3_0.o : iotk_base.o
iotk_attr+REAL3_0.o : iotk_error_interf.o
iotk_attr+REAL3_0.o : iotk_fmt_interf.o
iotk_attr+REAL3_0.o : iotk_misc_interf.o
iotk_attr+REAL3_0.o : iotk_str_interf.o
iotk_attr+REAL3_0.o : iotk_xtox_interf.o
iotk_attr+REAL3_3.o : iotk_attr_interf.o
iotk_attr+REAL3_3.o : iotk_base.o
iotk_attr+REAL3_3.o : iotk_error_interf.o
iotk_attr+REAL3_3.o : iotk_misc_interf.o
iotk_attr+REAL3_3.o : iotk_str_interf.o
iotk_attr+REAL3_6.o : iotk_attr_interf.o
iotk_attr+REAL3_6.o : iotk_base.o
iotk_attr+REAL3_6.o : iotk_error_interf.o
iotk_attr+REAL3_6.o : iotk_misc_interf.o
iotk_attr+REAL3_6.o : iotk_str_interf.o
iotk_attr+REAL4_0.o : iotk_attr_interf.o
iotk_attr+REAL4_0.o : iotk_base.o
iotk_attr+REAL4_0.o : iotk_error_interf.o
iotk_attr+REAL4_0.o : iotk_fmt_interf.o
iotk_attr+REAL4_0.o : iotk_misc_interf.o
iotk_attr+REAL4_0.o : iotk_str_interf.o
iotk_attr+REAL4_0.o : iotk_xtox_interf.o
iotk_attr+REAL4_3.o : iotk_attr_interf.o
iotk_attr+REAL4_3.o : iotk_base.o
iotk_attr+REAL4_3.o : iotk_error_interf.o
iotk_attr+REAL4_3.o : iotk_misc_interf.o
iotk_attr+REAL4_3.o : iotk_str_interf.o
iotk_attr+REAL4_6.o : iotk_attr_interf.o
iotk_attr+REAL4_6.o : iotk_base.o
iotk_attr+REAL4_6.o : iotk_error_interf.o
iotk_attr+REAL4_6.o : iotk_misc_interf.o
iotk_attr+REAL4_6.o : iotk_str_interf.o
iotk_copy.o : iotk_files_interf.o
iotk_copy.o : iotk_module.o
iotk_dat+CHARACTER1_0.o : iotk_attr_interf.o
iotk_dat+CHARACTER1_0.o : iotk_base.o
iotk_dat+CHARACTER1_0.o : iotk_dat_interf.o
iotk_dat+CHARACTER1_0.o : iotk_error_interf.o
iotk_dat+CHARACTER1_0.o : iotk_fmt_interf.o
iotk_dat+CHARACTER1_0.o : iotk_misc_interf.o
iotk_dat+CHARACTER1_0.o : iotk_scan_interf.o
iotk_dat+CHARACTER1_0.o : iotk_stream_interf.o
iotk_dat+CHARACTER1_0.o : iotk_str_interf.o
iotk_dat+CHARACTER1_0.o : iotk_unit_interf.o
iotk_dat+CHARACTER1_0.o : iotk_write_interf.o
iotk_dat+CHARACTER1_3.o : iotk_attr_interf.o
iotk_dat+CHARACTER1_3.o : iotk_base.o
iotk_dat+CHARACTER1_3.o : iotk_dat_interf.o
iotk_dat+CHARACTER1_3.o : iotk_error_interf.o
iotk_dat+CHARACTER1_3.o : iotk_fmt_interf.o
iotk_dat+CHARACTER1_3.o : iotk_misc_interf.o
iotk_dat+CHARACTER1_3.o : iotk_scan_interf.o
iotk_dat+CHARACTER1_3.o : iotk_str_interf.o
iotk_dat+CHARACTER1_3.o : iotk_unit_interf.o
iotk_dat+CHARACTER1_3.o : iotk_write_interf.o
iotk_dat+CHARACTER1_6.o : iotk_attr_interf.o
iotk_dat+CHARACTER1_6.o : iotk_base.o
iotk_dat+CHARACTER1_6.o : iotk_dat_interf.o
iotk_dat+CHARACTER1_6.o : iotk_error_interf.o
iotk_dat+CHARACTER1_6.o : iotk_fmt_interf.o
iotk_dat+CHARACTER1_6.o : iotk_misc_interf.o
iotk_dat+CHARACTER1_6.o : iotk_scan_interf.o
iotk_dat+CHARACTER1_6.o : iotk_str_interf.o
iotk_dat+CHARACTER1_6.o : iotk_unit_interf.o
iotk_dat+CHARACTER1_6.o : iotk_write_interf.o
iotk_dat+COMPLEX1_0.o : iotk_attr_interf.o
iotk_dat+COMPLEX1_0.o : iotk_base.o
iotk_dat+COMPLEX1_0.o : iotk_dat_interf.o
iotk_dat+COMPLEX1_0.o : iotk_error_interf.o
iotk_dat+COMPLEX1_0.o : iotk_fmt_interf.o
iotk_dat+COMPLEX1_0.o : iotk_misc_interf.o
iotk_dat+COMPLEX1_0.o : iotk_scan_interf.o
iotk_dat+COMPLEX1_0.o : iotk_stream_interf.o
iotk_dat+COMPLEX1_0.o : iotk_str_interf.o
iotk_dat+COMPLEX1_0.o : iotk_unit_interf.o
iotk_dat+COMPLEX1_0.o : iotk_write_interf.o
iotk_dat+COMPLEX1_3.o : iotk_attr_interf.o
iotk_dat+COMPLEX1_3.o : iotk_base.o
iotk_dat+COMPLEX1_3.o : iotk_dat_interf.o
iotk_dat+COMPLEX1_3.o : iotk_error_interf.o
iotk_dat+COMPLEX1_3.o : iotk_fmt_interf.o
iotk_dat+COMPLEX1_3.o : iotk_misc_interf.o
iotk_dat+COMPLEX1_3.o : iotk_scan_interf.o
iotk_dat+COMPLEX1_3.o : iotk_str_interf.o
iotk_dat+COMPLEX1_3.o : iotk_unit_interf.o
iotk_dat+COMPLEX1_3.o : iotk_write_interf.o
iotk_dat+COMPLEX1_6.o : iotk_attr_interf.o
iotk_dat+COMPLEX1_6.o : iotk_base.o
iotk_dat+COMPLEX1_6.o : iotk_dat_interf.o
iotk_dat+COMPLEX1_6.o : iotk_error_interf.o
iotk_dat+COMPLEX1_6.o : iotk_fmt_interf.o
iotk_dat+COMPLEX1_6.o : iotk_misc_interf.o
iotk_dat+COMPLEX1_6.o : iotk_scan_interf.o
iotk_dat+COMPLEX1_6.o : iotk_str_interf.o
iotk_dat+COMPLEX1_6.o : iotk_unit_interf.o
iotk_dat+COMPLEX1_6.o : iotk_write_interf.o
iotk_dat+COMPLEX2_0.o : iotk_attr_interf.o
iotk_dat+COMPLEX2_0.o : iotk_base.o
iotk_dat+COMPLEX2_0.o : iotk_dat_interf.o
iotk_dat+COMPLEX2_0.o : iotk_error_interf.o
iotk_dat+COMPLEX2_0.o : iotk_fmt_interf.o
iotk_dat+COMPLEX2_0.o : iotk_misc_interf.o
iotk_dat+COMPLEX2_0.o : iotk_scan_interf.o
iotk_dat+COMPLEX2_0.o : iotk_stream_interf.o
iotk_dat+COMPLEX2_0.o : iotk_str_interf.o
iotk_dat+COMPLEX2_0.o : iotk_unit_interf.o
iotk_dat+COMPLEX2_0.o : iotk_write_interf.o
iotk_dat+COMPLEX2_3.o : iotk_attr_interf.o
iotk_dat+COMPLEX2_3.o : iotk_base.o
iotk_dat+COMPLEX2_3.o : iotk_dat_interf.o
iotk_dat+COMPLEX2_3.o : iotk_error_interf.o
iotk_dat+COMPLEX2_3.o : iotk_fmt_interf.o
iotk_dat+COMPLEX2_3.o : iotk_misc_interf.o
iotk_dat+COMPLEX2_3.o : iotk_scan_interf.o
iotk_dat+COMPLEX2_3.o : iotk_str_interf.o
iotk_dat+COMPLEX2_3.o : iotk_unit_interf.o
iotk_dat+COMPLEX2_3.o : iotk_write_interf.o
iotk_dat+COMPLEX2_6.o : iotk_attr_interf.o
iotk_dat+COMPLEX2_6.o : iotk_base.o
iotk_dat+COMPLEX2_6.o : iotk_dat_interf.o
iotk_dat+COMPLEX2_6.o : iotk_error_interf.o
iotk_dat+COMPLEX2_6.o : iotk_fmt_interf.o
iotk_dat+COMPLEX2_6.o : iotk_misc_interf.o
iotk_dat+COMPLEX2_6.o : iotk_scan_interf.o
iotk_dat+COMPLEX2_6.o : iotk_str_interf.o
iotk_dat+COMPLEX2_6.o : iotk_unit_interf.o
iotk_dat+COMPLEX2_6.o : iotk_write_interf.o
iotk_dat+COMPLEX3_0.o : iotk_attr_interf.o
iotk_dat+COMPLEX3_0.o : iotk_base.o
iotk_dat+COMPLEX3_0.o : iotk_dat_interf.o
iotk_dat+COMPLEX3_0.o : iotk_error_interf.o
iotk_dat+COMPLEX3_0.o : iotk_fmt_interf.o
iotk_dat+COMPLEX3_0.o : iotk_misc_interf.o
iotk_dat+COMPLEX3_0.o : iotk_scan_interf.o
iotk_dat+COMPLEX3_0.o : iotk_stream_interf.o
iotk_dat+COMPLEX3_0.o : iotk_str_interf.o
iotk_dat+COMPLEX3_0.o : iotk_unit_interf.o
iotk_dat+COMPLEX3_0.o : iotk_write_interf.o
iotk_dat+COMPLEX3_3.o : iotk_attr_interf.o
iotk_dat+COMPLEX3_3.o : iotk_base.o
iotk_dat+COMPLEX3_3.o : iotk_dat_interf.o
iotk_dat+COMPLEX3_3.o : iotk_error_interf.o
iotk_dat+COMPLEX3_3.o : iotk_fmt_interf.o
iotk_dat+COMPLEX3_3.o : iotk_misc_interf.o
iotk_dat+COMPLEX3_3.o : iotk_scan_interf.o
iotk_dat+COMPLEX3_3.o : iotk_str_interf.o
iotk_dat+COMPLEX3_3.o : iotk_unit_interf.o
iotk_dat+COMPLEX3_3.o : iotk_write_interf.o
iotk_dat+COMPLEX3_6.o : iotk_attr_interf.o
iotk_dat+COMPLEX3_6.o : iotk_base.o
iotk_dat+COMPLEX3_6.o : iotk_dat_interf.o
iotk_dat+COMPLEX3_6.o : iotk_error_interf.o
iotk_dat+COMPLEX3_6.o : iotk_fmt_interf.o
iotk_dat+COMPLEX3_6.o : iotk_misc_interf.o
iotk_dat+COMPLEX3_6.o : iotk_scan_interf.o
iotk_dat+COMPLEX3_6.o : iotk_str_interf.o
iotk_dat+COMPLEX3_6.o : iotk_unit_interf.o
iotk_dat+COMPLEX3_6.o : iotk_write_interf.o
iotk_dat+COMPLEX4_0.o : iotk_attr_interf.o
iotk_dat+COMPLEX4_0.o : iotk_base.o
iotk_dat+COMPLEX4_0.o : iotk_dat_interf.o
iotk_dat+COMPLEX4_0.o : iotk_error_interf.o
iotk_dat+COMPLEX4_0.o : iotk_fmt_interf.o
iotk_dat+COMPLEX4_0.o : iotk_misc_interf.o
iotk_dat+COMPLEX4_0.o : iotk_scan_interf.o
iotk_dat+COMPLEX4_0.o : iotk_stream_interf.o
iotk_dat+COMPLEX4_0.o : iotk_str_interf.o
iotk_dat+COMPLEX4_0.o : iotk_unit_interf.o
iotk_dat+COMPLEX4_0.o : iotk_write_interf.o
iotk_dat+COMPLEX4_3.o : iotk_attr_interf.o
iotk_dat+COMPLEX4_3.o : iotk_base.o
iotk_dat+COMPLEX4_3.o : iotk_dat_interf.o
iotk_dat+COMPLEX4_3.o : iotk_error_interf.o
iotk_dat+COMPLEX4_3.o : iotk_fmt_interf.o
iotk_dat+COMPLEX4_3.o : iotk_misc_interf.o
iotk_dat+COMPLEX4_3.o : iotk_scan_interf.o
iotk_dat+COMPLEX4_3.o : iotk_str_interf.o
iotk_dat+COMPLEX4_3.o : iotk_unit_interf.o
iotk_dat+COMPLEX4_3.o : iotk_write_interf.o
iotk_dat+COMPLEX4_6.o : iotk_attr_interf.o
iotk_dat+COMPLEX4_6.o : iotk_base.o
iotk_dat+COMPLEX4_6.o : iotk_dat_interf.o
iotk_dat+COMPLEX4_6.o : iotk_error_interf.o
iotk_dat+COMPLEX4_6.o : iotk_fmt_interf.o
iotk_dat+COMPLEX4_6.o : iotk_misc_interf.o
iotk_dat+COMPLEX4_6.o : iotk_scan_interf.o
iotk_dat+COMPLEX4_6.o : iotk_str_interf.o
iotk_dat+COMPLEX4_6.o : iotk_unit_interf.o
iotk_dat+COMPLEX4_6.o : iotk_write_interf.o
iotk_dat+INTEGER1_0.o : iotk_attr_interf.o
iotk_dat+INTEGER1_0.o : iotk_base.o
iotk_dat+INTEGER1_0.o : iotk_dat_interf.o
iotk_dat+INTEGER1_0.o : iotk_error_interf.o
iotk_dat+INTEGER1_0.o : iotk_fmt_interf.o
iotk_dat+INTEGER1_0.o : iotk_misc_interf.o
iotk_dat+INTEGER1_0.o : iotk_scan_interf.o
iotk_dat+INTEGER1_0.o : iotk_stream_interf.o
iotk_dat+INTEGER1_0.o : iotk_str_interf.o
iotk_dat+INTEGER1_0.o : iotk_unit_interf.o
iotk_dat+INTEGER1_0.o : iotk_write_interf.o
iotk_dat+INTEGER1_3.o : iotk_attr_interf.o
iotk_dat+INTEGER1_3.o : iotk_base.o
iotk_dat+INTEGER1_3.o : iotk_dat_interf.o
iotk_dat+INTEGER1_3.o : iotk_error_interf.o
iotk_dat+INTEGER1_3.o : iotk_fmt_interf.o
iotk_dat+INTEGER1_3.o : iotk_misc_interf.o
iotk_dat+INTEGER1_3.o : iotk_scan_interf.o
iotk_dat+INTEGER1_3.o : iotk_str_interf.o
iotk_dat+INTEGER1_3.o : iotk_unit_interf.o
iotk_dat+INTEGER1_3.o : iotk_write_interf.o
iotk_dat+INTEGER1_6.o : iotk_attr_interf.o
iotk_dat+INTEGER1_6.o : iotk_base.o
iotk_dat+INTEGER1_6.o : iotk_dat_interf.o
iotk_dat+INTEGER1_6.o : iotk_error_interf.o
iotk_dat+INTEGER1_6.o : iotk_fmt_interf.o
iotk_dat+INTEGER1_6.o : iotk_misc_interf.o
iotk_dat+INTEGER1_6.o : iotk_scan_interf.o
iotk_dat+INTEGER1_6.o : iotk_str_interf.o
iotk_dat+INTEGER1_6.o : iotk_unit_interf.o
iotk_dat+INTEGER1_6.o : iotk_write_interf.o
iotk_dat+INTEGER2_0.o : iotk_attr_interf.o
iotk_dat+INTEGER2_0.o : iotk_base.o
iotk_dat+INTEGER2_0.o : iotk_dat_interf.o
iotk_dat+INTEGER2_0.o : iotk_error_interf.o
iotk_dat+INTEGER2_0.o : iotk_fmt_interf.o
iotk_dat+INTEGER2_0.o : iotk_misc_interf.o
iotk_dat+INTEGER2_0.o : iotk_scan_interf.o
iotk_dat+INTEGER2_0.o : iotk_stream_interf.o
iotk_dat+INTEGER2_0.o : iotk_str_interf.o
iotk_dat+INTEGER2_0.o : iotk_unit_interf.o
iotk_dat+INTEGER2_0.o : iotk_write_interf.o
iotk_dat+INTEGER2_3.o : iotk_attr_interf.o
iotk_dat+INTEGER2_3.o : iotk_base.o
iotk_dat+INTEGER2_3.o : iotk_dat_interf.o
iotk_dat+INTEGER2_3.o : iotk_error_interf.o
iotk_dat+INTEGER2_3.o : iotk_fmt_interf.o
iotk_dat+INTEGER2_3.o : iotk_misc_interf.o
iotk_dat+INTEGER2_3.o : iotk_scan_interf.o
iotk_dat+INTEGER2_3.o : iotk_str_interf.o
iotk_dat+INTEGER2_3.o : iotk_unit_interf.o
iotk_dat+INTEGER2_3.o : iotk_write_interf.o
iotk_dat+INTEGER2_6.o : iotk_attr_interf.o
iotk_dat+INTEGER2_6.o : iotk_base.o
iotk_dat+INTEGER2_6.o : iotk_dat_interf.o
iotk_dat+INTEGER2_6.o : iotk_error_interf.o
iotk_dat+INTEGER2_6.o : iotk_fmt_interf.o
iotk_dat+INTEGER2_6.o : iotk_misc_interf.o
iotk_dat+INTEGER2_6.o : iotk_scan_interf.o
iotk_dat+INTEGER2_6.o : iotk_str_interf.o
iotk_dat+INTEGER2_6.o : iotk_unit_interf.o
iotk_dat+INTEGER2_6.o : iotk_write_interf.o
iotk_dat+INTEGER3_0.o : iotk_attr_interf.o
iotk_dat+INTEGER3_0.o : iotk_base.o
iotk_dat+INTEGER3_0.o : iotk_dat_interf.o
iotk_dat+INTEGER3_0.o : iotk_error_interf.o
iotk_dat+INTEGER3_0.o : iotk_fmt_interf.o
iotk_dat+INTEGER3_0.o : iotk_misc_interf.o
iotk_dat+INTEGER3_0.o : iotk_scan_interf.o
iotk_dat+INTEGER3_0.o : iotk_stream_interf.o
iotk_dat+INTEGER3_0.o : iotk_str_interf.o
iotk_dat+INTEGER3_0.o : iotk_unit_interf.o
iotk_dat+INTEGER3_0.o : iotk_write_interf.o
iotk_dat+INTEGER3_3.o : iotk_attr_interf.o
iotk_dat+INTEGER3_3.o : iotk_base.o
iotk_dat+INTEGER3_3.o : iotk_dat_interf.o
iotk_dat+INTEGER3_3.o : iotk_error_interf.o
iotk_dat+INTEGER3_3.o : iotk_fmt_interf.o
iotk_dat+INTEGER3_3.o : iotk_misc_interf.o
iotk_dat+INTEGER3_3.o : iotk_scan_interf.o
iotk_dat+INTEGER3_3.o : iotk_str_interf.o
iotk_dat+INTEGER3_3.o : iotk_unit_interf.o
iotk_dat+INTEGER3_3.o : iotk_write_interf.o
iotk_dat+INTEGER3_6.o : iotk_attr_interf.o
iotk_dat+INTEGER3_6.o : iotk_base.o
iotk_dat+INTEGER3_6.o : iotk_dat_interf.o
iotk_dat+INTEGER3_6.o : iotk_error_interf.o
iotk_dat+INTEGER3_6.o : iotk_fmt_interf.o
iotk_dat+INTEGER3_6.o : iotk_misc_interf.o
iotk_dat+INTEGER3_6.o : iotk_scan_interf.o
iotk_dat+INTEGER3_6.o : iotk_str_interf.o
iotk_dat+INTEGER3_6.o : iotk_unit_interf.o
iotk_dat+INTEGER3_6.o : iotk_write_interf.o
iotk_dat+INTEGER4_0.o : iotk_attr_interf.o
iotk_dat+INTEGER4_0.o : iotk_base.o
iotk_dat+INTEGER4_0.o : iotk_dat_interf.o
iotk_dat+INTEGER4_0.o : iotk_error_interf.o
iotk_dat+INTEGER4_0.o : iotk_fmt_interf.o
iotk_dat+INTEGER4_0.o : iotk_misc_interf.o
iotk_dat+INTEGER4_0.o : iotk_scan_interf.o
iotk_dat+INTEGER4_0.o : iotk_stream_interf.o
iotk_dat+INTEGER4_0.o : iotk_str_interf.o
iotk_dat+INTEGER4_0.o : iotk_unit_interf.o
iotk_dat+INTEGER4_0.o : iotk_write_interf.o
iotk_dat+INTEGER4_3.o : iotk_attr_interf.o
iotk_dat+INTEGER4_3.o : iotk_base.o
iotk_dat+INTEGER4_3.o : iotk_dat_interf.o
iotk_dat+INTEGER4_3.o : iotk_error_interf.o
iotk_dat+INTEGER4_3.o : iotk_fmt_interf.o
iotk_dat+INTEGER4_3.o : iotk_misc_interf.o
iotk_dat+INTEGER4_3.o : iotk_scan_interf.o
iotk_dat+INTEGER4_3.o : iotk_str_interf.o
iotk_dat+INTEGER4_3.o : iotk_unit_interf.o
iotk_dat+INTEGER4_3.o : iotk_write_interf.o
iotk_dat+INTEGER4_6.o : iotk_attr_interf.o
iotk_dat+INTEGER4_6.o : iotk_base.o
iotk_dat+INTEGER4_6.o : iotk_dat_interf.o
iotk_dat+INTEGER4_6.o : iotk_error_interf.o
iotk_dat+INTEGER4_6.o : iotk_fmt_interf.o
iotk_dat+INTEGER4_6.o : iotk_misc_interf.o
iotk_dat+INTEGER4_6.o : iotk_scan_interf.o
iotk_dat+INTEGER4_6.o : iotk_str_interf.o
iotk_dat+INTEGER4_6.o : iotk_unit_interf.o
iotk_dat+INTEGER4_6.o : iotk_write_interf.o
iotk_dat_interf.o : iotk_base.o
iotk_dat+LOGICAL1_0.o : iotk_attr_interf.o
iotk_dat+LOGICAL1_0.o : iotk_base.o
iotk_dat+LOGICAL1_0.o : iotk_dat_interf.o
iotk_dat+LOGICAL1_0.o : iotk_error_interf.o
iotk_dat+LOGICAL1_0.o : iotk_fmt_interf.o
iotk_dat+LOGICAL1_0.o : iotk_misc_interf.o
iotk_dat+LOGICAL1_0.o : iotk_scan_interf.o
iotk_dat+LOGICAL1_0.o : iotk_stream_interf.o
iotk_dat+LOGICAL1_0.o : iotk_str_interf.o
iotk_dat+LOGICAL1_0.o : iotk_unit_interf.o
iotk_dat+LOGICAL1_0.o : iotk_write_interf.o
iotk_dat+LOGICAL1_3.o : iotk_attr_interf.o
iotk_dat+LOGICAL1_3.o : iotk_base.o
iotk_dat+LOGICAL1_3.o : iotk_dat_interf.o
iotk_dat+LOGICAL1_3.o : iotk_error_interf.o
iotk_dat+LOGICAL1_3.o : iotk_fmt_interf.o
iotk_dat+LOGICAL1_3.o : iotk_misc_interf.o
iotk_dat+LOGICAL1_3.o : iotk_scan_interf.o
iotk_dat+LOGICAL1_3.o : iotk_str_interf.o
iotk_dat+LOGICAL1_3.o : iotk_unit_interf.o
iotk_dat+LOGICAL1_3.o : iotk_write_interf.o
iotk_dat+LOGICAL1_6.o : iotk_attr_interf.o
iotk_dat+LOGICAL1_6.o : iotk_base.o
iotk_dat+LOGICAL1_6.o : iotk_dat_interf.o
iotk_dat+LOGICAL1_6.o : iotk_error_interf.o
iotk_dat+LOGICAL1_6.o : iotk_fmt_interf.o
iotk_dat+LOGICAL1_6.o : iotk_misc_interf.o
iotk_dat+LOGICAL1_6.o : iotk_scan_interf.o
iotk_dat+LOGICAL1_6.o : iotk_str_interf.o
iotk_dat+LOGICAL1_6.o : iotk_unit_interf.o
iotk_dat+LOGICAL1_6.o : iotk_write_interf.o
iotk_dat+LOGICAL2_0.o : iotk_attr_interf.o
iotk_dat+LOGICAL2_0.o : iotk_base.o
iotk_dat+LOGICAL2_0.o : iotk_dat_interf.o
iotk_dat+LOGICAL2_0.o : iotk_error_interf.o
iotk_dat+LOGICAL2_0.o : iotk_fmt_interf.o
iotk_dat+LOGICAL2_0.o : iotk_misc_interf.o
iotk_dat+LOGICAL2_0.o : iotk_scan_interf.o
iotk_dat+LOGICAL2_0.o : iotk_stream_interf.o
iotk_dat+LOGICAL2_0.o : iotk_str_interf.o
iotk_dat+LOGICAL2_0.o : iotk_unit_interf.o
iotk_dat+LOGICAL2_0.o : iotk_write_interf.o
iotk_dat+LOGICAL2_3.o : iotk_attr_interf.o
iotk_dat+LOGICAL2_3.o : iotk_base.o
iotk_dat+LOGICAL2_3.o : iotk_dat_interf.o
iotk_dat+LOGICAL2_3.o : iotk_error_interf.o
iotk_dat+LOGICAL2_3.o : iotk_fmt_interf.o
iotk_dat+LOGICAL2_3.o : iotk_misc_interf.o
iotk_dat+LOGICAL2_3.o : iotk_scan_interf.o
iotk_dat+LOGICAL2_3.o : iotk_str_interf.o
iotk_dat+LOGICAL2_3.o : iotk_unit_interf.o
iotk_dat+LOGICAL2_3.o : iotk_write_interf.o
iotk_dat+LOGICAL2_6.o : iotk_attr_interf.o
iotk_dat+LOGICAL2_6.o : iotk_base.o
iotk_dat+LOGICAL2_6.o : iotk_dat_interf.o
iotk_dat+LOGICAL2_6.o : iotk_error_interf.o
iotk_dat+LOGICAL2_6.o : iotk_fmt_interf.o
iotk_dat+LOGICAL2_6.o : iotk_misc_interf.o
iotk_dat+LOGICAL2_6.o : iotk_scan_interf.o
iotk_dat+LOGICAL2_6.o : iotk_str_interf.o
iotk_dat+LOGICAL2_6.o : iotk_unit_interf.o
iotk_dat+LOGICAL2_6.o : iotk_write_interf.o
iotk_dat+LOGICAL3_0.o : iotk_attr_interf.o
iotk_dat+LOGICAL3_0.o : iotk_base.o
iotk_dat+LOGICAL3_0.o : iotk_dat_interf.o
iotk_dat+LOGICAL3_0.o : iotk_error_interf.o
iotk_dat+LOGICAL3_0.o : iotk_fmt_interf.o
iotk_dat+LOGICAL3_0.o : iotk_misc_interf.o
iotk_dat+LOGICAL3_0.o : iotk_scan_interf.o
iotk_dat+LOGICAL3_0.o : iotk_stream_interf.o
iotk_dat+LOGICAL3_0.o : iotk_str_interf.o
iotk_dat+LOGICAL3_0.o : iotk_unit_interf.o
iotk_dat+LOGICAL3_0.o : iotk_write_interf.o
iotk_dat+LOGICAL3_3.o : iotk_attr_interf.o
iotk_dat+LOGICAL3_3.o : iotk_base.o
iotk_dat+LOGICAL3_3.o : iotk_dat_interf.o
iotk_dat+LOGICAL3_3.o : iotk_error_interf.o
iotk_dat+LOGICAL3_3.o : iotk_fmt_interf.o
iotk_dat+LOGICAL3_3.o : iotk_misc_interf.o
iotk_dat+LOGICAL3_3.o : iotk_scan_interf.o
iotk_dat+LOGICAL3_3.o : iotk_str_interf.o
iotk_dat+LOGICAL3_3.o : iotk_unit_interf.o
iotk_dat+LOGICAL3_3.o : iotk_write_interf.o
iotk_dat+LOGICAL3_6.o : iotk_attr_interf.o
iotk_dat+LOGICAL3_6.o : iotk_base.o
iotk_dat+LOGICAL3_6.o : iotk_dat_interf.o
iotk_dat+LOGICAL3_6.o : iotk_error_interf.o
iotk_dat+LOGICAL3_6.o : iotk_fmt_interf.o
iotk_dat+LOGICAL3_6.o : iotk_misc_interf.o
iotk_dat+LOGICAL3_6.o : iotk_scan_interf.o
iotk_dat+LOGICAL3_6.o : iotk_str_interf.o
iotk_dat+LOGICAL3_6.o : iotk_unit_interf.o
iotk_dat+LOGICAL3_6.o : iotk_write_interf.o
iotk_dat+LOGICAL4_0.o : iotk_attr_interf.o
iotk_dat+LOGICAL4_0.o : iotk_base.o
iotk_dat+LOGICAL4_0.o : iotk_dat_interf.o
iotk_dat+LOGICAL4_0.o : iotk_error_interf.o
iotk_dat+LOGICAL4_0.o : iotk_fmt_interf.o
iotk_dat+LOGICAL4_0.o : iotk_misc_interf.o
iotk_dat+LOGICAL4_0.o : iotk_scan_interf.o
iotk_dat+LOGICAL4_0.o : iotk_stream_interf.o
iotk_dat+LOGICAL4_0.o : iotk_str_interf.o
iotk_dat+LOGICAL4_0.o : iotk_unit_interf.o
iotk_dat+LOGICAL4_0.o : iotk_write_interf.o
iotk_dat+LOGICAL4_3.o : iotk_attr_interf.o
iotk_dat+LOGICAL4_3.o : iotk_base.o
iotk_dat+LOGICAL4_3.o : iotk_dat_interf.o
iotk_dat+LOGICAL4_3.o : iotk_error_interf.o
iotk_dat+LOGICAL4_3.o : iotk_fmt_interf.o
iotk_dat+LOGICAL4_3.o : iotk_misc_interf.o
iotk_dat+LOGICAL4_3.o : iotk_scan_interf.o
iotk_dat+LOGICAL4_3.o : iotk_str_interf.o
iotk_dat+LOGICAL4_3.o : iotk_unit_interf.o
iotk_dat+LOGICAL4_3.o : iotk_write_interf.o
iotk_dat+LOGICAL4_6.o : iotk_attr_interf.o
iotk_dat+LOGICAL4_6.o : iotk_base.o
iotk_dat+LOGICAL4_6.o : iotk_dat_interf.o
iotk_dat+LOGICAL4_6.o : iotk_error_interf.o
iotk_dat+LOGICAL4_6.o : iotk_fmt_interf.o
iotk_dat+LOGICAL4_6.o : iotk_misc_interf.o
iotk_dat+LOGICAL4_6.o : iotk_scan_interf.o
iotk_dat+LOGICAL4_6.o : iotk_str_interf.o
iotk_dat+LOGICAL4_6.o : iotk_unit_interf.o
iotk_dat+LOGICAL4_6.o : iotk_write_interf.o
iotk_dat+REAL1_0.o : iotk_attr_interf.o
iotk_dat+REAL1_0.o : iotk_base.o
iotk_dat+REAL1_0.o : iotk_dat_interf.o
iotk_dat+REAL1_0.o : iotk_error_interf.o
iotk_dat+REAL1_0.o : iotk_fmt_interf.o
iotk_dat+REAL1_0.o : iotk_misc_interf.o
iotk_dat+REAL1_0.o : iotk_scan_interf.o
iotk_dat+REAL1_0.o : iotk_stream_interf.o
iotk_dat+REAL1_0.o : iotk_str_interf.o
iotk_dat+REAL1_0.o : iotk_unit_interf.o
iotk_dat+REAL1_0.o : iotk_write_interf.o
iotk_dat+REAL1_3.o : iotk_attr_interf.o
iotk_dat+REAL1_3.o : iotk_base.o
iotk_dat+REAL1_3.o : iotk_dat_interf.o
iotk_dat+REAL1_3.o : iotk_error_interf.o
iotk_dat+REAL1_3.o : iotk_fmt_interf.o
iotk_dat+REAL1_3.o : iotk_misc_interf.o
iotk_dat+REAL1_3.o : iotk_scan_interf.o
iotk_dat+REAL1_3.o : iotk_str_interf.o
iotk_dat+REAL1_3.o : iotk_unit_interf.o
iotk_dat+REAL1_3.o : iotk_write_interf.o
iotk_dat+REAL1_6.o : iotk_attr_interf.o
iotk_dat+REAL1_6.o : iotk_base.o
iotk_dat+REAL1_6.o : iotk_dat_interf.o
iotk_dat+REAL1_6.o : iotk_error_interf.o
iotk_dat+REAL1_6.o : iotk_fmt_interf.o
iotk_dat+REAL1_6.o : iotk_misc_interf.o
iotk_dat+REAL1_6.o : iotk_scan_interf.o
iotk_dat+REAL1_6.o : iotk_str_interf.o
iotk_dat+REAL1_6.o : iotk_unit_interf.o
iotk_dat+REAL1_6.o : iotk_write_interf.o
iotk_dat+REAL2_0.o : iotk_attr_interf.o
iotk_dat+REAL2_0.o : iotk_base.o
iotk_dat+REAL2_0.o : iotk_dat_interf.o
iotk_dat+REAL2_0.o : iotk_error_interf.o
iotk_dat+REAL2_0.o : iotk_fmt_interf.o
iotk_dat+REAL2_0.o : iotk_misc_interf.o
iotk_dat+REAL2_0.o : iotk_scan_interf.o
iotk_dat+REAL2_0.o : iotk_stream_interf.o
iotk_dat+REAL2_0.o : iotk_str_interf.o
iotk_dat+REAL2_0.o : iotk_unit_interf.o
iotk_dat+REAL2_0.o : iotk_write_interf.o
iotk_dat+REAL2_3.o : iotk_attr_interf.o
iotk_dat+REAL2_3.o : iotk_base.o
iotk_dat+REAL2_3.o : iotk_dat_interf.o
iotk_dat+REAL2_3.o : iotk_error_interf.o
iotk_dat+REAL2_3.o : iotk_fmt_interf.o
iotk_dat+REAL2_3.o : iotk_misc_interf.o
iotk_dat+REAL2_3.o : iotk_scan_interf.o
iotk_dat+REAL2_3.o : iotk_str_interf.o
iotk_dat+REAL2_3.o : iotk_unit_interf.o
iotk_dat+REAL2_3.o : iotk_write_interf.o
iotk_dat+REAL2_6.o : iotk_attr_interf.o
iotk_dat+REAL2_6.o : iotk_base.o
iotk_dat+REAL2_6.o : iotk_dat_interf.o
iotk_dat+REAL2_6.o : iotk_error_interf.o
iotk_dat+REAL2_6.o : iotk_fmt_interf.o
iotk_dat+REAL2_6.o : iotk_misc_interf.o
iotk_dat+REAL2_6.o : iotk_scan_interf.o
iotk_dat+REAL2_6.o : iotk_str_interf.o
iotk_dat+REAL2_6.o : iotk_unit_interf.o
iotk_dat+REAL2_6.o : iotk_write_interf.o
iotk_dat+REAL3_0.o : iotk_attr_interf.o
iotk_dat+REAL3_0.o : iotk_base.o
iotk_dat+REAL3_0.o : iotk_dat_interf.o
iotk_dat+REAL3_0.o : iotk_error_interf.o
iotk_dat+REAL3_0.o : iotk_fmt_interf.o
iotk_dat+REAL3_0.o : iotk_misc_interf.o
iotk_dat+REAL3_0.o : iotk_scan_interf.o
iotk_dat+REAL3_0.o : iotk_stream_interf.o
iotk_dat+REAL3_0.o : iotk_str_interf.o
iotk_dat+REAL3_0.o : iotk_unit_interf.o
iotk_dat+REAL3_0.o : iotk_write_interf.o
iotk_dat+REAL3_3.o : iotk_attr_interf.o
iotk_dat+REAL3_3.o : iotk_base.o
iotk_dat+REAL3_3.o : iotk_dat_interf.o
iotk_dat+REAL3_3.o : iotk_error_interf.o
iotk_dat+REAL3_3.o : iotk_fmt_interf.o
iotk_dat+REAL3_3.o : iotk_misc_interf.o
iotk_dat+REAL3_3.o : iotk_scan_interf.o
iotk_dat+REAL3_3.o : iotk_str_interf.o
iotk_dat+REAL3_3.o : iotk_unit_interf.o
iotk_dat+REAL3_3.o : iotk_write_interf.o
iotk_dat+REAL3_6.o : iotk_attr_interf.o
iotk_dat+REAL3_6.o : iotk_base.o
iotk_dat+REAL3_6.o : iotk_dat_interf.o
iotk_dat+REAL3_6.o : iotk_error_interf.o
iotk_dat+REAL3_6.o : iotk_fmt_interf.o
iotk_dat+REAL3_6.o : iotk_misc_interf.o
iotk_dat+REAL3_6.o : iotk_scan_interf.o
iotk_dat+REAL3_6.o : iotk_str_interf.o
iotk_dat+REAL3_6.o : iotk_unit_interf.o
iotk_dat+REAL3_6.o : iotk_write_interf.o
iotk_dat+REAL4_0.o : iotk_attr_interf.o
iotk_dat+REAL4_0.o : iotk_base.o
iotk_dat+REAL4_0.o : iotk_dat_interf.o
iotk_dat+REAL4_0.o : iotk_error_interf.o
iotk_dat+REAL4_0.o : iotk_fmt_interf.o
iotk_dat+REAL4_0.o : iotk_misc_interf.o
iotk_dat+REAL4_0.o : iotk_scan_interf.o
iotk_dat+REAL4_0.o : iotk_stream_interf.o
iotk_dat+REAL4_0.o : iotk_str_interf.o
iotk_dat+REAL4_0.o : iotk_unit_interf.o
iotk_dat+REAL4_0.o : iotk_write_interf.o
iotk_dat+REAL4_3.o : iotk_attr_interf.o
iotk_dat+REAL4_3.o : iotk_base.o
iotk_dat+REAL4_3.o : iotk_dat_interf.o
iotk_dat+REAL4_3.o : iotk_error_interf.o
iotk_dat+REAL4_3.o : iotk_fmt_interf.o
iotk_dat+REAL4_3.o : iotk_misc_interf.o
iotk_dat+REAL4_3.o : iotk_scan_interf.o
iotk_dat+REAL4_3.o : iotk_str_interf.o
iotk_dat+REAL4_3.o : iotk_unit_interf.o
iotk_dat+REAL4_3.o : iotk_write_interf.o
iotk_dat+REAL4_6.o : iotk_attr_interf.o
iotk_dat+REAL4_6.o : iotk_base.o
iotk_dat+REAL4_6.o : iotk_dat_interf.o
iotk_dat+REAL4_6.o : iotk_error_interf.o
iotk_dat+REAL4_6.o : iotk_fmt_interf.o
iotk_dat+REAL4_6.o : iotk_misc_interf.o
iotk_dat+REAL4_6.o : iotk_scan_interf.o
iotk_dat+REAL4_6.o : iotk_str_interf.o
iotk_dat+REAL4_6.o : iotk_unit_interf.o
iotk_dat+REAL4_6.o : iotk_write_interf.o
iotk_error_interf.o : iotk_base.o
iotk_error.o : iotk_base.o
iotk_error.o : iotk_error_interf.o
iotk_error.o : iotk_misc_interf.o
iotk_error.o : iotk_xtox_interf.o
iotk_files_interf.o : iotk_base.o
iotk_files.o : iotk_attr_interf.o
iotk_files.o : iotk_base.o
iotk_files.o : iotk_error_interf.o
iotk_files.o : iotk_files_interf.o
iotk_files.o : iotk_misc_interf.o
iotk_files.o : iotk_scan_interf.o
iotk_files.o : iotk_str_interf.o
iotk_files.o : iotk_unit_interf.o
iotk_files.o : iotk_write_interf.o
iotk_fmt_interf.o : iotk_base.o
iotk_fmt.o : iotk_base.o
iotk_fmt.o : iotk_fmt_interf.o
iotk_fmt.o : iotk_misc_interf.o
iotk_fmt.o : iotk_str_interf.o
iotk_fmt.o : iotk_xtox_interf.o
iotk_misc_interf.o : iotk_base.o
iotk_misc.o : iotk_attr_interf.o
iotk_misc.o : iotk_base.o
iotk_misc.o : iotk_dat_interf.o
iotk_misc.o : iotk_error_interf.o
iotk_misc.o : iotk_misc_interf.o
iotk_misc.o : iotk_scan_interf.o
iotk_misc.o : iotk_str_interf.o
iotk_misc.o : iotk_unit_interf.o
iotk_misc.o : iotk_write_interf.o
iotk_misc.o : iotk_xtox_interf.o
iotk_module.o : iotk_attr_interf.o
iotk_module.o : iotk_base.o
iotk_module.o : iotk_dat_interf.o
iotk_module.o : iotk_error_interf.o
iotk_module.o : iotk_files_interf.o
iotk_module.o : iotk_fmt_interf.o
iotk_module.o : iotk_misc_interf.o
iotk_module.o : iotk_scan_interf.o
iotk_module.o : iotk_tool_interf.o
iotk_module.o : iotk_unit_interf.o
iotk_module.o : iotk_write_interf.o
iotk_module.o : iotk_xtox_interf.o
iotk.o : iotk_base.o
iotk.o : iotk_error_interf.o
iotk.o : iotk_module.o
iotk_print_kinds.o : iotk_module.o
iotk_scan_interf.o : iotk_base.o
iotk_scan.o : iotk_attr_interf.o
iotk_scan.o : iotk_base.o
iotk_scan.o : iotk_error_interf.o
iotk_scan.o : iotk_files_interf.o
iotk_scan.o : iotk_misc_interf.o
iotk_scan.o : iotk_scan_interf.o
iotk_scan.o : iotk_stream_interf.o
iotk_scan.o : iotk_str_interf.o
iotk_scan.o : iotk_unit_interf.o
iotk_stream_interf.o : iotk_base.o
iotk_stream.o : iotk_base.o
iotk_stream.o : iotk_error_interf.o
iotk_stream.o : iotk_stream_interf.o
iotk_str_interf.o : iotk_base.o
iotk_str.o : iotk_base.o
iotk_str.o : iotk_error_interf.o
iotk_str.o : iotk_misc_interf.o
iotk_str.o : iotk_str_interf.o
iotk_tool.o : iotk_base.o
iotk_tool.o : iotk_error_interf.o
iotk_tool.o : iotk_files_interf.o
iotk_tool.o : iotk_misc_interf.o
iotk_tool.o : iotk_str_interf.o
iotk_tool.o : iotk_tool_interf.o
iotk_tool.o : iotk_xtox_interf.o
iotk_unit_interf.o : iotk_base.o
iotk_unit_list.o : iotk_base.o
iotk_unit.o : iotk_base.o
iotk_unit.o : iotk_error_interf.o
iotk_unit.o : iotk_misc_interf.o
iotk_unit.o : iotk_str_interf.o
iotk_unit.o : iotk_unit_interf.o
iotk_unit.o : iotk_unit_list.o
iotk_write_interf.o : iotk_base.o
iotk_write.o : iotk_base.o
iotk_write.o : iotk_error_interf.o
iotk_write.o : iotk_files_interf.o
iotk_write.o : iotk_misc_interf.o
iotk_write.o : iotk_str_interf.o
iotk_write.o : iotk_unit_interf.o
iotk_write.o : iotk_write_interf.o
iotk_xtox_interf.o : iotk_base.o
iotk_xtox.o : iotk_base.o
iotk_xtox.o : iotk_misc_interf.o
