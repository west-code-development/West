# Makefile for Tools

include ../../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Libraries/Forpy \
          $(MOD_FLAG)../Libraries/Json \
          $(MOD_FLAG)../Libraries/Base64 \
          $(MOD_FLAG)../Modules \
          $(MOD_FLAG).
IFLAGS=


TOOLS_OBJS = \
io_push.o \
fetch_input.o \
set_dirs.o \
west_print_clocks.o \
exx_go.o \
exx_ungo.o \
report_dynamical_memory.o \
parse_command_arguments.o \
class_bz_grid.o \
types_bz_grid.o \
do_setup.o \
linear_algebra_kernel.o \
human_readable_time.o \
destroy_pw_arrays.o \
init_pw_arrays.o \
clean_scratchfiles.o \
set_nbndocc.o \
get_alpha_pv.o \
set_iks_l2g.o \
v_x.o \
gradcorr_x.o \
west_environment.o \
set_npwq.o \
set_eprec.o \
read_pwout.o \
bar.o  

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../Libraries/Forpy/libforpy.a ../Libraries/Json/libjson.a ../Libraries/Base64/libbase64.a

TLDEPS= bindir mods libs pw

all : title tldeps libtools.a

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libtools.a : $(TOOLS_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f libtools.a

title :
	@echo "              "
	@echo "##############"
	@echo "### Tools ####"
	@echo "##############"
	@echo "              "

include make.depend
