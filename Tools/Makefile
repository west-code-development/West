# Makefile for Tools

include ../../make.inc

# location of needed modules
MODFLAGS = $(BASEMOD_FLAGS) \
           $(MOD_FLAG)../../PW/src \
           $(MOD_FLAG)../Libraries/Base64 \
           $(MOD_FLAG)../Libraries/Forpy \
           $(MOD_FLAG)../Libraries/Json \
           $(MOD_FLAG)../Modules \
           $(MOD_FLAG).

TOOLS_OBJS = \
bar.o \
class_bz_grid.o \
clean_scratchfiles.o \
conversions.o \
destroy_pw_arrays.o \
do_setup.o \
exx_go.o \
exx_ungo.o \
fetch_input.o \
get_alpha_pv.o \
gradcorr_x.o \
heapsort.o \
human_readable_time.o \
init_pw_arrays.o \
io_push.o \
linear_algebra_kernel.o \
my_mkdir.o \
parse_command_arguments.o \
report_dynamical_memory.o \
set_dirs.o \
set_eprec.o \
set_nbndocc.o \
set_npwq.o \
types_bz_grid.o \
v_x.o \
west_environment.o \
west_print_clocks.o

TLDEPS = bindir mods libs pw

all : title tldeps libtools.a

tldeps :
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libtools.a : $(TOOLS_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f libtools.a

title :
	@echo
	@echo "##############"
	@echo "### Tools ####"
	@echo "##############"
	@echo

include make.depend
