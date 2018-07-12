# Makefile
#

default: all 

all: \
test001_all 

run: \
test001_run 

clean: \
test001_clean 

test001_all:
	if test -d test001 ; then \
	( cd test001 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi	

test001_run:
	if test -d test001 ; then \
	( cd test001 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) run; \
	else $(MAKE) $(MFLAGS) run ; fi ) ; fi	

test001_clean:
	if test -d test001 ; then \
	( cd test001 ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

title :
	@echo "             "
	@echo "#############"
	@echo "## Testing ##"
	@echo "#############"
	@echo "             "


