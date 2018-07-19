# Makefile
#

default: title all 

all: \
test001_all \
test002_all \
test003_all \
test004_all

clean: \
test001_clean \
test002_clean \
test003_clean \
test004_clean 

test001_all:
	if test -d test001 ; then \
	( cd test001 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi	
test001_clean:
	if test -d test001 ; then \
	( cd test001 ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

test002_all:
	if test -d test002 ; then \
	( cd test002 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi	
test002_clean:
	if test -d test002 ; then \
	( cd test002 ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

test003_all:
	if test -d test003 ; then \
	( cd test003 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi	
test003_clean:
	if test -d test003 ; then \
	( cd test003 ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

test004_all:
	if test -d test004 ; then \
	( cd test004 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi	
test004_clean:
	if test -d test004 ; then \
	( cd test004 ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

title :
	@echo "                "
	@echo "################"
	@echo "## Test-suite ##"
	@echo "################"
	@echo "                "


