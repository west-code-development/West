# Makefile for Libraries

default: all

all: json_do

json_do:
	if test -d Json ; then \
	( cd Json ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

clean: json_undo 

json_undo:
	if test -d Json ; then \
	( cd Json ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

distclean: clean
