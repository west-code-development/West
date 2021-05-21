# Makefile for Libraries

default: all

all: title forpy_do json_do base64_do

test: json_test_do

forpy_do:
	if test -d Forpy ; then \
	( cd Forpy ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

json_do:
	if test -d Json ; then \
	( cd Json ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

base64_do:
	if test -d Base64 ; then \
	( cd Base64 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

json_test_do:
	if test -d Json-test ; then \
	( cd Json-test ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

clean: forpy_undo json_undo base64_undo json_test_undo

forpy_undo:
	if test -d Forpy ; then \
	( cd Forpy ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

json_undo:
	if test -d Json ; then \
	( cd Json ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

base64_undo:
	if test -d Base64 ; then \
	( cd Base64 ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

json_test_undo:
	if test -d Json-test ; then \
	( cd Json-test ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

distclean: clean

title :
	@echo "              "
	@echo "##############"
	@echo "## Libraries #"
	@echo "##############"
	@echo "              "
