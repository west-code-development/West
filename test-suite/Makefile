# Makefile

testdirs:= test001 test002 test003 test004 test005 test006 test007

default: title all

all:
	for d in $(testdirs); do cd $$d; $(MAKE); cd ..; done
	python3 check.py

clean:
	for d in $(testdirs); do cd $$d; $(MAKE) clean; cd ..; done

title :
	@echo
	@echo "################"
	@echo "## Test-suite ##"
	@echo "################"
	@echo

