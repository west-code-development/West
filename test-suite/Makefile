# Makefile

###############
# Directories # Note that the directories to check are listed in main_test.py
###############

testdirs:= test001 test002 test003 test004 test005 test006 test007 test008 test009 test010 test011 test012 test013 test014 test015 test016 test017 test018 test019 test020 test021

default: all

all:
	@echo
	@echo "#################"
	@echo "## RUNNING...  ##"
	@echo "#################"
	@echo
	for d in $(testdirs); do cd $$d; $(MAKE); cd ..; done
	@echo
	@echo "#################"
	@echo "## CHECKING... ##"
	@echo "#################"
	@echo
	python3 -m pytest -v --capture=tee-sys --durations=0

clean:
	for d in $(testdirs); do cd $$d; $(MAKE) clean; cd ..; done
	- /bin/rm -rf __pycache__ .pytest_cache > /dev/null 2>&1
