# Makefile
#
testdirs:= test001 test002 test003 test003

default: title all
all:
	for d in $(testdirs); do (cd $$d; $(MAKE)); done
clean:
	for d in $(testdirs); do (cd $$d; $(MAKE) clean); done

title :
	@echo "                "
	@echo "################"
	@echo "## Test-suite ##"
	@echo "################"
	@echo "                "

