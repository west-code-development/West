# Makefile
#
testdirs:= test001 test002 test003 test004

default: title all
all:
	for d in $(testdirs); do cd $$d; $(MAKE); [[ $$? != 0 ]] && exit -1; echo 'done'; cd ..; done
clean:
	for d in $(testdirs); do cd $$d; $(MAKE) clean; cd ..; done

title :
	@echo "                "
	@echo "################"
	@echo "## Test-suite ##"
	@echo "################"
	@echo "                "

