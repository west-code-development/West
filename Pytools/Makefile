#PYT=python3

include ../west_make.inc
include ../VERSION

all: \
title \
install 
  
install:
	cat setup.py.in | sed "s,checkit,${WEST_VERSION_NUMBER}," > setup.py
	- ${PYT} setup.py install --user --prefix= --record files.txt

clean:
	- rm -rf west.egg-info
	- rm -rf build 
	- rm -rf dist
	- rm files.txt

title :
	@echo "              "
	@echo "##############"
	@echo "### Pytools ##"
	@echo "##############"
	@echo "              "


