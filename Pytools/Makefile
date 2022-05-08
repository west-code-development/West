#PYT=python3

include ../west_make.inc

all: \
title \
install

install:
	- ${PYT} setup.py install --user --prefix= --record files.txt

clean:
	- rm -rf west.egg-info
	- rm -rf build
	- rm -rf dist
	- rm -f files.txt

title :
	@echo
	@echo "##############"
	@echo "### Pytools ##"
	@echo "##############"
	@echo
