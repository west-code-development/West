# Makefile

include ../maketest.inc

default: title description all

all: prepare_inputs run

prepare_inputs:
	bash prepare_inputs.sh

run: pw westpp

pw:
	${PARA_PREFIX_QE} ${BINDIR}/pw.x -i pw.in > pw.out 2> pw.err

westpp: pw
	${PARA_PREFIX} ${BINDIR}/westpp.x -nimage ${NIMAGE} -i westpp.in > westpp.out 2> westpp.err

clean:
	- /bin/rm -rf *out *err *in *xml *upf *tab *restart *save *wfc* CRASH > /dev/null 2>&1

title:
	@echo
	@echo "###############"
	@echo "### test008 ###"
	@echo "###############"
	@echo

description:
	@echo
	@echo "MgO gamma_only test localization"
	@echo
