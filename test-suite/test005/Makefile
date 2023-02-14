# Makefile

include ../maketest.inc

default: title description all

all: prepare_inputs run

prepare_inputs:
	bash prepare_inputs.sh

run: pw wstat westpp

pw:
	${PARA_PREFIX_QE} ${BINDIR}/pw.x -i pw.in > pw.out 2> pw.err

wstat: pw
	${PARA_PREFIX} ${BINDIR}/wstat.x -nimage ${NIMAGE} -nband ${NBAND} -i wstat.in > wstat.out 2> wstat.err

westpp: wstat
	${PARA_PREFIX} ${BINDIR}/westpp.x -nimage ${NIMAGE} -i westpp.in > westpp.out 2> westpp.err

clean:
	- /bin/rm -rf *out *err *in *xml *upf *tab *restart *save *wfc* CRASH > /dev/null 2>&1

title:
	@echo
	@echo "###############"
	@echo "### test005 ###"
	@echo "###############"
	@echo

description:
	@echo
	@echo "SiH4 molecule gamma_only + westpp"
	@echo
