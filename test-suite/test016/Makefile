# Makefile

include ../maketest.inc

default: title description all

all: prepare_inputs run

prepare_inputs:
	bash prepare_inputs.sh

run: pw wstat wbse_init wbse

pw:
	${PARA_PREFIX_QE} ${BINDIR}/pw.x -i pw.in > pw.out 2> pw.err

wstat: pw
	${PARA_PREFIX} ${BINDIR}/wstat.x -nimage ${NIMAGE} -nband ${NBAND} -i wstat.in > wstat.out 2> wstat.err

wbse_init: wstat
	${PARA_PREFIX} ${BINDIR}/wbse_init.x -nimage ${NIMAGE} -nband ${NBAND} -i wbse_init.in > wbse_init.out 2> wbse_init.err

wbse: wbse_init
	${PARA_PREFIX} ${BINDIR}/wbse.x -nimage ${NIMAGE} -i wbse.in > wbse.out 2> wbse.err

clean:
	- /bin/rm -rf *out *err *in *xml *upf *tab *restart *save *wfc* CRASH > /dev/null 2>&1

title:
	@echo
	@echo "###############"
	@echo "### test016 ###"
	@echo "###############"
	@echo

description:
	@echo
	@echo "SiH4 molecule test BSE Lanczos"
	@echo
