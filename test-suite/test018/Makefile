# Makefile

include ../maketest.inc

default: title description all

all: prepare_inputs run

prepare_inputs:
	bash prepare_inputs.sh

run: pw wbse_init wbse

pw:
	${PARA_PREFIX_QE} ${BINDIR}/pw.x -i pw.in > pw.out 2> pw.err

wbse_init: pw
	${PARA_PREFIX} ${BINDIR}/wbse_init.x -nimage ${NIMAGE} -nband ${NBAND} -i wbse_init.in > wbse_init.out 2> wbse_init.err

wbse: wbse_init
	${PARA_PREFIX} ${BINDIR}/wbse.x -nband ${NBAND} -i wbse.in > wbse.out 2> wbse.err

clean:
	- /bin/rm -rf *out *err *in *xml *upf *tab *restart *save *wfc* CRASH > /dev/null 2>&1

title:
	@echo
	@echo "###############"
	@echo "### test018 ###"
	@echo "###############"
	@echo

description:
	@echo
	@echo "SiH4 molecule test hybrid TDDFT Lanczos"
	@echo
