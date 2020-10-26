# Makefile

include ../maketest.inc

default: title description all

all: prepare_inputs run 

prepare_inputs: 
	bash prepare_inputs.sh

run: pw wstat wfreq

pw:
	${PARA_PREFIX} ${BINDIR}/pw.x -i pw.in > pw.out
	tail -5 pw.out
wstat:
	${PARA_PREFIX} ${BINDIR}/wstat.x -nimage ${NIMAGE} -i wstat.in > wstat.out
	tail -5 wstat.out
wfreq:
	${PARA_PREFIX} ${BINDIR}/wfreq.x -nimage ${NIMAGE} -i wfreq.in > wfreq.out
	tail -5 wfreq.out

clean:
	- /bin/rm -rf *out *in *upf *tab *restart *save >& /dev/null

title:
	@echo "               "
	@echo "###############"
	@echo "### test001 ###"
	@echo "###############"
	@echo "               "

description:
	@echo " "
	@echo "SiH4 molecule gamma_only"
	@echo " "
