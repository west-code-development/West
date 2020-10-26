# Makefile

include ../maketest.inc

default: title description all

all: prepare_inputs run 

prepare_inputs: 
	bash prepare_inputs.sh

run: pw wstat1 wstat2 wfreq 

pw:  
	${PARA_PREFIX} ${BINDIR}/pw.x -i pw.in > pw.out
	tail -5 pw.out
wstat1:
	${PARA_PREFIX} ${BINDIR}/wstat.x -nimage ${NIMAGE} -i wstat1.in > wstat1.out
	tail -5 wstat1.out
wstat2:
	${PARA_PREFIX} ${BINDIR}/wstat.x -nimage ${NIMAGE} -i wstat2.in > wstat2.out
	tail -5 wstat2.out
wfreq:
	${PARA_PREFIX} ${BINDIR}/wfreq.x -nimage ${NIMAGE} -i wfreq.in > wfreq.out
	tail -5 wfreq.out

clean:
	- /bin/rm -rf *out *in *upf *tab *restart *save >& /dev/null

title:
	@echo "               "
	@echo "###############"
	@echo "### test002 ###"
	@echo "###############"
	@echo "               "

description:
	@echo " "
	@echo "CH4 molecule gamma_only test n_read_from_file"
	@echo " "
