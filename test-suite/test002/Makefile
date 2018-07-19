# Makefile

include ../maketest.inc

default: title all

all: prepare_inputs run 

prepare_inputs: 
	-bash prepare_inputs.sh

run: 
	- ${PARA_PREFIX} ${BINDIR}/pw.x -i pw.in > pw.out
	- tail -5 pw.out
	- ${PARA_PREFIX} ${BINDIR}/wstat.x -nimage ${NIMAGE} -i wstat.in > wstat.out
	- tail -5 wstat.out
	- ${PARA_PREFIX} ${BINDIR}/wfreq.x -nimage ${NIMAGE} -i wfreq.in > wfreq.out
	- tail -5 wfreq.out

clean:
	- /bin/rm -rf *out *in *upf *save >& /dev/null

title:
	@echo "               "
	@echo "###############"
	@echo "### test002 ###"
	@echo "###############"
	@echo "               "
