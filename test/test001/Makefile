# Makefile

include ../maketest.inc

default: title all

all: prepare_input run 

prepare_input: 
	-wget -N -q http://www.quantum-simulation.org/potentials/sg15_oncv/upf/H_ONCV_PBE-1.0.upf
	-wget -N -q http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Si_ONCV_PBE-1.1.upf 
	-wget -N -q http://www.west-code.org/doc/training/silane/pw.in
	-wget -N -q http://www.west-code.org/doc/training/silane/wstat.in
	-wget -N -q http://www.west-code.org/doc/training/silane/wfreq.in

run: 
	- ${PARA_PREFIX} ${BINDIR}/pw.x -i pw.in > pw.out
	- ${PARA_PREFIX} ${BINDIR}/wstat.x -nimage ${NIMAGE} -i wstat.in > wstat.out
	- ${PARA_PREFIX} ${BINDIR}/wfreq.x -nimage ${NIMAGE} -i wfreq.in > wfreq.out

clean:
	- /bin/rm -rf *out *in *upf *save >& /dev/null

title:
	@echo "               "
	@echo "###############"
	@echo "### test001 ###"
	@echo "###############"
	@echo "               "
