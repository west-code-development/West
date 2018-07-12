
NP = 2
NI = 1
NT = 1

NIMAGE = ${NI}
OMP_NUM_THREADS = ${NT}
PARA_PREFIX = mpirun -np ${NP} --allow-run-as-root

WGET = wget
BINDIR = ../../../bin
