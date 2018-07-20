#
# Parallel execution knobs
#

NP = 2  # Number of MPI processes 
NI = 1  # NUmber of images 
NT = 1  # Number of OPENMP threads 

#
# Commands needed by the scripts 
#

# Directory where the executables (*.x) are located  
BINDIR = ../../../bin

# How to run in a parallel executable with MPI
PARA_PREFIX = mpirun -np ${NP} --allow-run-as-root  

# How to download files
WGET = wget -N -q

###### DO NOT TOUCH BELOW ######

NIMAGE = ${NI}
OMP_NUM_THREADS = ${NT}

