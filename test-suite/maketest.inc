#
# Parallel execution knobs
#

export NP=2  # Number of MPI processes 
export NI=1  # NUmber of images 
export NT=1  # Number of OPENMP threads 

#
# Commands needed by the scripts 
#

# Directory where the executables (*.x) are located  
export BINDIR = ../../../bin

# How to run in a parallel executable with MPI
export PARA_PREFIX = mpirun -np ${NP} 

# How to download files
export WGET=wget -N -q

###### DO NOT TOUCH BELOW ######

export NIMAGE = ${NI}
export OMP_NUM_THREADS = ${NT}

