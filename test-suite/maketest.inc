#
# Parallel execution knobs
#

export NP=2 # Number of MPI processes
export NI=1 # Number of images
export NB=1 # Number of band groups
export NT=1 # Number of OpenMP threads

#
# Commands needed by the scripts
#

# Directory where the executables (*.x) are located
export BINDIR=../../../bin

# How to run West parallel executables with MPI
export PARA_PREFIX=mpirun -np ${NP}

# How to run QE parallel executables with MPI
export PARA_PREFIX_QE=mpirun -np 2

# Whether to test GPU version
export USEGPU=no

# How to download files
export WGET=wget -N -q

###### DO NOT TOUCH BELOW ######

export NIMAGE=${NI}
export NBAND=${NB}
export OMP_NUM_THREADS=${NT}

