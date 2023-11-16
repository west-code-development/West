#
# Parallel execution knobs
#

export NP=8 # Number of MPI processes
export NI=1 # Number of images
export NB=1 # Number of band groups
export NK=1 # Number of pools
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

# How to download files
export WGET=wget -N -q

###### DO NOT TOUCH BELOW ######

export NIMAGE=${NI}
export NBAND=${NB}
export NPOOL=${NK}
export OMP_NUM_THREADS=${NT}

