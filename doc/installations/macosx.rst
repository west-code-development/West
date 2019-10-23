.. _macosx:

======
MacOSX
======

The following instructions have been tested on MacOSX 10.14.6.

Requirements: 

- Gcc/Gfortran 9
- MPICH 
- Blas/Lapack/Scalapack
- FFTW3
- Python3

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script: 

.. code-block:: bash 

   $ cat build_west.sh
   #!/bin/bash

   export MY_LIB_PATH=/Users/myname/LIBRARIES

   export CPP='cpp-9'
   export CC='gcc-9'
   export CFLAGS='-m64 -Wall -Wextra'
   export F77='mpif77'
   export FFLAGS='-m64 -fopenmp'
   export FC='mpif90'
   export F90='mpif90'
   export FCFLAGS='-m64 -fopenmp -Wall -fbacktrace -fbounds-check'
   export BLAS_LIBS=${MY_LIB_PATH}/BLAS/libblas.a
   export LAPACK_LIBS=${MY_LIB_PATH}/LAPACK/liblapack.a
   export SCALAPACK_LIBS=${MY_LIB_PATH}/SCALAPACK/libscalapack.a
   export FFT_LIBS="${MY_LIB_PATH}/FFTW3/lib/libfftw3.a ${MY_LIB_PATH}/FFTW3/lib/libfftw3_omp.a"
   export PYT=python3

   ./configure --with-scalapack --enable-openmp LD_LIBS="`${PYT}-config --ldflags`"
   
   make -j 4 pw
   
   cd West
   make

To use the script do: 

.. code-block:: bash 

   $ bash build_west.sh


Running WEST
~~~~~~~~~~~~

We can run the `wstat.x` WEST executables on 2 cores using the following command:

.. code-block:: bash 

   $ export OMP_NUM_THREADS=1
   $ mpirun -np 2 ./wstat.x -i wstat.in > wstat.out
