.. _macos:

=====
macOS
=====

The following instructions have been tested on macOS 11.3.

Requirements:

- gcc/gfortran (e.g. GCC 9)
- MPICH
- BLAS/LAPACK/ScaLAPACK
- FFTW3
- Python3

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script:

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   export MY_LIB_PATH=/Users/myname/LIBRARIES

   export CPP='cpp-10'
   export CC='gcc-10'
   export CFLAGS='-m64 -Wall -Wextra'
   export F77='mpif77'
   export FFLAGS='-m64 -fopenmp'
   export FC='mpif90'
   export MPIF90='mpif90'
   export F90='mpif90'
   export FCFLAGS='-m64 -fopenmp -Wall -fbacktrace -fbounds-check'
   export BLAS_LIBS=${MY_LIB_PATH}/BLAS/libblas.a
   export LAPACK_LIBS=${MY_LIB_PATH}/LAPACK/liblapack.a
   export SCALAPACK_LIBS=${MY_LIB_PATH}/SCALAPACK/libscalapack.a
   export FFT_LIBS="${MY_LIB_PATH}/FFTW3/lib/libfftw3.a ${MY_LIB_PATH}/FFTW3/lib/libfftw3_omp.a"

   ./configure --enable-openmp=yes --enable-parallel=yes --enable-shared=yes --with-scalapack --with-hdf5=no

   make -j 4 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   make all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh


Running WEST
~~~~~~~~~~~~

We can run the `wstat.x` WEST executables on 2 cores using the following command:

.. code-block:: bash

   $ export OMP_NUM_THREADS=1
   $ mpirun -np 2 ./wstat.x -i wstat.in > wstat.out
