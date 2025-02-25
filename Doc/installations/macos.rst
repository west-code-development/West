.. _macos:

=====
macOS
=====

The following instructions have been tested on macOS 15.3 (with Apple M4).

Requirements:

- C and Fortran compilers (e.g. gcc/gfortran in `GCC <https://gcc.gnu.org/>`_)
- MPI (e.g. `OpenMPI <https://www.open-mpi.org/>`_)
- BLAS/LAPACK (e.g. `Apple Accelerate <https://developer.apple.com/documentation/accelerate/>`_)
- `FFTW3 <https://www.fftw.org/>`_
- Python3

The dependencies can be installed via `Homebrew <https://brew.sh/>`_ or compiled manually from source.

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script:

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   export MPIF90=mpif90
   export F90=gfortran
   export CC=gcc
   export BLAS_LIBS="-framework Accelerate"
   export LAPACK_LIBS="-framework Accelerate"
   export LIBDIRS="-L/PATH/TO/FFTW3/lib"

   ./configure

   # Add -I/PATH/TO/FFTW3/include to IFLAGS if needed
   make -j 4 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   make -j 4 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh


Running WEST
~~~~~~~~~~~~

We can run the `wstat.x` WEST executables on 2 cores using the following command:

.. code-block:: bash

   $ export OMP_NUM_THREADS=1
   $ mpirun -np 2 ./wstat.x -i wstat.in > wstat.out
