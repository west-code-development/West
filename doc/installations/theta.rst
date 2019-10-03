.. _theta:

=================
Theta-ALCF (XC40)
=================

Theta is a Cray XC40 located at Argonne National Laboratory, maintained by `ALCF <https://www.alcf.anl.gov/>`_. 

.. code-block:: bash 

   $ ssh -XY <username>@theta.alcf.anl.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script: 

.. code-block:: bash 

   $ cat build_west.sh
   #!/bin/bash

   module load miniconda-3.6/conda-4.5.12

   export BLAS_LIBS="-L$MKLROOT/intel64/lib -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -Wl,--end-group"
   export SCALAPACK_LIBS="-L$MKLROOT/intel64/lib -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -Wl,--end-group"
   export FFT_LIBS=""
   export MPIF90="ftn -g -mkl -dynamic"
   export CC="cc -g -mkl -dynamic"
   export F77="ftn -g -mkl -dynamic"
   export FFLAGS="-xMIC-AVX512 -align array64byte -fp-model fast=2 -no-prec-div -assume byterecl -dynamic"
   
   export CRAYPE_LINK_TYPE=dynamic

   ./install/configure --host=x86_64-build-linux-gnu --build=x86_64-target-linux-gnu --enable-parallel --with-scalapack --enable-openmp LD_LIBS="`python3-config --ldflags`"

   make pw -j 16

   cd West
   make

To use the script do: 

.. code-block:: bash 

   $ bash build_west.sh

Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Theta with 64 MPI ranks per node. The <project_name> must be replaced with an active project allocation.

.. code-block:: bash 

   $ cat run_west.sh
   #!/bin/bash
   #COBALT -n 2 -t 10 -q debug-cache-quad -A <project_name> -O WEST
   module load miniconda-3.6/conda-4.5.12 
   aprun -n 128 -N 64 -d 1 --cc depth -e OMP_NUM_THREADS=1 -j 1 ./wstat.x -i wstat.in > wstat.out

Make the script executable: 

.. code-block:: bash 

   $ chmod +x run_west.sh

Job submission is done with the following: 

.. code-block:: bash 

   $ qsub run_west.sh

.. seealso::
   For more information, visit the ALCF user guide (`https://www.alcf.anl.gov/user-guides/xc40-system-overview <https://www.alcf.anl.gov/user-guides/xc40-system-overview/>`_).
