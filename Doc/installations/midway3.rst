.. _midway3:

================
Midway3-UChicago
================

Midway3 is the HPC cluster of the University of Chicago, maintained by UChicago's `RCC <https://rcc.uchicago.edu/>`_.

.. code-block:: bash

   $ ssh <username>@midway3.rcc.uchicago.edu

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script:

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module load intel/19.1.1
   module load intelmpi/2019.up7+intel-19.1.1
   module load mkl/2020.up1
   module load python/anaconda-2020.11

   export MPIF90=mpiifort
   export F90=ifort
   export F77=ifort
   export CC=icc
   export SCALAPACK_LIBS="-lmkl_scalapack_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -Wl,--end-group"

   ./configure --enable-parallel --with-scalapack --enable-openmp
   make -j 8 pw

   cd West
   make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   sed -i 's/-L.*config-3.8-x86_64-linux-gnu //' west_make.inc
   make all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh


Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Midway3 with 32 MPI ranks per node. The <project_name> and <account_name> must be replaced with an active project allocation.

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash
   #SBATCH --time=00:20:00
   #SBATCH --partition=<partition_name>
   #SBATCH --account=<account_name>
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=48
   #SBATCH --cpus-per-task=1

   module load intel/19.1.1
   module load intelmpi/2019.up7+intel-19.1.1
   module load mkl/2020.up1
   module load python/anaconda-2020.11

   export I_MPI_PMI_LIBRARY=/software/slurm-current-$DISTARCH/lib/libpmi.so
   export LD_LIBRARY_PATH=/software/python-anaconda-2020.11-el8-x86_64/lib:$LD_LIBRARY_PATH
   export OMP_NUM_THREADS=1

   srun -n 96 -N 2 ./wstat.x -i wstat.in > wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the `RCC user guide <https://rcc.uchicago.edu/docs/>`_.
