.. _bebop:

==============
ANL-LCRC-Bebop
==============

Bebop is an HPC cluster maintained by the `Laboratory Computing Resource Center (LCRC) <https://www.lcrc.anl.gov/>`_ at Argonne National Laboratory.

.. code-block:: bash

   $ ssh <username>@bebop.lcrc.anl.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on December 21, 2022):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module purge
   module load git/2.31.1-6p7naeb
   module load intel-oneapi/2021.4.0.3422
   module load anaconda3/2021.05

   export MPIF90=mpiifort
   export F90=ifort
   export CC=icc
   export SCALAPACK_LIBS="-lmkl_scalapack_lp64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -Wl,--end-group"

   ./configure --with-scalapack=intel --enable-openmp
   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="-L/gpfs/fs1/home/software/anaconda3/2021.05/lib -lpython3.9"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh


Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Bebop (Broadwell partition) with 36 MPI ranks per node. The <account_name> must be replaced with an active project allocation.

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash
   #SBATCH --time=00:20:00
   #SBATCH --partition=bdwall
   #SBATCH --account=<account_name>
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=36
   #SBATCH --cpus-per-task=1

   module purge
   module load intel-oneapi/2021.4.0.3422
   module load anaconda3/2021.05

   export LD_LIBRARY_PATH=/gpfs/fs1/home/software/anaconda3/2021.05/lib:$LD_LIBRARY_PATH
   export OMP_NUM_THREADS=1

   ulimit -s unlimited

   srun -n 2 -N 72 ./wstat.x -i wstat.in > wstat.out

To run on the KNL partition, use the following flags:

.. code-block:: bash

   #SBATCH --partition=knlall
   #SBATCH --constraint knl,quad,cache
   #SBATCH --ntasks-per-node=64

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the `LCRC user guide <https://www.lcrc.anl.gov/for-users/using-lcrc/running-jobs/running-jobs-on-bebop/>`_.
