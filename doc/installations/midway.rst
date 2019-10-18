.. _midway:

===============
Midway-UChicago
===============

Midway is the HPC cluster of the University of Chicago, maintained by UChicago's `RCC <https://rcc.uchicago.edu/>`_. 

.. code-block:: bash 

   $ ssh -Y <username>@midway.rcc.uchicago.edu

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script: 

.. code-block:: bash 

   $ cat build_west.sh
   #!/bin/bash
   
   module load intelmpi/5.1+intel-16.0 mkl/2017.up4 Anaconda3/5.1.0 
   
   export F77=mpiifort
   export CC=mpiicc
   export MPIF90=mpiifort
   export FC=mpiifort
   export CFLAGS="-O3 -xHost -fno-alias -ansi-alias -g -mkl -Bdynamic"
   export FFLAGS="-O3 -xHost -fno-alias -ansi-alias -g -mkl -Bdynamic"
   export BLAS_LIBS_SWITCH="external"
   export BLAS_LIBS=" -lmkl_intel_lp64  -lmkl_sequential -lmkl_core"
   export LAPACK_LIBS_SWITCH="external"
   export LAPACK_LIBS=" "
   export SCALAPACK_LIBS=" -lmkl_scalapack_lp64 -Wl,--start-group  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -Wl,--end-group"
   
   ./configure --enable-parallel --with-scalapack --enable-openmp LD_LIBS="`python3-config --ldflags`"
   make -j 6 pw
   
   cd West
   make

To use the script do: 

.. code-block:: bash 

   $ bash build_west.sh


Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of the `broadwl` partition with 28 MPI ranks per node.

.. code-block:: bash 

   $ cat run_west.sh
   #!/bin/bash
   #SBATCH --time=00:30:00
   #SBATCH --partition=broadwl
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=28
   #SBATCH --cpus-per-task=1

   module load intelmpi/5.1+intel-16.0 mkl/2017.up4 Anaconda3/5.1.0

   export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
   NTASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES))

   mpirun -np $NTASKS ./wstat.x -i wstat.in > wstat.out

Job submission is done with the following: 

.. code-block:: bash 

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the RCC user guide (`https://rcc.uchicago.edu/docs/ <https://rcc.uchicago.edu/docs/>`_).
