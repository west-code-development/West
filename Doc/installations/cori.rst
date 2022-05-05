.. _cori:

=================
Cori-NERSC (XC40)
=================

Cori is a Cray XC40 located at National Energy Research Scientific Computing Center (`NERSC <https://www.nersc.gov/>`_).

.. code-block:: bash

   $ ssh <username>@cori.nersc.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on Aug 1, 2021):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module unload cray-libsci
   module load cray-python/3.8.2.1

   export CRAYPE_LINK_TYPE=dynamic
   export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2019.3.199/linux/compiler/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/python/3.8.2.1/lib:$LD_LIBRARY_PATH
   export MPIF90=ftn
   export F77=ftn
   export CC=cc
   export SCALAPACK_LIBS="$MKLROOT/lib/intel64/libmkl_scalapack_lp64.so -Wl,--start-group $MKLROOT/lib/intel64/libmkl_intel_lp64.so $MKLROOT/lib/intel64/libmkl_intel_thread.so $MKLROOT/lib/intel64/libmkl_core.so $MKLROOT/lib/intel64/libmkl_blacs_intelmpi_lp64.so -Wl,--end-group"

   ./configure --enable-openmp=yes --enable-parallel=yes --enable-shared=yes --with-scalapack=intel --with-hdf5=no

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Cori (Haswell partition) with 32 MPI ranks per node. The <project_name> must be replaced with an active project allocation.

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash

   #SBATCH --job-name=WEST
   #SBATCH --time=00:20:00
   #SBATCH --account=<project_name>
   #SBATCH --constraint=haswell
   #SBATCH --qos=debug
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=32
   #SBATCH --cpus-per-task=2

   module unload cray-libsci
   module load cray-python/3.8.2.1

   export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2019.3.199/linux/compiler/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/python/3.8.2.1/lib:$LD_LIBRARY_PATH

   export OMP_NUM_THREADS=1
   export OMP_PLACE=threads
   export OMP_PROC_BIND=spread
   export MKL_NUM_THREADS=$OMP_NUM_THREADS

   NTASKS=$(($SLURM_NTASKS_PER_NODE * $SLURM_JOB_NUM_NODES))

   srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK ./wstat.x -i wstat.in &> wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the `NERSC user guide <https://docs.nersc.gov/systems/cori/>`_.
