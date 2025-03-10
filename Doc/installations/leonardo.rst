.. _leonardo:

===============
CINECA-Leonardo
===============

Leonardo is a GPU-accelerated supercomputer located at `CINECA <https://www.cineca.it/en/>`_.

.. code-block:: bash

   $ ssh <username>@login.leonardo.cineca.it

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on March 10, 2025):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash
   
   module load anaconda3/2023.09-0
   module load nvhpc/23.11
   module load openmpi/4.1.6--nvhpc--23.11
   module load fftw/3.3.10--openmpi--4.1.6--nvhpc--23.11
   module load openblas/0.3.24--nvhpc--23.11
   
   export MPIF90=mpif90
   export F90=nvfortran
   export CC=nvc
   export BLAS_LIBS="-L$OPENBLAS_LIB -lopenblas"
   export LAPACK_LIBS="-L$OPENBLAS_LIB -lopenblas"
   
   ./configure --with-cuda=/leonardo/prod/opt/compilers/cuda/12.3/none --with-cuda-runtime=12.3 --with-cuda-cc=80 --with-cuda-mpi=yes
   
   make -j 8 pw
   
   cd West
   
   make conf PYT=python3 PYT_LDFLAGS="$ANACONDA3_LIB/libpython3.11.so"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh


Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on one node of Leonardo with 4 MPI ranks and 4 GPUs. The <account_name> must be replaced with an active project allocation.

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash
   #SBATCH --time=00:20:00
   #SBATCH --account=<account_name>
   #SBATCH --partition=boost_usr_prod
   #SBATCH --qos=boost_qos_dbg
   #SBATCH --nodes=1
   #SBATCH --ntasks-per-node=4
   #SBATCH --gres=gpu:4
   #SBATCH --cpus-per-task=8
   #SBATCH --threads-per-core=1

   module load anaconda3/2023.09-0
   module load nvhpc/23.11
   module load openmpi/4.1.6--nvhpc--23.11
   module load fftw/3.3.10--openmpi--4.1.6--nvhpc--23.11
   module load openblas/0.3.24--nvhpc--23.11

   export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
   export OMP_PLACES=cores; export OMP_PROC_BIND=close

   mpirun -np $SLURM_NTASKS --map-by socket:PE=$SLURM_CPUS_PER_TASK --rank-by core \
        wstat.x -i wstat.in > wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the `CINECA user guide <https://wiki.u-gov.it/confluence/display/SCAIUS/UG3.2%3A+LEONARDO+UserGuide>`_.
