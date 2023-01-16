.. _perlmutter:

================
NERSC-Perlmutter
================

Perlmutter (Phase I) is a GPU-accelerated supercomputer located at National Energy Research Scientific Computing Center (`NERSC <https://www.nersc.gov/>`_).

.. code-block:: bash

   $ ssh <username>@saul-p1.nersc.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on January 13, 2023):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module load PrgEnv-nvidia
   module load nvidia/22.9
   module load cudatoolkit/11.7
   module load craype-accel-nvidia80
   module load cray-python/3.9.13.1
   module unload darshan

   ./configure --with-cuda=$CUDA_HOME --with-cuda-runtime=11.7 --with-cuda-cc=80

   # Manually edit make.inc:

   # DFLAGS = -D__PGI -D__CUDA -D__FFTW -D__MPI
   # MPIF90 = ftn
   # F90 = ftn
   # CC = cc
   # LD = ftn
   # BLAS_LIBS = # leave blank
   # LAPACK_LIBS = # leave blank

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Perlmutter with 4 MPI ranks and 4 GPUs per node. The <project_name> must be replaced with an active project allocation.

**Important**: The following environment variable is needed to work around a bug in ROMIO, Cray MPICH.

.. code-block:: bash

   export ROMIO_FSTYPE_FORCE="ufs:"

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash

   #SBATCH --job-name=WEST
   #SBATCH --time=00:20:00
   #SBATCH --account=<project_name>
   #SBATCH --constraint=gpu
   #SBATCH --qos=debug
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=4
   #SBATCH --gpus-per-node=4
   #SBATCH --cpus-per-task=32

   module load PrgEnv-nvidia
   module load nvidia/22.9
   module load cudatoolkit/11.7
   module load craype-accel-nvidia80
   module load cray-python/3.9.13.1
   module unload darshan

   export OMP_NUM_THREADS=1
   export SLURM_CPU_BIND=cores
   export MPICH_MPIIO_HINTS=*:romio_cb_write=enable:romio_ds_write=disable
   export ROMIO_FSTYPE_FORCE="ufs:"

   srun -N 2 -n 8 -c 32 -G 8 ./wstat.x -i wstat.in &> wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the `NERSC user guide <https://docs.nersc.gov/systems/perlmutter/>`_.
