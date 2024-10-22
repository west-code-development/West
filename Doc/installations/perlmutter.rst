.. _perlmutter:

================
NERSC-Perlmutter
================

Perlmutter is an HPE Cray EX supercomputer located at National Energy Research Scientific Computing Center (`NERSC <https://www.nersc.gov/>`_). Perlmutter has both GPU-accelerated nodes and CPU-only nodes.

.. code-block:: bash

   $ ssh <username>@saul-p1.nersc.gov

Building WEST (GPU)
~~~~~~~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on October 22, 2024):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module unload darshan
   module load PrgEnv-nvidia
   module load nvidia/23.9
   module load cudatoolkit/12.2
   module load craype-accel-nvidia80
   module load cray-python/3.11.5

   ./configure --with-cuda=$CUDA_HOME --with-cuda-runtime=12.2 --with-cuda-cc=80 --with-cuda-mpi=yes

   # Manually edit make.inc:

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

Running WEST Jobs (GPU)
~~~~~~~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two GPU nodes of Perlmutter with 4 MPI ranks and 4 GPUs per node. The <project_name> must be replaced with an active project allocation.

**Important**: The following environment variable is needed to work around a bug in ROMIO, Cray MPICH.

.. code-block:: bash

   export ROMIO_FSTYPE_FORCE="ufs:"

**Important**: It is recommended to run the calculation from the Lustre file system (`$PSCRATCH` instead of `/home`).

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

   module unload darshan
   module load PrgEnv-nvidia
   module load nvidia/23.9
   module load cudatoolkit/12.2
   module load craype-accel-nvidia80
   module load cray-python/3.11.5

   export LD_LIBRARY_PATH=/opt/cray/pe/python/3.11.5/lib:$LD_LIBRARY_PATH
   export OMP_NUM_THREADS=1
   export SLURM_CPU_BIND=cores
   export MPICH_GPU_SUPPORT_ENABLED=1
   export ROMIO_FSTYPE_FORCE="ufs:"

   srun -n 8 ./wstat.x -i wstat.in &> wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

Building WEST (CPU)
~~~~~~~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on October 22, 2024):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module unload darshan
   module load cpu
   module load cray-fftw/3.3.10.6
   module load cray-python/3.11.5

   export CRAYPE_LINK_TYPE=dynamic
   export MPIF90=ftn
   export F90=ftn
   export CC=cc

   ./configure --enable-openmp --with-scalapack

   # Manually edit make.inc:

   # DFLAGS = -D__FFTW3 -D__MPI -D__SCALAPACK
   # IFLAGS = -I. -I$(TOPDIR)/include -I$(TOPDIR)/FoX/finclude -I/opt/cray/pe/fftw/3.3.10.6/x86_milan/include
   # BLAS_LIBS = # leave blank
   # LAPACK_LIBS = # leave blank

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs (CPU)
~~~~~~~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two CPU nodes of Perlmutter with 128 MPI ranks per node. The <project_name> must be replaced with an active project allocation.

**Important**: The following environment variable is needed to work around a bug in ROMIO, Cray MPICH.

.. code-block:: bash

   export ROMIO_FSTYPE_FORCE="ufs:"

**Important**: It is recommended to run the calculation from the Lustre file system (`$PSCRATCH` instead of `/home`).

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash

   #SBATCH --job-name=WEST
   #SBATCH --time=00:20:00
   #SBATCH --account=<project_name>
   #SBATCH --constraint=cpu
   #SBATCH --qos=debug
   #SBATCH --nodes=2
   #SBATCH --ntasks-per-node=128
   #SBATCH --cpus-per-task=2

   module unload darshan
   module load cpu
   module load cray-fftw/3.3.10.6
   module load cray-python/3.11.5

   export LD_LIBRARY_PATH=/opt/cray/pe/python/3.11.5/lib:$LD_LIBRARY_PATH
   export OMP_NUM_THREADS=1
   export SLURM_CPU_BIND=cores
   export ROMIO_FSTYPE_FORCE="ufs:"

   srun -n 256 ./wstat.x -i wstat.in &> wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the `NERSC user guide <https://docs.nersc.gov/systems/perlmutter/>`_.
