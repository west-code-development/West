.. _polaris:

============
ALCF-Polaris
============

Polaris is a GPU-accelerated supercomputer located at Argonne National Laboratory, maintained by `ALCF <https://www.alcf.anl.gov/>`_.

.. code-block:: bash

   $ ssh <username>@polaris.alcf.anl.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on May 8, 2024):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module load craype-accel-nvidia80
   module load nvhpc/23.9
   module load cray-libsci/23.12.5
   module load cray-python/3.11.5

   export MPICH_GPU_SUPPORT_ENABLED=1

   ./configure --with-cuda=$NVIDIA_PATH/cuda/12.2 --with-cuda-runtime=12.2 --with-cuda-cc=80 --with-cuda-mpi=yes

   # Manually edit make.inc:

   # MPIF90 = ftn
   # F90 = ftn
   # CC = cc
   # LD = ftn
   # BLAS_LIBS = # leave blank
   # LAPACK_LIBS = # leave blank

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="$PYTHON_PATH/lib/libpython3.11.so"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Polaris with 4 MPI ranks and 4 GPUs per node. The <project_name> must be replaced with an active project allocation.

**Important**: The following environment variable is needed to work around a bug in ROMIO, Cray MPICH.

.. code-block:: bash

   export ROMIO_FSTYPE_FORCE="ufs:"

**Important**: It is recommended to run the calculation from one of the Lustre file systems (`/grand` or `/eagle` instead of `/home`).

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash -l
   #PBS -l select=1:system=polaris
   #PBS -l place=scatter
   #PBS -l walltime=0:20:00
   #PBS -l filesystems=home:grand
   #PBS -j oe
   #PBS -q debug
   #PBS -A <project_name>
   #PBS -N job_name

   module load craype-accel-nvidia80
   module load nvhpc/23.9
   module load cray-libsci/23.12.5
   module load cray-python/3.11.5

   export MPICH_GPU_SUPPORT_ENABLED=1
   export ROMIO_FSTYPE_FORCE="ufs:"
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHON_PATH/lib

   NNODES=`wc -l < $PBS_NODEFILE`
   NRANKS_PER_NODE=$(nvidia-smi -L | wc -l)
   NDEPTH=8
   NTHREADS=1
   NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))

   cd ${PBS_O_WORKDIR}

   mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --depth=${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS=${NTHREADS} -env OMP_PLACES=threads ./wstat.x -i wstat.in &> wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ qsub run_west.sh

.. seealso::
   For more information, visit the `ALCF user guide <https://docs.alcf.anl.gov/polaris/getting-started/>`_.
