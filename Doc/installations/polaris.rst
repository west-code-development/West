.. _polaris:

============
ALCF-Polaris
============

Polaris is a GPU-accelerated supercomputer located at Argonne National Laboratory, maintained by `ALCF <https://www.alcf.anl.gov/>`_.

.. code-block:: bash

   $ ssh <username>@polaris.alcf.anl.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on September 17, 2023):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module load nvhpc/23.3
   module load cray-libsci/21.08.1.2
   module load cray-python/3.9.12.1

   export MPICH_GPU_SUPPORT_ENABLED=1
   export CUDA_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/cuda/11.8

   ./configure --with-cuda=$CUDA_HOME --with-cuda-runtime=11.8 --with-cuda-cc=80 --with-cuda-mpi=yes

   # Manually edit make.inc:

   # MPIF90 = ftn
   # F90 = ftn
   # CC = cc
   # LD = ftn
   # BLAS_LIBS = # leave blank
   # LAPACK_LIBS = # leave blank

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="-L/opt/cray/pe/python/3.9.12.1/lib -lpython3.9"
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

   module load nvhpc/23.3
   module load cray-libsci/21.08.1.2
   module load cray-python/3.9.12.1

   export MPICH_GPU_SUPPORT_ENABLED=1
   export ROMIO_FSTYPE_FORCE="ufs:"

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
