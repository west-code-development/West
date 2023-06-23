.. _summit:

===========
OLCF-Summit
===========

Summit is a GPU-accelerated supercomputer located at Oak Ridge National Laboratory, maintained by `OLCF <https://www.olcf.ornl.gov/>`_.

.. code-block:: bash

   $ ssh <username>@summit.olcf.ornl.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on August 9, 2022):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module load nvhpc/21.9
   module load cuda/11.0.3
   module load spectrum-mpi/10.4.0.3-20210112
   module load essl/6.3.0
   module load netlib-lapack/3.9.1
   module load python/3.8-anaconda3
   module load git/2.36.1
   module unload darshan-runtime

   export BLAS_LIBS="$OLCF_ESSL_ROOT/lib64/libessl.so"
   export LAPACK_LIBS="$OLCF_ESSL_ROOT/lib64/libessl.so $OLCF_NETLIB_LAPACK_ROOT/lib64/liblapack.a"

   ./configure --with-cuda=$OLCF_CUDA_ROOT --with-cuda-runtime=11.0 --with-cuda-cc=70

   # Manually edit make.inc: add -D__GPU_MPI to DFLAGS

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   sed -i 's/-L.*/-L\/autofs\/nccs-svm1_sw\/summit\/python\/3.8\/anaconda3\/2020.07-rhel8\/lib -lpython3.8/' west_make.inc
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs
~~~~~~~~~~~~~~~~~

Summit uses the `jsrun <https://docs.olcf.ornl.gov/systems/summit_user_guide.html#job-launcher-jsrun>`_ job manager. The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Summit with 6 MPI ranks and 6 GPUs per node. The <project_name> must be replaced with an active project allocation.

**Important**: It is recommended to run the calculation from the IBM Spectrum Scale file system (`$MEMBERWORK` instead of `/home`).

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash
   #BSUB -P <project_name>
   #BSUB -W 0:20
   #BSUB -nnodes 2
   #BSUB -J jobname
   #BSUB -o jobname.%J
   #BSUB -N
   #BSUB -q debug

   module load nvhpc/21.9
   module load cuda/11.0.3
   module load spectrum-mpi/10.4.0.3-20210112
   module load essl/6.3.0
   module load netlib-lapack/3.9.1
   module load python/3.8-anaconda3
   module unload darshan-runtime

   export OMP_NUM_THREADS=1

   # The following env vars improve MPI I/O performance

   export PAMI_ENABLE_STRIPING=1
   export PAMI_IBV_ADAPTER_AFFINITY=1
   export PAMI_IBV_DEVICE_NAME="mlx5_0:1,mlx5_3:1"
   export PAMI_IBV_DEVICE_NAME_1="mlx5_3:1,mlx5_0:1"

   export OMPI_MCA_io=romio321
   export ROMIO_HINTS=/path/to/romio_hints

   jsrun -n 4 -a 3 -c 3 -g 3 -r 2 --smpiargs="-gpu" ./wstat.x -i wstat.in &> wstat.out

The value of `-n` should be two times the number of nodes. When running QE and WEST, usually there is no need to change `-a`, `-c`, `-g`, and `-r`.

`romio_hints` is a text file with the following content:

.. code-block::

   romio_cb_write enable
   romio_ds_write enable
   cb_buffer_size 16777216
   cb_nodes 2

Job submission is done with the following:

.. code-block:: bash

   $ bsub run_west.sh

.. seealso::
   For more information, visit the `OLCF user guide <https://docs.olcf.ornl.gov/systems/summit_user_guide.html>`_.
