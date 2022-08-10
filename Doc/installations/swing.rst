.. _swing:

=========
Swing-ANL
=========

Swing is an HPC cluster maintained by the `Laboratory Computing Resource Center (LCRC) <https://www.lcrc.anl.gov/>`_ at Argonne National Laboratory.

.. code-block:: bash

   $ ssh <username>@swing.lcrc.anl.gov

Building WEST
~~~~~~~~~~~~~

To download and install NVIDIA HPC SDK, do:

.. code-block:: bash

   $ wget https://developer.download.nvidia.com/hpc-sdk/21.9/nvhpc_2021_219_Linux_x86_64_cuda_multi.tar.gz
   $ tar xpzf nvhpc_2021_219_Linux_x86_64_cuda_multi.tar.gz
   $ nvhpc_2021_219_Linux_x86_64_cuda_multi/install

WEST executables can be compiled using the following script (tested on Aug 10, 2022):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   export LD_LIBRARY_PATH=/path/to/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/path/to/nvidia/hpc_sdk/Linux_x86_64/21.9/comm_libs/openmpi4/openmpi-4.0.5/lib:$LD_LIBRARY_PATH
   export PATH=/path/to/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/bin:$PATH
   export PATH=/path/to/nvidia/hpc_sdk/Linux_x86_64/21.9/comm_libs/openmpi4/openmpi-4.0.5/bin:$PATH
   export SCALAPACK_LIBS=/path/to/nvidia/hpc_sdk/Linux_x86_64/21.9/comm_libs/openmpi4/openmpi-4.0.5/lib/libscalapack.a

   ./configure --with-cuda=/path/to/nvidia/hpc_sdk/Linux_x86_64/21.9/cuda/11.0 --with-cuda-cc=80 --with-cuda-runtime=11.0

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   sed -i 's/ -L\/usr\/lib .*/ -lpython3.8/' ./west_make.inc
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh


Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on one node of Swing with 8 MPI ranks and 8 GPUs. The <account_name> must be replaced with an active project allocation.

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash
   #SBATCH --time=00:20:00
   #SBATCH --account=<account_name>
   #SBATCH --nodes=1
   #SBATCH --gres=gpu:8

   export LD_LIBRARY_PATH=/path/to/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/lib:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/path/to/nvidia/hpc_sdk/Linux_x86_64/21.9/comm_libs/openmpi4/openmpi-4.0.5/lib:$LD_LIBRARY_PATH
   export PATH=/path/to/nvidia/hpc_sdk/Linux_x86_64/21.9/comm_libs/openmpi4/openmpi-4.0.5/bin:$PATH

   export OMP_NUM_THREADS=1

   mpirun -n 8 ./wstat.x -i wstat.in > wstat.out

Job submission is done with the following:

.. code-block:: bash

   $ sbatch run_west.sh

.. seealso::
   For more information, visit the `LCRC user guide <https://www.lcrc.anl.gov/for-users/using-lcrc/running-jobs/running-jobs-on-swing/>`_.
