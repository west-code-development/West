.. _dgx:

===============
NVIDIA DGX A100
===============

The following instructions have been tested on an NVIDIA DGX A100 machine.

Requirements:

- NVIDIA HPC SDK (21.9 or newer)
- Python3

To download and install NVIDIA HPC SDK, do:

.. code-block:: bash

   $ wget https://developer.download.nvidia.com/hpc-sdk/21.9/nvhpc_2021_219_Linux_x86_64_cuda_multi.tar.gz
   $ tar xpzf nvhpc_2021_219_Linux_x86_64_cuda_multi.tar.gz
   $ nvhpc_2021_219_Linux_x86_64_cuda_multi/install

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script:

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

   make conf PYT=python3 PYT_LDFLAGS="-L/usr/lib/python3.8/config-3.8-x86_64-linux-gnu -lpython3.8"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST
~~~~~~~~~~~~

We can run the `wstat.x` WEST executables on 2 GPUs using the following command:

.. code-block:: bash

   $ export OMP_NUM_THREADS=1
   $ mpirun -np 2 ./wstat.x -i wstat.in > wstat.out
