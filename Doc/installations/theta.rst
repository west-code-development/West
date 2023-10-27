.. _theta:

==========
ALCF-Theta
==========

Theta (**retired** on December 31, 2023) was a Cray XC40 located at Argonne National Laboratory, maintained by `ALCF <https://www.alcf.anl.gov/>`_.

.. code-block:: bash

   $ ssh <username>@theta.alcf.anl.gov

Building WEST
~~~~~~~~~~~~~

Start an interactive job and ssh to a compute node. The <project_name> must be replaced with an active project allocation.

.. code-block:: bash

   # Start interactive job
   $ qsub -I -n 1 -t 1:00:00 --attrs enable_ssh=1 -q debug-cache-quad -A <project_name>
     Connecting to thetamom1 for interactive qsub...
     Job routed to queue "debug-cache-quad".
     Memory mode set to cache quad for queue debug-cache-quad
     Wait for job 266815 to start...
     Opening interactive session to 3835

   # Get compute node number
   $ echo $COBALT_PARTNAME
     3835

   # Full name of compute node is nid + 5-digit node number
   $ ssh nid03835

   # Set up proxy for internet access
   $ export HTTP_PROXY=http://theta-proxy.tmi.alcf.anl.gov:3128
   $ export HTTPS_PROXY=http://theta-proxy.tmi.alcf.anl.gov:3128
   $ export http_proxy=http://theta-proxy.tmi.alcf.anl.gov:3128
   $ export https_proxy=http://theta-proxy.tmi.alcf.anl.gov:3128

**Important**: The `ELPA eigensolver library <https://elpa.mpcdf.mpg.de/>`_ is found to greatly improve the performance of Quantum ESPRESSO on Theta. ELPA can be installed following these steps (tested on February 3, 2023):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module unload cray-libsci
   module load cray-python/3.8.2.1

   export CRAYPE_LINK_TYPE=dynamic
   export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2020.0.166/linux/compiler/lib/intel64:$LD_LIBRARY_PATH
   export FC=ftn
   export CC=cc
   export CXX=CC

   wget https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/2022.11.001/elpa-2022.11.001.tar.gz
   tar zxf elpa-2022.11.001.tar.gz
   cd elpa-2022.11.001

   mkdir build
   cd build

   ../configure --prefix=$(pwd) LDFLAGS="$MKLROOT/lib/intel64/libmkl_scalapack_lp64.so -Wl,--start-group $MKLROOT/lib/intel64/libmkl_intel_lp64.so $MKLROOT/lib/intel64/libmkl_sequential.so $MKLROOT/lib/intel64/libmkl_core.so $MKLROOT/lib/intel64/libmkl_blacs_intelmpi_lp64.so -Wl,--end-group" --disable-sse-assembly --disable-sse --disable-avx512 --enable-c-tests=no

   make -j 8
   make install

WEST executables can be compiled using the following script (tested on February 3, 2023):

.. code-block:: bash

   $ cat build_west.sh
   #!/bin/bash

   module unload cray-libsci
   module load cray-python/3.8.2.1

   export CRAYPE_LINK_TYPE=dynamic
   export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2020.0.166/linux/compiler/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/python/3.8.2.1/lib:$LD_LIBRARY_PATH
   export MPIF90=ftn
   export F90=ftn
   export CC=cc
   export DFLAGS="-D__DFTI -D__MPI -D__SCALAPACK -D__ELPA"
   export BLAS_LIBS="-Wl,--start-group $MKLROOT/lib/intel64/libmkl_intel_lp64.so $MKLROOT/lib/intel64/libmkl_intel_thread.so $MKLROOT/lib/intel64/libmkl_core.so -Wl,--end-group"
   export LAPACK_LIBS="-Wl,--start-group $MKLROOT/lib/intel64/libmkl_intel_lp64.so $MKLROOT/lib/intel64/libmkl_intel_thread.so $MKLROOT/lib/intel64/libmkl_core.so -Wl,--end-group"
   export SCALAPACK_LIBS="$MKLROOT/lib/intel64/libmkl_scalapack_lp64.so $MKLROOT/lib/intel64/libmkl_blacs_intelmpi_lp64.so"

   # Edit ELPA installation path
   ./configure --enable-openmp --with-elpa-include=/path/to/elpa-2022.11.001/build/include/elpa-2022.11.001/modules --with-elpa-lib=/path/to/elpa-2022.11.001/build/lib/libelpa.a

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="-L/opt/python/3.8.2.1/lib -lpython3.8"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Theta with 64 MPI ranks per node. The <project_name> must be replaced with an active project allocation.

**Important**: The following environment variable is needed to work around a bug in ROMIO, Cray MPICH.

.. code-block:: bash

   export ROMIO_FSTYPE_FORCE="ufs:"

**Important**: It is recommended to run the calculation from one of the Lustre file systems (`/grand` or `/eagle` instead of `/home`).

.. code-block:: bash

   $ cat run_west.sh
   #!/bin/bash
   #COBALT -n 2
   #COBALT -t 00:20:00
   #COBALT -q debug-cache-quad
   #COBALT -A <project_name>
   #COBALT -O WEST

   MPIRANKS_PERNODE=64
   MPIRANKS=$((COBALT_PARTSIZE * MPIRANKS_PERNODE))
   NTHREADS=1
   HT=1

   module unload cray-libsci
   module load cray-python/3.8.2.1

   export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2020.0.166/linux/compiler/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/python/3.8.2.1/lib:$LD_LIBRARY_PATH

   export ROMIO_FSTYPE_FORCE="ufs:"

   export OMP_NUM_THREADS=$NTHREADS
   aprun -n $MPIRANKS -N $MPIRANKS_PERNODE -cc depth -d $NTHREADS -j $HT ./wstat.x -i wstat.in &> wstat.out

Make the script executable:

.. code-block:: bash

   $ chmod 755 run_west.sh

Job submission is done with the following:

.. code-block:: bash

   $ qsub run_west.sh

.. seealso::
   For more information, visit the `ALCF user guide <https://docs.alcf.anl.gov/theta/hardware-overview/machine-overview/>`_.
