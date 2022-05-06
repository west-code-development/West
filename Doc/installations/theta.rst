.. _theta:

=================
Theta-ALCF (XC40)
=================

Theta is a Cray XC40 located at Argonne National Laboratory, maintained by `ALCF <https://www.alcf.anl.gov/>`_.

.. code-block:: bash

   $ ssh <username>@theta.alcf.anl.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script (tested on May 6, 2022):

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
   export SCALAPACK_LIBS="$MKLROOT/lib/intel64/libmkl_scalapack_lp64.so -Wl,--start-group $MKLROOT/lib/intel64/libmkl_intel_lp64.so $MKLROOT/lib/intel64/libmkl_intel_thread.so $MKLROOT/lib/intel64/libmkl_core.so $MKLROOT/lib/intel64/libmkl_blacs_intelmpi_lp64.so -Wl,--end-group"

   ./configure --enable-openmp --with-scalapack=intel

   make -j 8 pw

   cd West

   make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   make -j 8 all

To use the script do:

.. code-block:: bash

   $ bash build_west.sh

Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Theta with 64 MPI ranks per node. The <project_name> must be replaced with an active project allocation.

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

   echo "Running Cobalt Job $COBALT_JOBID."

   export OMP_NUM_THREADS=$NTHREADS
   aprun -n $MPIRANKS -N $MPIRANKS_PERNODE -cc depth -d $NTHREADS -j $HT ./wstat.x -i wstat.in &> wstat.out

Make the script executable:

.. code-block:: bash

   $ chmod 755 run_west.sh

Job submission is done with the following:

.. code-block:: bash

   $ qsub run_west.sh

.. seealso::
   For more information, visit the `ALCF user guide <https://www.alcf.anl.gov/user-guides/xc40-system-overview/>`_.
