.. _theta:

=================
Theta-ALCF (XC40)
=================

Theta is a Cray XC40 located at Argonne National Laboratory, maintained by `ALCF <https://www.alcf.anl.gov/>`_. 

.. code-block:: bash 

   $ ssh -XY <username>@theta.alcf.anl.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script: 

.. code-block:: bash 

   $ cat build_west.sh
   #!/bin/bash

   module unload cray-libsci
   module load cray-python/3.6.5.3

   export CRAYPE_LINK_TYPE=dynamic
   export LD_LIBRARY_PATH=/opt/python/3.6.5.3/lib:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64:$LD_LIBRARY_PATH
   
   ./configure MPIF90=ftn CC=cc --enable-openmp --with-scalapack=intel LD_LIBS="`python3-config --ldflags`" SCALAPACK_LIBS="${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group"

   make pw -j 16

   cd West
   make

To use the script do: 

.. code-block:: bash 

   $ bash build_west.sh

Running WEST Jobs
~~~~~~~~~~~~~~~~~

The following is an example executable script `run_west.sh` to run the `wstat.x` WEST executable on two nodes of Theta with 64 MPI ranks per node. The <project_name> must be replaced with an active project allocation.

.. code-block:: bash 

   $ cat run_west.sh
   #!/bin/bash
   #COBALT -n 2 -t 10 -q debug-cache-quad -A <project_name> -O WEST

   module unload cray-libsci
   module load cray-python/3.6.5.3

   export CRAYPE_LINK_TYPE=dynamic
   export LD_LIBRARY_PATH=/opt/python/3.6.5.3/lib:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64:$LD_LIBRARY_PATH

   aprun -n 128 -N 64 -d 1 --cc depth -e OMP_NUM_THREADS=1 -j 1 ./wstat.x -i wstat.in > wstat.out

Make the script executable: 

.. code-block:: bash 

   $ chmod +x run_west.sh

Job submission is done with the following: 

.. code-block:: bash 

   $ qsub run_west.sh

.. seealso::
   For more information, visit the ALCF user guide (`https://www.alcf.anl.gov/user-guides/xc40-system-overview <https://www.alcf.anl.gov/user-guides/xc40-system-overview/>`_).
