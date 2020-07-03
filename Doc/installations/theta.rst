.. _theta:

=================
Theta-ALCF (XC40)
=================

Theta is a Cray XC40 located at Argonne National Laboratory, maintained by `ALCF <https://www.alcf.anl.gov/>`_. 

.. code-block:: bash 

   $ ssh -Y <username>@theta.alcf.anl.gov

Building WEST
~~~~~~~~~~~~~

WEST executables can be compiled using the following script: 

.. code-block:: bash 

   $ cat build_west.sh
   #!/bin/bash

   module unload cray-libsci
   module load cray-python/3.6.5.3
   export CRAYPE_LINK_TYPE=dynamic
   export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/python/3.6.5.3/lib:$LD_LIBRARY_PATH

   ./configure MPIF90=ftn CC=cc CXX=CC --enable-openmp=yes --enable-parallel=yes --enable-shared=yes --with-scalapack=intel SCALAPACK_LIBS="${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.so -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.so ${MKLROOT}/lib/intel64/libmkl_intel_thread.so ${MKLROOT}/lib/intel64/libmkl_core.so ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.so -Wl,--end-group" FFLAGS=" -xMIC-AVX512 -qopenmp -align array64byte -fp-model fast=2 -no-prec-div -assume byterecl" --with-hdf5=no CFLAGS=" -xMIC-AVX512" LDFLAGS=" -shared-intel -qopenmp" 

   make pw -j 16

   cd West
   make conf PYT=python3
   make all 

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
   module load cray-python/3.6.5.3
   export CRAYPE_LINK_TYPE=dynamic
   export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018.0.128/linux/compiler/lib/intel64:$LD_LIBRARY_PATH
   export LD_LIBRARY_PATH=/opt/python/3.6.5.3/lib:$LD_LIBRARY_PATH

   echo "Running Cobalt Job $COBALT_JOBID."

   export OMP_NUM_THREADS=$NTHREADS
   aprun -n $MPIRANKS -N $MPIRANKS_PERNODE -cc depth -d $NTHREADS -j $HT ./wstat.x -i wstat.in &> wstat.out

Make the script executable: 

.. code-block:: bash 

   $ chmod +x run_west.sh

Job submission is done with the following: 

.. code-block:: bash 

   $ qsub run_west.sh

.. seealso::
   For more information, visit the ALCF user guide (`https://www.alcf.anl.gov/user-guides/xc40-system-overview <https://www.alcf.anl.gov/user-guides/xc40-system-overview/>`_).
