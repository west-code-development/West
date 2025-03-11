.. _installation:

============
Installation
============

In order to install WEST you need to download `Quantum ESPRESSO 7.4.1 <https://gitlab.com/QEF/q-e/-/tags>`_.

`Quantum ESPRESSO <http://www.quantum-espresso.org/>`_ (QE) is an integrated suite of open-source computer codes for electronic-structure calculations and materials modeling at the nanoscale, based on density-functional theory (DFT), plane waves (PW), and pseudopotentials (PP).

QE can be installed with HDF5 support. Currently the installation of QE with CMake is not supported by WEST. Configure QE by running the ``configure`` script that comes with the QE distribution. WEST requires `MPI <https://en.wikipedia.org/?title=Message_Passing_Interface>`_ support. `HDF5 <https://www.hdfgroup.org/solutions/hdf5/>`_ is recommended, but optional. For large-scale calculations, `ScaLAPACK <http://www.netlib.org/scalapack/>`_ and `ELPA <https://elpa.mpcdf.mpg.de/>`_ are also recommended. If all the environment variables (compilers, libraries etc.) have been set according to the QE configure guide, this would simply be:

.. code-block:: bash

   $ git clone -b 'qe-7.4.1' --single-branch --depth 1 https://gitlab.com/QEF/q-e.git QEdir
   $ cd QEdir
   $ git clone -b 'v6.1.0' --single-branch --depth 1 https://github.com/west-code-development/West.git West
   $ ./configure

.. note::
   Note that since v4.0.0 WEST requires dynamic linking and Python3.

It's now time to create the ``pw.x``, ``wstat.x``, ``wfreq.x``, and ``westpp.x`` executables by doing:

.. code-block:: bash

   $ cd QEdir
   $ make pw
   $ cd QEdir/West
   $ make conf PYT=python3 PYT_LDFLAGS="`python3-config --ldflags --embed`"
   $ make all

You have succefully installed QE and WEST if you see the executables ``pw.x``, ``wstat.x``, ``wfreq.x``, ``wbse.x``, and ``westpp.x`` created in the QEdir/bin directory.

.. code-block:: bash

   $ ls QEdir/bin/
   pw.x
   wstat.x
   wfreq.x
   wbse.x
   westpp.x
   ... (other content) ...

Congratulations, you are all set for running QE and WEST!


Suggested configuration options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   installations/polaris.rst
   installations/theta.rst
   installations/bebop.rst
   installations/swing.rst
   installations/leonardo.rst
   installations/macos.rst
   installations/cori.rst
   installations/perlmutter.rst
   installations/dgx.rst
   installations/summit.rst
   installations/midway2.rst
   installations/midway3.rst
