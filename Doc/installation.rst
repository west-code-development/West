.. _installation:

============
Installation
============

In order to install WEST you need to download `Quantum ESPRESSO 6.8 <https://gitlab.com/QEF/q-e/-/archive/qe-6.8/q-e-qe-6.8.tar>`_.

To compute absorption spectra (BSE), you also need to download and install `Qbox <http://qboxcode.org>`_.

`Quantum ESPRESSO <http://www.quantum-espresso.org/>`_ (QE) is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale, based on density-functional theory (DFT), plane waves (PW), and pseudopotentials (PP).
Configure Quantum ESPRESSO by running the ``configure`` script that comes with the QE distribution. WEST requires `MPI <https://en.wikipedia.org/?title=Message_Passing_Interface>`_ support (`ScaLAPACK <http://www.netlib.org/scalapack/>`_ and `OpenMP <http://openmp.org/>`_ support is also recommended, but optional). If all the environment variables (compilers, libraries etc.) have been set according to the QE configure guide, this would simply be:

.. code-block:: bash

   $ git clone -b 'qe-6.8' --single-branch --depth 1 https://gitlab.com/QEF/q-e.git QEdir
   $ cd QEdir
   $ git clone -b 'v4.3.0' --single-branch --depth 1 http://greatfire.uchicago.edu/west-public/West.git West
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

You have succefully installed Quantum ESPRESSO and WEST if you see the executables ``pw.x``, ``wstat.x``, ``wfreq.x``, and ``westpp.x`` created in the QEdir/bin directory.

.. code-block:: bash

   $ ls QEdir/bin/
   pw.x
   wstat.x
   wfreq.x
   westpp.x
   ... (other content) ...

Congratulations, you are all set for running Quantum ESPRESSO and WEST!


Suggested configuration options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
   :maxdepth: 1

   installations/cori.rst
   installations/theta.rst
   installations/midway3.rst
   installations/macos.rst
