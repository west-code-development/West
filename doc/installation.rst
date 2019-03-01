.. _installation:

============
Installation
============

In order to install WEST you need to download the `QuantumEspresso 6.1 <https://gitlab.com/QEF/q-e/-/archive/qe-6.1.0/q-e-qe-6.1.0.tar>`_.

`QuantumEspresso <http://www.quantum-espresso.org/>`_ (QE) is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale, based on density-functional theory (DFT), plane waves (PW), and pseudopotentials (PP).
Configure QuantumEspresso by running the ``configure`` script that comes with the QE distribution. WEST requires `MPI <https://en.wikipedia.org/?title=Message_Passing_Interface>`_ support (`Scalapack <http://www.netlib.org/scalapack/>`_ and `OpenMP <http://openmp.org/>`_ support is also recommended, but optional). If all the environment variables (compilers, libraries etc.) have been set according to the QE configure guide, this would simply be:

.. code-block:: bash 

   $ git clone -b 'qe-6.1.0' --single-branch --depth 1 https://gitlab.com/QEF/q-e.git QEdir
   $ cd QEdir
   $ git clone -b 'v3.1.1' --single-branch --depth 1 http://greatfire.uchicago.edu/west-public/West.git West
   $ ./configure

It's now time to create the ``pw.x``, ``wstat.x``, ``wfreq.x``, and ``westpp.x`` executables by doing:

.. code-block:: bash 

   $ cd QEdir
   $ make pw
   $ cd QEdir/West
   $ make

You have succefully installed QuantumEspresso and WEST if you see the executables ``pw.x``, ``wstat.x``, ``wfreq.x``, and ``westpp.x`` created in the QEdir/bin directory.

.. code-block:: bash 

   $ ls QEdir/bin/
   pw.x
   wstat.x
   wfreq.x
   westpp.x
   ... (other content) ...

Congratulations, you are all set for running QuantumEspresso and WEST! 
