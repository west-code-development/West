.. _overview:

Overview
========

**WEST** (Without Empty STates) is a massively parallel software for large scale electronic structure calculations within many-body perturbation theory.

Features:

   - GW self-energy (full-frequency) [*GPU enabled*]
   - Absorption spectra and energies with TDDFT or BSE [*GPU enabled*]
   - Excited states forces with TDDFT [*GPU enabled*]
   - CI-in-DFT quantum defect embedding theory (QDET) [*GPU enabled*]
   - Electron-phonon *under development*

WEST uses as starting point the charge density and Kohn-Sham orbitals obtained with hybrid or semilocal DFT. WEST uses the plane-wave pseudopotential method.

.. seealso::
   **WESTpy** is a Python package, designed to assist users of the WEST code in pre- and post-process massively parallel calculations. Click `here <https://west-code.org/doc/westpy/latest/>`_ to know more.
