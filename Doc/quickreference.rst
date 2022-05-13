.. _quickreference:

Quick Reference
===============

These are quick references for **WEST** input file examples.

.. seealso::
   All input keywords are referenced and explained in the :ref:`manual`.

wstat.x
~~~~~~~

This is a typical input for ``wstat.x``.

.. code-block:: yaml

   input_west:
      qe_prefix: silane
      west_prefix: silane
      outdir: ./

   wstat_control:
      wstat_calculation: S
      n_pdep_eigen: 50

wfreq.x
~~~~~~~

This is a typical input for ``wfreq.x``.

.. code-block:: yaml

   input_west:
      qe_prefix: silane
      west_prefix: silane
      outdir: ./

   wstat_control:
      wstat_calculation: S
      n_pdep_eigen: 50

   wfreq_control:
      wfreq_calculation: XWGQ
      n_pdep_eigen_to_use: 50
      qp_bandrange: [1,5]
      ecut_refreq: 2.0
      n_refreq: 300

westpp.x
~~~~~~~~

This is a typical input for ``westpp.x``.

.. code-block:: yaml

   input_west:
      qe_prefix: silane
      west_prefix: silane
      outdir: ./

   wstat_control:
      wstat_calculation: S
      n_pdep_eigen: 50

   westpp_control:
      westpp_calculation: E
      westpp_range: [1,2]
      westpp_format: C
      westpp_sign: True
