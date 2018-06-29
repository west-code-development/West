.. _quickreference:

Quick Reference
===============

These are quick references for **WEST** input file examples. 

wstat.x
~~~~~~~

This is a typical input for ``wstat.x``. 

.. code-block:: bash 

   {
     "input_west": {
       "qe_prefix": "silane",
       "west_prefix": "silane",
       "outdir": "./"
     },
     "wstat_control": {
       "wstat_calculation": "S",
       "n_pdep_eigen": 50
     }
   }

.. note:: 
   Indentation is not important. 

.. seealso:: 
   All input keywords are referenced and explained in the Manual. 

wfreq.x
~~~~~~~

This is a typical input for ``wfreq.x``. 

.. code-block:: bash 

   {
     "input_west": {
       "qe_prefix": "silane",
       "west_prefix": "silane",
       "outdir": "./"
     },
     "wstat_control": {
       "wstat_calculation": "S",
       "n_pdep_eigen": 50
     },
     "wfreq_control": {
        "wfreq_calculation": "XWGQ",
        "n_pdep_eigen_to_use": 50,
        "qp_bandrange": [1,5],
        "ecut_refreq": 2.0,
        "n_refreq": 300
     }
   }

.. note:: 
   Indentation is not important. 

.. seealso:: 
   All input keywords are referenced and explained in the Manual. 
