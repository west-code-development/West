.. _Manual:

Manual
======

The complete **WEST** reference for input parameters. 

.. contents:: :local:
              :depth: 1

.. note:: 
   Not all input parameters listed below are mandatory. Check :ref:`quickreference` or :ref:`tutorial` pages to see examples of input files. 

.. seealso::
   **WESTpy** is a Python package, designed to assist users of the WEST code in pre- and post-process massively parallel calculations. Click `here <http://www.west-code.org/doc/westpy/latest/>`_ to know more.

.. seealso::
   The input file is given according to the YAML Notation (`https://yaml.org/ <https://yaml.org//>`_).  

|

----------
input_west
----------

.. data:: qe_prefix

   :type: string 
   :default: "calc"
   :description: Prefix prepended to the QuantumEspresso save folder. 

   
.. data:: west_prefix

   :type: string 
   :default: "west"
   :description: Prefix prepended to the WEST save and restart folder. 
 
.. data:: outdir 

   :type: string 
   :default: Value of the ESPRESSO_TMPDIR environment variable if set; current directory ("./") otherwise
   :description: Input, temporary, output files are found in this directory.

|

-------------
wstat_control
-------------

.. data:: wstat_calculation

   :type: string 
   :default: "S"
   :description: Available options are:
 
      - "S" : Start from scratch
      - "R" : Restart from an interrupted run. You should restart with the same number of cores, and images. 
      - "E" : Calculation is outsourced to a server
   
.. data:: n_pdep_eigen

   :type: int 
   :default: 1
   :description: Number of PDEP eigenpotentials. 
   
.. data:: n_pdep_times

   :type: int 
   :default: 4
   :description: Maximum dimension of the search space = n_pdep_eigen * n_pdep_times. 
   
.. data:: n_pdep_maxiter

   :type: int 
   :default: 100
   :description: Maximum number of iterations in PDEP. 
   
.. data:: n_dfpt_maxiter

   :type: int 
   :default: 250
   :description: Maximum number of iterations in DFPT. 
   
.. data:: n_pdep_read_from_file

   :type: int 
   :default: 0
   :description: Number of PDEP eigenpotentials that can be read from file.  
   
.. data:: trev_pdep

   :type: float 
   :default: 0.001
   :description: Absolute convergence threshold in the PDEP eigenvalues.  
   
.. data:: trev_pdep_rel

   :type: float 
   :default: 0.1
   :description: Relative convergence threshold in the PDEP eigenvalues.  
   
.. data:: tr2_dfpt

   :type: float 
   :default: 1e-12
   :description: Convergence threshold in DFPT. Note that in the first PDEP iterations a reduced threshold for DFPT could be used by the code in order to speed up the computation.
   
.. data:: l_minimize_exx_if_active

   :type: boolean
   :default: False
   :description: If (True), then the exact-exchange term in the Hamiltonian is computed with the cutoff of the wavefunction.
   
.. data:: l_kinetic_only

   :type: boolean
   :default: False
   :description: If (True), then only the kinetic term in the Hamiltonian is kept.
   
.. data:: l_use_ecutrho 

   :type: boolean
   :default: False
   :description: If (True), then the eigenpotentials are represented with ecutrho instead of ecutwfc.
   
.. data:: qlist 

   :type: list of int
   :default: [1,2,...,number of q-points]
   :description: List of q-points to compute.

|

-------------
wfreq_control
-------------

.. data:: wfreq_calculation

   :type: string 
   :default: "XWGQ"
   :description: Available options are:
 
      - "XWGQ" : Compute the QP corrections.
      - "XwGQ" : Compute the QP corrections, restart from an interrupted / just read W run.
      - "XwgQ" : Compute the QP corrections, restart from an interrupted / just read G run.
      - "X" : Compute the HF corrections.
      - "XWO" : Compute the optical properties.
      - "XWGQP" : Compute the QP corrections, and plot spectral functions.
      - "XWGQOP" : Compute all.
                    
.. data:: n_pdep_eigen_to_use

   :type: int
   :default: 2
   :description: Number of PDEP eigenvectors to use in Wfreq. They are read from previous Wstat run. This value cannot exceed n_pdep_eigen (defined in wstat_control) and is used to check the convergence of the calculation.

.. data:: qp_bandrange

   :type: list of int
   :default: [1,2]
   :description: Compute the QP corrections from band qp_bandrange[0] to band qp_bandrange[1].

.. data:: macropol_calculation

   :type: string
   :default: "N"
   :description: Available options are:
   
      - "N" : None. Choice valid for isolated systems.
      - "C" : Include long-wavelength limit. Choice valid for condensed systems.

.. data:: n_lanczos

   :type: int
   :default: 30
   :description: Number of Lanczos chains.

.. data:: n_imfreq

   :type: int
   :default: 128
   :description: Number of frequecies used to sample the imaginary frequency axis in the range [0,ecut_imfreq].

.. data:: n_refreq

   :type: int
   :default: 10
   :description: Number of frequecies used to sample the real frequency axis in the range [0,ecut_refreq].

.. data:: ecut_imfreq

   :type: float
   :default: Cut of the density, read from the ground state
   :description: Cutoff for the imaginary frequencies (in Ry).

.. data:: ecut_refreq

   :type: float
   :default: 2.0
   :description: Cutoff for the real frequencies (in Ry).

.. data:: wfreq_eta

   :type: float
   :default: 0.003675
   :description: Energy shift of the poles (in Ry). 

.. data:: n_secant_maxiter

   :type: int
   :default: 1
   :description: Maximum number of iterations in the secant solver.

.. data:: trev_secant

   :type: float
   :default: 0.003675
   :description: Convergence energy threshold (in Ry) for the secant solver.

.. data:: l_enable_lanczos

   :type: boolean
   :default: True
   :description: If (False), then Lanczos solvers are turned off.

.. data:: l_enable_gwetot

   :type: boolean
   :default: False
   :description: Deprecated parameter.

.. data:: o_restart_time

   :type: float
   :default: 0.0
   :description: Available options are:

      - If ( o_restart_time == 0 ) A checkpoint is written at every iteration of the W and G loops.
      - If ( o_restart_time >  0 ) A checkpoint is written every o_restart_time minutes in the W and G loops.
      - If ( o_restart_time <  0 ) A checkpoint is NEVER written in the W and G loops. Restart will not be possible.

.. data:: ecut_spectralf

   :type: list of float
   :default: [-2.0,2.0]
   :description: Energy cutoff (in Ry) for the real frequencies. Used when wfreq_caculation contains the runlevel "P".

.. data:: n_spectralf

   :type: int
   :default: 10
   :description: Number of frequecies used to plot the spectral function (runlevel "P"), sampling the interval [-ecut_spectralf[0],ecut_spectralf[1]].

|

--------------
westpp_control
--------------

.. data:: westpp_calculation

   :type: string 
   :default: "R"
   :description: Available options are:

      - "R" : Output rho, the electronic density.
      - "W" : Output the electronic wavefunctions.
      - "E" : Output the eigenpotentials.
      - "S" : Output the screened exchange constant.

.. data:: westpp_range

   :type: list of int 
   :default: [1,2]
   :description: Range for W, E, and S run.

.. data:: westpp_format

   :type: string 
   :default: "C"
   :description: Available options for the output fortmat are:
          
      - "c" : Cube.
      - "x" : Planar average yz.
      - "y" : Planar average xz.
      - "z" : Planar average xy.
      - "s" : Spherical average.

.. data:: westpp_sign

   :type: boolean
   :default: False
   :description: If (True), then the sign of the wavefunction/eigenpotential is kept in the output file.

.. data:: westpp_n_pdep_eigen_to_use

   :type: int
   :default: 1
   :description: Number PDEP eigenpotentials to read/use.

.. data:: westpp_r0

   :type: 3-dim list of floats (a vector)
   :default: [0.0, 0.0, 0.0]
   :description: Position of the center (in a.u.) for spherical average plot.

.. data:: westpp_nr

   :type: int
   :default: 100
   :description: Number of points in the spherical average plot.

.. data:: westpp_rmax

   :type: float
   :default: 1.0
   :description: Max radius (in a.u.) for the spherical average plot.

.. data:: westpp_epsinfty

   :type: float
   :default: 1.0
   :description: Macroscopic relative dielectric constant. Used in the "S" runlevel.
            

