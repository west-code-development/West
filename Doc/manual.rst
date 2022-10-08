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

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "pwscf"
   * - **Description**
     - Prefix prepended to the Quantum ESPRESSO save folder.


.. data:: west_prefix

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "west"
   * - **Description**
     - Prefix prepended to the WEST save and restart folders.


.. data:: outdir

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "./"
   * - **Description**
     - Directory for: input, temporary, and output files.


|

-------------
wstat_control
-------------

.. data:: wstat_calculation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "S"
   * - **Description**
     - Available options are:

       - "S" : Start from scratch.
       - "R" : Restart from an interrupted run. You should restart with the same number of cores, and images.
       - "E" : Calculation of the response is external, i.e. outsourced to a server.


.. data:: n_pdep_eigen

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - dynamically set to match the number of electrons
   * - **Description**
     - Number of PDEP eigenpotentials.


.. data:: n_pdep_times

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 4
   * - **Description**
     - Maximum dimension of the search space = n_pdep_eigen * n_pdep_times.


.. data:: n_pdep_maxiter

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 100
   * - **Description**
     - Maximum number of iterations in PDEP.


.. data:: n_dfpt_maxiter

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 250
   * - **Description**
     - Maximum number of iterations in DFPT.


.. data:: n_pdep_read_from_file

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 0
   * - **Description**
     - Number of PDEP eigenpotentials that can be read from file.


.. data:: n_steps_write_restart

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 1
   * - **Description**
     - Available options are:

       - If ( n_steps_write_restart >  0 ) A checkpoint is written every n_steps_write_restart iterations in the PDEP loop.
       - If ( n_steps_write_restart <= 0 ) A checkpoint is NEVER written in the PDEP loop. Restart will not be possible.


.. data:: trev_pdep

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.001
   * - **Description**
     - Absolute convergence threshold for PDEP eigenvalues.


.. data:: trev_pdep_rel

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.1
   * - **Description**
     - Relative convergence threshold for PDEP eigenvalues.


.. data:: tr2_dfpt

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1.e-12
   * - **Description**
     - Convergence threshold in DFPT. Note that in the first PDEP iterations a reduced threshold for DFPT could be used by the code in order to speed up the computation.


.. data:: l_kinetic_only

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then only the kinetic term in the Hamiltonian is kept.


.. data:: l_minimize_exx_if_active

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then the exact-exchange term in the Hamiltonian is computed with the cutoff of the wavefunction.


.. data:: l_use_ecutrho

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then the eigenpotentials are represented with ecutrho instead of ecutwfc.


.. data:: qlist

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - list of int
   * - **Default**
     - [1,2,...,number of q-points]
   * - **Description**
     - List of q-points to compute.


|

-------------
wfreq_control
-------------


.. data:: wfreq_calculation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "XWGQ"
   * - **Description**
     - Available options are:

       - "XWGQ" : Compute the QP corrections.
       - "XwGQ" : Compute the QP corrections, restart from an interrupted / just read W run.
       - "XwgQ" : Compute the QP corrections, restart from an interrupted / just read G run.
       - "XWGQH" : Compute the QP corrections and parameters of QDET effective Hamiltonian. Only available for Gamma-point sampling.
       - "XwGQH" : Compute the QP corrections and parameters of QDET effective Hamiltonian, restart from interrupted / just read W run. Only available for Gamma-point sampling.
       - "X" : Compute the HF corrections.
       - "XWO" : Compute the optical properties.
       - "XWGQP" : Compute the QP corrections, and plot spectral functions.
       - "XWGQOP" : Compute all.


.. data:: n_pdep_eigen_to_use

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - dynamically set to match the number of electrons
   * - **Description**
     - Number of PDEP eigenvectors to use in Wfreq. They are read from previous Wstat run. This value cannot exceed n_pdep_eigen (defined in wstat_control) and is used to check the convergence of the calculation.


.. data:: qp_bandrange

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - 2-dim list of int
   * - **Default**
     - [1,2]
   * - **Description**
     - Compute the QP corrections from band qp_bandrange[0] to band qp_bandrange[1]. Used only when qp_bands is not set. If qp_bands is set, the value of qp_bandrange is discarded.


.. data:: qp_bands

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - list of int
   * - **Default**
     - [0]
   * - **Description**
     - List of bands to compute the QP corrections. If qp_bands is not set, qp_bands is determined from qp_bandrange: qp_bands = [qp_bandrange(1), qp_bandrange(1)+1, ..., qp_bandrange(2)].


.. data:: macropol_calculation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "N"
   * - **Description**
     - Available options are:

       - "N" : None. Choice valid for isolated systems.
       - "C" : Include long-wavelength limit. Choice valid for condensed systems.


.. data:: n_lanczos

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 30
   * - **Description**
     - Number of Lanczos chains.


.. data:: n_imfreq

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 128
   * - **Description**
     - Number of frequecies used to sample the imaginary frequency axis in the range [0,ecut_imfreq].


.. data:: n_refreq

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 272
   * - **Description**
     - Number of frequecies used to sample the real frequency axis in the range [0,ecut_refreq].


.. data:: ecut_imfreq

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - dynamically set to the cutoff energy of the density, read from the ground state
   * - **Description**
     - Cutoff for the imaginary frequencies (in Ry).


.. data:: ecut_refreq

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 2.0
   * - **Description**
     - Cutoff for the real frequencies (in Ry).


.. data:: wfreq_eta

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.05 / 13.6056980659
   * - **Description**
     - Energy shift of the poles (in Ry).


.. data:: n_secant_maxiter

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 1
   * - **Description**
     - Maximum number of iterations in the secant solver.


.. data:: trev_secant

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.05 / 13.6056980659
   * - **Description**
     - Convergence energy threshold (in Ry) for the secant solver.


.. data:: l_enable_lanczos

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - True
   * - **Description**
     - If (False), then Lanczos solvers are turned off.


.. data:: l_enable_off_diagonal

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     -
       - If (False) then only the diagonal matrix elements of the :math:`{G_0 W_0}` self-energy are evaluated (i.e., same band).
       - If (True) then both the diagonal and off-diagonal matrix elements of the :math:`{G_0 W_0}` self-energy are evaluated (mixing different bands). In this case the upper triangular part of the self-energy matrix is calculated and written to file according to :math:`{  {\left[ \Sigma \right]}_{ij} = \frac{1}{2} \mathrm{Re} \; \left[ {\left[ \Sigma \right]}_{ij} (\epsilon^{\mathrm{QP}}_i) + {\left[ \Sigma \right]}_{ij}(\epsilon^{\mathrm{QP}}_j) \right] }`. l_enable_off_diagonal can be set to True only when the Brillouin Zone is sampled at the :math:`{\Gamma}`-point.


.. data:: l_enable_gwetot

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - Deprecated parameter.


.. data:: o_restart_time

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.0
   * - **Description**
     - Available options are:

       - If ( o_restart_time == 0 ) A checkpoint is written at every iteration of the W and G loops.
       - If ( o_restart_time >  0 ) A checkpoint is written every o_restart_time minutes in the W and G loops.
       - If ( o_restart_time <  0 ) A checkpoint is NEVER written in the W and G loops. Restart will not be possible.


.. data:: ecut_spectralf

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - 2-dim list of float
   * - **Default**
     - [-2.0,1.0]
   * - **Description**
     - Energy cutoff (in Ry) for the real frequencies. Used when wfreq_caculation contains the runlevel "P".


.. data:: n_spectralf

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 204
   * - **Description**
     - Number of frequecies used to plot the spectral function (runlevel "P"), sampling the interval [ecut_spectralf[0],ecut_spectralf[1]].

|

--------------
westpp_control
--------------


.. data:: westpp_calculation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "R"
   * - **Description**
     - Available options are:

       - "R" : Output rho, the electronic density.
       - "W" : Output the electronic wavefunctions.
       - "E" : Output the eigenpotentials.
       - "S" : Output the screened exchange constant.
       - "D" : Output the dipole matrix elements.
       - "L" : Output the localization factor and the inverse participation ratio.


.. data:: westpp_range

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - list of int
   * - **Default**
     - [1,2]
   * - **Description**
     - Range for W, E, S, and D run.


.. data:: westpp_format

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "C"
   * - **Description**
     - Available options for the output fortmat are:

       - "c" : Cube.
       - "x" : Planar average yz.
       - "y" : Planar average xz.
       - "z" : Planar average xy.
       - "s" : Spherical average.


.. data:: westpp_sign

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then the sign of the wavefunction/eigenpotential is kept in the output file.


.. data:: westpp_n_pdep_eigen_to_use

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 1
   * - **Description**
     - Number PDEP eigenpotentials to read/use.


.. data:: westpp_r0

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - 3-dim list of floats (a vector)
   * - **Default**
     - [0.0, 0.0, 0.0]
   * - **Description**
     - Position of the center (in a.u.) for spherical average plot.


.. data:: westpp_nr

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 100
   * - **Description**
     - Number of points in the spherical average plot.


.. data:: westpp_rmax

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1.0
   * - **Description**
     - Max radius (in a.u.) for the spherical average plot.


.. data:: westpp_epsinfty

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1.0
   * - **Description**
     - Macroscopic relative dielectric constant. Used in the "S" runlevel.


.. data:: westpp_box

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - 6-dim list of floats (a vector)
   * - **Default**
     - [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
   * - **Description**
     - Box [x_0, x_1, y_0, y_1, z_0, z_1] (in a.u.) within which the localization factor is computed (the "L" runlevel).

|

--------------
server_control
--------------

.. data:: document

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - jsonizable object
   * - **Default**
     - "{}"
   * - **Description**
     - The document is serialized into a JSON string and passed to the server (see `West/Pytools/west_clientserver.py`).

|


--------------
wbse_init_control
--------------


.. data:: wbse_init_calculation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "S"
   * - **Description**
     - Available options are:

       - "S" : Start from scratch
       - "R" : Restart from an interrupted run. You should restart with the same number of cores and images.

.. data:: localization

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "N"
   * - **Description**
     - Available options are:

       - "N" : Kohn-Sham orbitals are not localized
       - "B" : Bisected orbitals are used

.. data:: wfc_from_qbox

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "qb.xml"
   * - **Description**
     - Name of the Qbox sample file.

.. data:: bisection_info

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "info.bis"
   * - **Description**
     - Name of the file that contains info about bisection.

.. data:: chi_kernel

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "CHI"
   * - **Description**
     - Available options are:

       - "CHI" : W = vc + vc X vc
       - "CHI_RPA" : W = vc + vc X_RPA vc
       - "XC_CHI_RPA" : W = vc + (vc+fxc) X vc


|

--------------
wbse_control
--------------


.. data:: wbse_calculation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "D"
   * - **Description**
     - Available options are:

       - "D" : Diagonalize the Liouville super-operator with the Davidson iterative solver.
       - "d" : Restart the calculation for wbse_calculation = "D" from an interrupted run. You should restart with the same number of cores and images.
       - "L" : Compute the absorption spectrum with the Lanczos method.
       - "l" : Restart the calculation for wbse_calculation = "L" from an interrupted run. You should restart with the same number of cores and images.


.. data:: qp_correction

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "None"
   * - **Description**
     - Available options are:

       - "None" : Quasiparticle corrections are not added.
       - specify the name of the file from which quasiparticle corrections are read (in eV).


.. data:: scissor_operator

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.0
   * - **Description**
     - Value of the scissor operator (in eV).

.. data:: solver

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "BSE"
   * - **Description**
     - Available options are:

       - "BSE" : Bethe-Salpeter equation
       - "TDDFT" : Time-dependent density functional theory

.. data:: n_liouville_eigen

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 1
   * - **Description**
     - Number of Liouville eigenvectors and eigenvalues.


.. data:: n_liouville_times

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 4
   * - **Description**
     - Maximum dimension of the search space = n_liouville_eigen * n_liouville_times.


.. data:: n_liouville_maxiter

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 100
   * - **Description**
     - Maximum number of iterations of the Davidson method.


.. data:: n_liouville_read_from_file

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 0
   * - **Description**
     - Number of Liouville eigenvectors that can be read from file.


.. data:: trev_liouville

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.001
   * - **Description**
     - Absolute convergence threshold for Liouville eigenvalues.


.. data:: trev_liouville_rel

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.1
   * - **Description**
     - Relative convergence threshold for Liouville eigenvalues.


.. data:: n_lanczos

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 10
   * - **Description**
     - Number of Lanczos iterations to be performed.


.. data:: ipol_input

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "XX"
   * - **Description**
     - Determine which components of the polarizability tensor (alpha) are computed:
       - "XX": alpha_xx
       - "YY": alpha_yy
       - "ZZ": alpha_zz
        -"XYZ": three Lanczos chains are sequentially performed and the full polarizability tensor and the absorption coefficient are computed.


.. data:: epsinfty

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1.0
   * - **Description**
     - Macroscopic relative dielectric constant.


.. data:: macropol_calculation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "N"
   * - **Description**
     - Available options are:

       - "N" : None. Choice valid for isolated systems.
       - "C" : Include long-wavelength limit. Choice valid for condensed systems.

.. data:: spin_excitation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "S"
   * - **Description**
     - Available options are:

       - "S" : Singlet.
       - "T" : Triplet.

.. data:: l_preconditioning

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - Preconditioning for the Davidson methos, for finite systems.
