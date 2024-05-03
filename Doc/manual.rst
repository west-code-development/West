.. _Manual:

Manual
======

The complete **WEST** reference for input parameters.

.. contents:: :local:
              :depth: 1

.. note::
   Not all input parameters listed below are mandatory. Check :ref:`quickreference` or :ref:`tutorial` pages to see examples of input files.

.. seealso::
   **WESTpy** is a Python package, designed to assist users of the WEST code in pre- and post-processing massively parallel calculations. Click `here <https://west-code.org/doc/westpy/latest/>`_ to know more.

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
       - "R" : Restart from an interrupted run. You should restart with the same number of MPI tasks and parallelization flags.
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

       - If ( n_steps_write_restart >  0 ) A checkpoint is written every n_steps_write_restart iterations in the Davidson loop.
       - If ( n_steps_write_restart <= 0 ) A checkpoint is NEVER written in the Davidson loop. Restart will not be possible.

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
     - If (True), then the exact-exchange term in the Hamiltonian is computed with the cutoff of the wavefunction. Used only when n_exx_lowrank = 0.

.. data:: n_exx_lowrank

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - dynamically set to match the number of bands, read from the ground state
   * - **Description**
     - If ( n_exx_lowrank > 0 ), then the exact-exchange is computed with a low-rank approximation of rank n_exx_lowrank.

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
     - List of bands to compute the QP corrections. If nspin = 2, two distinct sets of bands can be specified for the two spin channels, for example, [[1, 3, 4], [2, 3, 4]]. If nspin = 2 and qp_bands specifies only one set of bands, then the same bands are used for both channels. If qp_bands is not set, it is determined from qp_bandrange as [qp_bandrange(1), qp_bandrange(1)+1, ..., qp_bandrange(2)].

.. data:: macropol_calculation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "C"
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
     - 21
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

.. data:: l_qdet_verbose

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - Controls what two-body terms of the QDET effective Hamiltonian are written to file.

       - If (False), then only the partially screened two-body terms are written to file.
       - If (True), then the fully screened, partially screened, and bare two-body terms are written to file.

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
       - If (False), then only the diagonal matrix elements of the :math:`{G_0 W_0}` self-energy are evaluated (i.e., same band).
       - If (True), then both the diagonal and off-diagonal matrix elements of the :math:`{G_0 W_0}` self-energy are evaluated (mixing different bands). In this case the upper triangular part of the self-energy matrix is calculated and written to file according to :math:`{ {\left[ \Sigma \right]}_{ij} = \frac{1}{2} \mathrm{Re} \; \left[ {\left[ \Sigma \right]}_{ij} (\epsilon^{\mathrm{QP}}_i) + {\left[ \Sigma \right]}_{ij}(\epsilon^{\mathrm{QP}}_j) \right] }`. l_enable_off_diagonal can be set to True only when the Brillouin Zone is sampled at the :math:`{\Gamma}`-point.

.. data:: n_pdep_eigen_off_diagonal

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 0
   * - **Description**
     - If ( n_pdep_eigen_off_diagonal > 0 ), then the off-diagonal matrix elements of the :math:`{G_0 W_0}` self-energy are computed using n_pdep_eigen_off_diagonal PDEPs. This is to reduce file system usage in large-scale QDET calculations. The diagonal matrix elements are always computed using n_pdep_eigen_to_use PDEPs. Used only when l_enable_off_diagonal is True.

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
     - Number of frequecies used to plot the spectral function (the "P" runlevel), sampling the interval [ecut_spectralf[0],ecut_spectralf[1]].

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
       - "X" : Output the exciton state.
       - "P" : Output the density response to exciton state.
       - "B" : Output the unitary transformation matrix of Boys/Wannier localization.
       - "C" : Output the decomposition of BSE/TDDFT excited state, and the transition dipole moments.
       - "M" : Output the spin multiplicity of the BSE/TDDFT excited state (<S^2>, nspin = 2).

.. data:: westpp_range

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - list of int
   * - **Default**
     - [1,2]
   * - **Description**
     - Range for W, E, S, D, L, X, and P run.

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

       - "C" : Cube.
       - "X" : Planar average yz.
       - "Y" : Planar average xz.
       - "Z" : Planar average xy.
       - "S" : Spherical average.

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
     - Number of PDEP eigenpotentials to read/use.

.. data:: westpp_r0

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - 3-dim list of floats (a vector)
   * - **Default**
     - [0.0, 0.0, 0.0]
   * - **Description**
     - Position of the center (in a.u.) for spherical average plot or localization factor in a sphere.

.. data:: westpp_nr

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 100
   * - **Description**
     - Number of points in the spherical average plot or localization factor in a sphere.

.. data:: westpp_rmax

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1.0
   * - **Description**
     - Max radius (in a.u.) for the spherical average plot or localization factor in a sphere.

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

.. data:: westpp_l_spin_flip

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then a spin-flip calculation is performed. Used only when westpp_calculation = "C" or "M" and nspin = 2.

.. data:: westpp_wannier_tr_rel

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1e-6
   * - **Description**
     - Relative convergence threshold for Wannier localization. Used only when westpp_calculation = "B".

.. data:: westpp_l_compute_tdm

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then the transition dipole moment is computed. Used only when westpp_calculation = "C".

.. data:: westpp_l_dipole_realspace

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - Controls how the dipole is computed. Used only when westpp_calculation = "C".

       - If (False), then the dipole is computed in the reciprocal space by computing [H,r]. Choice valid for isolated and condensed systems.
       - If (True), then the dipole is computed in the real space. Choice valid for isolated systems only.

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


-----------------
wbse_init_control
-----------------

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

       - "S" : Start from scratch.
       - "R" : Restart from an interrupted run. You should restart with the same number of MPI tasks and parallelization flags.

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

       - "BSE" : Bethe-Salpeter equation.
       - "TDDFT" : Time-dependent density-functional theory.

.. data:: bse_method

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "PDEP"
   * - **Description**
     - Available options are:

       - "PDEP" : Use the PDEP eigenpotentials to compute screened exchange integrals.
       - "FF_QBOX" : Use the finite field method with Qbox coupling to compute screened exchange integrals.

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

       - "N" : Kohn-Sham orbitals are not localized.
       - "B" : Bisected orbitals are used. Valid only when bse_method = "FF_QBOX".
       - "W" : Wannier orbitals are used.

.. data:: wannier_tr_rel

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1e-6
   * - **Description**
     - Relative convergence threshold for Wannier localization. Used only when localization = "W".

.. data:: wfc_from_qbox

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "qb_wfc"
   * - **Description**
     - Name of the file that contains Qbox wavefunctions. Used only when bse_method = "FF_QBOX" and localization = "B".

.. data:: bisection_info

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "bis_info"
   * - **Description**
     - Name of the file that contains info about bisection. Used only when bse_method = "FF_QBOX" and localization = "B".

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

       - "CHI" : :math:`{W = v_c + v_c \chi v_c}`
       - "XC_CHI" : :math:`{W = v_c + (v_c+f_{xc}) \chi v_c}`

       :math:`{W}` and :math:`{v_c}` are the screened and bare Coulomb interactions, respectively, :math:`{\chi}` is the density-density response function, :math:`{f_{xc}}` is the exchange-correlation potential.

       In addition to :math:`{\chi}`, :math:`{\chi_{\mathrm{RPA}}}` or :math:`{\chi_{\mathrm{IPA}}}` may be requested by specifying "approximation: RPA" or "approximation: IPA" in the document keyword of the server_control section (see also `West/Pytools/west_clientserver.py`).

.. data:: overlap_thr

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.0
   * - **Description**
     - If the overlap between two orbitals is below this threshold, the corresponding screened exchange integral is not computed. Used only when localization = "B" or "W".

.. data:: n_trunc_bands

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 0
   * - **Description**
     - If n_trunc_bands > 0, then the n_trunc_bands lowest occupied bands are not considered when summing over occupied bands.

|


------------
wbse_control
------------

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
       - "d" : Restart the calculation for wbse_calculation = "D" from an interrupted run. You should restart with the same number of MPI tasks and parallelization flags.
       - "L" : Compute the absorption spectrum with the Lanczos method.
       - "l" : Restart the calculation for wbse_calculation = "L" from an interrupted run. You should restart with the same number of MPI tasks and parallelization flags.

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
       - Specify the name of the Wfreq output file (in JSON format) from which quasiparticle corrections are read.

.. data:: scissor_ope

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.0
   * - **Description**
     - Value of the scissor operator (in Ry).

.. data:: n_liouville_eigen

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 1
   * - **Description**
     - Number of Liouville eigenvectors and eigenvalues. Used only when wbse_calculation = "D" or "d".

.. data:: n_liouville_times

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 4
   * - **Description**
     - Maximum dimension of the search space = n_liouville_eigen * n_liouville_times. Used only when wbse_calculation = "D" or "d".

.. data:: n_liouville_maxiter

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 100
   * - **Description**
     - Maximum number of iterations of the Davidson method. Used only when wbse_calculation = "D" or "d".

.. data:: n_liouville_read_from_file

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 0
   * - **Description**
     - Number of Liouville eigenvectors that can be read from file. Used only when wbse_calculation = "D" or "d".

.. data:: trev_liouville

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.001
   * - **Description**
     - Absolute convergence threshold for Liouville eigenvalues. Used only when wbse_calculation = "D" or "d".

.. data:: trev_liouville_rel

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 0.1
   * - **Description**
     - Relative convergence threshold for Liouville eigenvalues. Used only when wbse_calculation = "D" or "d".

.. data:: n_lanczos

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 1000
   * - **Description**
     - Number of Lanczos iterations to be performed. Used only when wbse_calculation = "L" or "l".

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

       - If ( n_steps_write_restart >  0 ) A checkpoint is written every n_steps_write_restart iterations in the Davidson or Lanczos loop.
       - If ( n_steps_write_restart <= 0 ) A checkpoint is NEVER written in the Davidson or Lanczos loop. Restart will not be possible.

.. data:: wbse_ipol

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "XX"
   * - **Description**
     - Controls which components of the polarizability tensor (alpha) are computed:

       - "XX": alpha_xx.
       - "YY": alpha_yy.
       - "ZZ": alpha_zz.
       - "XYZ": three Lanczos chains are sequentially performed and the full polarizability tensor and the absorption coefficient are computed.

.. data:: l_dipole_realspace

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - Controls how the dipole is computed. Used only when wbse_calculation = "L" or "l".

       - If (False), then the dipole is computed in the reciprocal space by computing [H,r]. Choice valid for isolated and condensed systems.
       - If (True), then the dipole is computed in the real space. Choice valid for isolated systems only.

.. data:: wbse_epsinfty

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1.0
   * - **Description**
     - Macroscopic relative dielectric constant.

.. data:: spin_excitation

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - string
   * - **Default**
     - "S"
   * - **Description**
     - Used only when nspin = 1. If nspin = 2, use instead l_spin_flip to choose spin-conserving or spin-flip excitations. Available options are:

       - "S" : Singlet.
       - "T" : Triplet. Valid only when solver = "BSE".

.. data:: l_preconditioning

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then preconditioning is used. Used only when wbse_calculation = "D" or "d". Should be set to True in most cases.

.. data:: l_pre_shift

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then the preconditioner is shifted by the corresponding Kohn-Sham orbital energy. Used only when wbse_calculation = "D" or "d" and l_preconditioning is True. Should be set to True for isolated systems and False for perodic systems.

.. data:: l_spin_flip

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then a spin-flip calculation is performed. Used only when wbse_calculation = "D" or "d" and nspin = 2.

.. data:: l_spin_flip_kernel

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then the spin-flip kernel is used in the spin-flip calculations. Used only in spin-flip TDDFT calculations. Should be set to False for spin-flip BSE calculations.

.. data:: l_spin_flip_alda0

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then the ALDA0 approximation is used in the spin-flip kernel, i.e. the gradient correction to the exchange-correlation potential is discarded in the spin-flip kernel. Used only in spin-flip TDDFT calculations using GGA type exchange-correlation functionals and when l_spin_flip_kernel is True.

.. data:: l_print_spin_flip_kernel

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then the spin-flip kernel is written to a cube file. Used only in spin-flip TDDFT calculations and when l_spin_flip_kernel is True.

.. data:: spin_flip_cut

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1e-3
   * - **Description**
     - Spin-flip cutoff to prevent divergence by setting values to zero on a grid point if the density on this point is smaller than spin_flip_cut. Used only in spin-flip TDDFT calculations.

.. data:: l_forces

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then analytical forces are computed for the forces_state excited state. Used only when wbse_calculation = “D” or “d”.

.. data:: forces_state

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 1
   * - **Description**
     - Excited state for which analytical forces are computed. Used only when l_forces is True.

.. data:: forces_zeq_cg_tr

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1e-10
   * - **Description**
     - Convergence threshold in the CG method that solves the Z vector equation. Used only when l_forces is True.

.. data:: forces_zeq_n_cg_maxiter

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 500
   * - **Description**
     - Maximum number of iterations in the CG solver of the Z vector equation. Used only when l_forces is True.

.. data:: ddvxc_fd_coeff

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1e-2
   * - **Description**
     - Finite difference step size to compute the second derivative of the local part of the exchange-correlation kernel. Used only when l_forces is True.

.. data:: forces_inexact_krylov

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - 0
   * - **Description**
     - Apply the inexact krylov subspace approximation in the CG solver of the Z vector equation. Used only when l_forces is True.

       - 0 : Do not apply the approximation.
       - 1 : Skip the exact-exchange term once the CG solver converges to forces_inexact_krylov_tr.
       - 2 : Skip the K_1d term once the CG solver converges to forces_inexact_krylov_tr.
       - 3 : Skip the K_2d term once the CG solver converges to forces_inexact_krylov_tr.
       - 4 : Skip the K_1d and K_2d terms once the CG solver converges to forces_inexact_krylov_tr.
       - 5 : Skip the exact-exchange, K_1d, and K_2d terms once the CG solver converges to forces_inexact_krylov_tr.

.. data:: forces_inexact_krylov_tr

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - float
   * - **Default**
     - 1e-16
   * - **Description**
     - Apply the inexact krylov subspace approximation if the norm of the residual vector is smaller than this threshold. Used only when l_forces is True.

.. data:: l_minimize_exx_if_active

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - False
   * - **Description**
     - If (True), then the exact-exchange term in the Hamiltonian is computed with the cutoff of the wavefunction. Used only when n_exx_lowrank = 0.

.. data:: n_exx_lowrank

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - int
   * - **Default**
     - dynamically set to match the number of bands, read from the ground state
   * - **Description**
     - If ( n_exx_lowrank > 0 ), then the exact-exchange is computed with a low-rank approximation of rank n_exx_lowrank.

.. data:: l_reduce_io

.. list-table::
   :widths: 10 90
   :stub-columns: 0

   * - **Type**
     - bool
   * - **Default**
     - True
   * - **Description**
     - Deprecated parameter.
