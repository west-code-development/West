MODULE wbseppcom
  !
  USE kinds,                ONLY : dp
  !
  LOGICAL :: l_meg
  LOGICAL :: l_eig_decomp
  LOGICAL :: l_lz_spec
  LOGICAL :: l_exc_plot
  LOGICAL :: l_exc_rho_res_plot
  !
  INTEGER :: wbsepp_type
  !
  ! lzc part
  !  
  INTEGER :: itermax, itermax0, ipol, sym_op, units, verbosity
  INTEGER :: spin_channel
  ! 
  CHARACTER(len=60) :: extrapolation
  REAL(dp) :: start,end,increment 
  REAL(dp) :: epsil
  !
  ! exc plot part
  !
  REAL(dp) :: r0_input(3)
  INTEGER  :: iexc_plot
  ! 
END MODULE
