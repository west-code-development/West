!
! Copyright (C) 2015-2017 M. Govoni 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file: 
! Matteo Gerosa
!
!-----------------------------------------------------------------------
MODULE types_bz_grid
  !-----------------------------------------------------------------------
  !
  USE class_bz_grid,   ONLY : bz_grid
  !
  IMPLICIT NONE
  !
  SAVE
  !
  TYPE(bz_grid) :: k_grid
  TYPE(bz_grid) :: q_grid
  TYPE(bz_grid) :: kmq_grid
  TYPE(bz_grid) :: kpq_grid
  !
END MODULE
