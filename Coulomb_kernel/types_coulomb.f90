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
! Marco Govoni
!
!-----------------------------------------------------------------------
MODULE types_coulomb
!-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE class_coulomb,   ONLY : coulomb 
  !
  IMPLICIT NONE
  !
  SAVE
  !
  TYPE(coulomb) :: pot3D 
  !
END MODULE
