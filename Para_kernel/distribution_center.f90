!
! Copyright (C) 2015-2016 M. Govoni 
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
MODULE distribution_center
  !-----------------------------------------------------------------------
  !
  USE class_idistribute, ONLY : idistribute
  !
  IMPLICIT NONE
  !
  SAVE
  !
  TYPE(idistribute) :: pert
  TYPE(idistribute) :: macropert
  TYPE(idistribute) :: ifr
  TYPE(idistribute) :: rfr
  TYPE(idistribute) :: aband
  !
END MODULE
