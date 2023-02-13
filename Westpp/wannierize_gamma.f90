!-----------------------------------------------------------------------
SUBROUTINE wannierize_gamma(m,psi,phi)
  !-----------------------------------------------------------------------
  !
  USE kinds,                  ONLY : DP
  USE fft_base,               ONLY : dfftp,dffts
  USE fft_at_gamma,           ONLY : SINGLEBAND_INVFFT,SINGLEBAND_FWFFT,DOUBLEBAND_INVFFT
  USE gvect,                  ONLY : mill
  USE mp_global,              ONLY : me_bgrp,root_bgrp,intra_bgrp_comm
  USE mp,                     ONLY : mp_sum,mp_bcast
  USE constants,              ONLY : pi, tpi
  USE cell_base,              ONLY : at,alat,bg
  USE westcom,                ONLY : wanu,wanc,wantot,proj,wann_nearestatom
  USE pwcom,                  ONLY : npw,npwx
  USE wavefunctions_module,   ONLY : psic
  USE io_push,                ONLY : io_push_title
  USE bar,                    ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: m
  COMPLEX(DP),INTENT(IN) :: psi(npwx,m)
  COMPLEX(DP),INTENT(OUT) :: phi(npwx,m)
  !
  ! Workspace
  !
  INTEGER :: i,j,l,k,irr,il
  !
  REAL(DP) :: melement(6)
  REAL(DP),ALLOCATABLE :: amatrix(:,:,:)
  COMPLEX(DP) :: umatrix_complex(m,m)
  REAL(DP) :: c,s
  REAL(DP) :: sigma,sigma_old
  REAL(DP) :: diff, d2, tempi, tempj
  REAL(DP) :: alp,bet,thetone,theta
  REAL(DP) :: wan_center(3)
  REAL(DP) :: wan_center2(3)
  INTEGER :: iter
  LOGICAL :: conv
  REAL(DP) :: dist
  CHARACTER(LEN=3) :: catm
  !
  TYPE(bar_type) :: barra
  !
  CALL io_push_title("Wann...")
  !
  wantot = m
  IF(ALLOCATED(wanu)) DEALLOCATE(wanu)
  ALLOCATE(wanu(m,m))
  IF(ALLOCATED(wanc)) DEALLOCATE(wanc)
  ALLOCATE(wanc(3,m))
  IF(ALLOCATED(amatrix)) DEALLOCATE(amatrix)
  ALLOCATE(amatrix(m,m,6))
  IF(ALLOCATED(proj)) DEALLOCATE(proj)
  ALLOCATE(proj(dffts%nnr,6))
  CALL define_proj( )
  !
  CALL start_bar_type( barra, 'wann', INT(REAL(m*(m+1))/2.))
  !
  amatrix = 0._DP
  DO i = 1, m
     DO j = i, m
        !
        CALL DOUBLEBAND_invfft(npw,psi(1,i),psi(1,j),npwx,psic,dffts%nnr)
        !
        melement = 0._DP
        DO il = 1, 6
           DO irr = 1,dffts%nnr
              melement(il) = melement(il) + REAL( psic(irr), DP) * DIMAG( psic(irr) ) * proj( irr, il)
           ENDDO
        ENDDO
        !
        amatrix(i,j,1:6) = melement(1:6)
        IF( i /= j ) amatrix(j,i,1:6) = melement(1:6)
        !
        CALL update_bar_type( barra,'wann', 1 )
     ENDDO
  ENDDO
  CALL stop_bar_type( barra, 'wann' )
  !
  CALL mp_sum(amatrix,intra_bgrp_comm)
  !
  DEALLOCATE(proj)
  !
  CALL joint_d( m, amatrix, 6, wanu )
  !
  WRITE(stdout,'(/,7X,"[WAN] Wannier centers")')
  !
  DO i=1,m
     !
     wan_center(1)=DIMAG( CDLOG( CMPLX( amatrix(i,i,1) , amatrix(i,i,2), KIND=DP )) )/tpi
     wan_center(2)=DIMAG( CDLOG( CMPLX( amatrix(i,i,3) , amatrix(i,i,4), KIND=DP )) )/tpi
     wan_center(3)=DIMAG( CDLOG( CMPLX( amatrix(i,i,5) , amatrix(i,i,6), KIND=DP )) )/tpi
     CALL fold_to_cube( wan_center, wan_center2)
     !
     wanc(:,i) = wan_center2(1)*at(:,1) + wan_center2(2)*at(:,2) + wan_center2(3)*at(:,3)
     IF( wann_nearestatom ) THEN
        CALL find_closest_atom( wanc(:,i), dist, catm)
        WRITE(stdout,'(7X,"[WAN] ",i14," : ",3f14.6,a5,f14.6)') i, wanc(:,i)*alat, catm, dist*alat
     ELSE
        WRITE(stdout,'(7X,"[WAN] ",i14," : ",3f14.6)') i, wanc(:,i)*alat
     ENDIF
     !
  ENDDO
  !
  DO j=1,m
     DO i=1,m
        umatrix_complex(i,j) = CMPLX( wanu(i,j), 0._DP)
     ENDDO
  ENDDO
  !
  CALL ZGEMM('N','N',npw,m,m,(1._DP,0._DP),psi,npwx,umatrix_complex,m,(0._DP,0._DP),phi,npwx)
  !
  DEALLOCATE(amatrix)
  !
END SUBROUTINE
!
SUBROUTINE fold_to_cube( rin, rout)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP),INTENT(IN) :: rin(3)
  REAL(DP),INTENT(OUT) :: rout(3)
  INTEGER :: i
  !
  DO i = 1, 3
     rout(i) = MODULO( rin(i), 1._DP)
  ENDDO
  !
END SUBROUTINE
!
SUBROUTINE find_closest_atom( rin, d, label)
  !
  USE ions_base,              ONLY : nat, tau, ityp, atm, nsp
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP),INTENT(IN) :: rin(3)
  REAL(DP),INTENT(OUT) :: d
  CHARACTER(LEN=3),INTENT(OUT) :: label
  !
  INTEGER :: ia
  REAL(DP) :: dnow, cx, cy, cz, x(3)
  !
  d = 2._DP
  label = "xx"
  !
  DO ia = 1, nat
     CALL fold_to_cube( tau(:,ia), x )
     cx = ABS( rin(1) - x(1) )
     cy = ABS( rin(2) - x(2) )
     cz = ABS( rin(3) - x(3) )
     IF( cx > 0.5_DP ) cx = 1._DP - cx
     IF( cy > 0.5_DP ) cy = 1._DP - cy
     IF( cz > 0.5_DP ) cz = 1._DP - cz
     dnow = SQRT( cx**2 + cy**2 + cz**2 )
     IF( dnow < d ) THEN
        d = dnow
        label = atm(ityp(ia))
     ENDIF
  ENDDO
  !
END SUBROUTINE
