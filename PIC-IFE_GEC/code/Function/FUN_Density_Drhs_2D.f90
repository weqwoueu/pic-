SUBROUTINE FUN_Density_Drhs_2D(rho_, phi_, region, f_, n)

USE IFE_MAIN_PARAM
USE IFE_RHS_PARAM
USE Object_Data_2D

IMPLICIT NONE

INTEGER									n
REAL(8), DIMENSION(n), INTENT(IN)	::	rho_, phi_
INTEGER, INTENT(IN)					::	region
REAL(8), DIMENSION(n), INTENT(OUT)	::	f_


f_ = -1 !0 ! -1.0

!INTEGER	i
! NO object BUG FIX
!INTEGER vacuum_backup

!vacuum_backup = vacuum
!vacuum =-1 

!IF (region==vacuum) THEN
!	DO i=1,SIZE(phi_)
!		IF ( phi_(i)<=phi_bkgd(-vacuum) ) THEN 
!			f_(i) = - Rho_bkgd(-vacuum) *	&
!					DEXP( (phi_(i)-phi_bkgd(-vacuum))/Te_bkgd(-vacuum) )/Te_bkgd(-vacuum)
!		ELSE
!			f_(i) = - Rho_bkgd(-vacuum) *	&
!					( One+(phi_(i)-phi_bkgd(-vacuum))/Te_bkgd(-vacuum) )/Te_bkgd(-vacuum)
!		END IF
!	END DO
!ELSE
!	DO i=1,SIZE(phi_)
!		IF ( phi_(i)<=phi_bkgd(-vacuum) ) THEN 
!			f_(i) = Zero
!		ELSE
!			f_(i) = Zero
!		END IF
!	END DO
!ENDIF

! NO object BUG FIX
!vacuum = vacuum_backup

!print*, '-vacuum=',-vacuum, ' region=',region
!print*, phi_bkgd(-vacuum),Rho_bkgd(-vacuum),Te_bkgd(-vacuum) 
!print*, 'maxval(f_)=', maxval(f_),' minval(f_)=',minval(f_)
!print*
!pause

END