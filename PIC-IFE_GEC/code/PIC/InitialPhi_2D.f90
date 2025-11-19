SUBROUTINE InitialPhi_2D

! Updated:	11/28/2004 10:45 PM
! Purpose:	Initialize potential array.

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D
USE TimeControl

!USE Object

IMPLICIT NONE

INTEGER		i, j
REAL(8)		X, Y
	
! Initialize Phi
WRITE(6,*)
WRITE(6,*) 'InitialPhi'

DO i=1,nx
	DO j=1,ny
			X = VertX(1,i,j)
			Y = VertX(2,i,j)
			IF ( i_grid_flag(i,j)==1 ) THEN
!!				IF (rho(i,j,k) < 1.E-7*den0_ref) THEN
!					Phi(i,j) = X + Y
!!				ELSE
!!					Phi(i,j) = phi0_ref !+ Te_ref*LOG( rho(i,j)/den0_ref )
!!				END IF
			END IF
	END DO
END DO


END
