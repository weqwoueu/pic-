SUBROUTINE RhoBC_2D(delta)

! Updated:	4/27/2005 01:30 PM
! Purpose:	Set charge density boundary for periodic, absorptive and reflective BC

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D
USE Particle_2D

IMPLICIT NONE

INTEGER		i, j, delta

!WRITE(6,*) 'RhoBC'

! Check for periodic particle BC
! at the i=0, and nz+1 boundary points
! face 1,  i=0
IF(f_periodic(1)) THEN
		DO  j=0, ny+1
			rho(0,   j)		= rho(nx,j)
			rho(nx+1,j)		= rho(1,j)

			rho_s(0,   j,1:ispe_tot)	= rho_s(nx,j,1:ispe_tot)
			rho_s(nx+1,j,1:ispe_tot)	= rho_s(1, j,1:ispe_tot)
		END DO
END IF

! at the j=0, and ny+1 boundary points
! face 3, j=0
IF(f_periodic(3)) THEN
		DO i=0, nx+1
			rho(i,0   )		= rho(i,ny)
            rho(i,ny+1)		= rho(i,1)

			rho_s(i,0,   1:ispe_tot)	= rho_s(i,ny,1:ispe_tot)
            rho_s(i,ny+1,1:ispe_tot)	= rho_s(i,1, 1:ispe_tot)
		END DO
END IF


! Check for reflective particle BC
! face 1,  i = 1
IF(preflect(1)) THEN
		DO  j=1, ny+1
			rho(1,j) = Two*rho(1,j)

			rho_s(1,j,1:ispe_tot) = Two*rho_s(1,j,1:ispe_tot)
		END DO
END IF

! face 2,  i = nx
IF(preflect(2)) THEN
		DO  j=1, ny+1
             rho(nx,j) = Two*rho(nx,j)

             rho_s(nx,j,1:ispe_tot) = Two*rho_s(nx,j,1:ispe_tot)
		END DO
END IF

! face 3, j = 1
IF(preflect(3)) THEN
	IF(delta == 0) THEN
		DO i=1, nx+1
			rho(i,1)   = Two*rho(i,1)

			rho_s(i,1,1:ispe_tot)   = Two*rho_s(i,1,1:ispe_tot)
		END DO
	ENDIF
END IF

! face 4, j = ny
IF(preflect(4)) THEN
		DO i=1, nx+1
            rho(i,ny)= Two*rho(i,ny)

            rho_s(i,ny,1:ispe_tot)= Two*rho_s(i,ny,1:ispe_tot)
		END DO
END IF


END
