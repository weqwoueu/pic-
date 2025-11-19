SUBROUTINE PhiBC_2D

! Updated:	11/28/2004 03:00 AM
! Purpose:	Set potential boundary for periodic and Neumann BC
! it acts when f_periodic or f_Neumann is true
! NO effect IF they are untrue
! same version as 990214 in Optics and Plume
! dphi/dx = 0 at all outer surface

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D

IMPLICIT NONE

INTEGER		i, j

! ***********I check for f_periodic BC******************
! at the i=0, and nz+1 boundary points
! face 1,  i=0
IF(f_periodic(1)) THEN
		DO  j=0, ny+1
			phi(0,j) = phi(nx,j)
			phi(nx+1,j) = phi(1,j)
		END DO
END IF

! at the j=0, and ny+1 boundary points
! face 3, j=0
IF(f_periodic(3)) THEN
		DO i=0, nx+1
			phi(i,0)   = phi(i,ny)
            phi(i,ny+1)= phi(i,1)
		END DO
END IF


!**********II check for Neumann BC*****************
! at the i=0, and nz+1 boundary points
! face 1,  i=0
IF(f_zeroe(1)) THEN
		DO  j=0, ny+1
             phi(0,j) = phi(2,j)
		END DO
END IF

! face 2,  i=nx+1
IF(f_zeroe(2)) THEN
		DO  j=0, ny+1
             phi(nx+1,j) = phi(nx-1,j)
		END DO
END IF

! at the j=0, and ny+1 boundary points
! face 3, j=0
IF(f_zeroe(3)) THEN
		DO i=0, nx+1
			phi(i,0)   = phi(i,2)
		END DO
END IF

! face 4, j=ny+1
IF(f_zeroe(4)) THEN
		DO i=0, nx+1
            phi(i,ny+1)= phi(i,ny-1)
		END DO
END IF

!DO j=1,ny
!	DO i=1,nx
!        Aver_phi(i,j)=Aver_phi(i,j)+phi(i,j) 
!        Aver_rho(i,j)=Aver_rho(i,j)+rho(i,j) 
!        Aver_rho_s(i,j,1:2)=Aver_rho_s(i,j,1:2)+rho_s(i,j,1:2)
!    END DO
!END DO


END