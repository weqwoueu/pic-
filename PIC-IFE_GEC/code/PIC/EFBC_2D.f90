SUBROUTINE EFBC_2D

! Updated:	11/28/2004 11:30 PM
! Purpose:	Electric field BC
! originally, this subroutine enforces E field for Dirchelet 
! and E=0 on i=0 and nx+1 points.
! now:  modified to enforce both Dirchlet and E=0 BC at i=1 and nx points
! the original one is EfBC.f.1
! the version that enforces E BC at both i=0 and i=1 ; i=nx, nx+1 is EfBC.f.2
! this one is simplified from EfBC.f.2

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D

IMPLICIT NONE

INTEGER		i, j

! for boundary grids,
!  use eq(3-29) and (3-30) on p.44 of Anderson, Tannehill, Pletcher CFD 
! to get 2nd order accuracy


! I components that can be calculated in a regular way
! this has been taken care of in GetEfield.f
! part II changes the component that has been affed

! II special treatment for components at boundarys
!  Ex at i=1
IF(f_periodic(1)) THEN
! f_periodic boundary
		DO j=0, ny+1
			efx(0,j)=efx(nx,j)
		END DO
ELSEIF(f_zeroe(1)) THEN
! neumann boundary
		DO j=0, ny+1
			efx(1,j)=Zero
		END DO
ELSE
! fixed boundary, 1 side derivative
		DO j=0, ny+1
!          efx(0,j,k)= -1.*(-3.*phi(0,j,k)+4.*phi(1,j,k)-phi(2,j,k))
!     #                 /2.*deltxm1
			efx(1,j)= -1.*(-3.*phi(1,j)+4.*phi(2,j)-phi(3,j))/Two*hx(1)
		END DO
END IF

!  Ex at i=nx
IF(f_periodic(2)) THEN
! f_periodic boundary
		DO j=0, ny+1
			efx(nx+1,j)=efx(1,j)
		END DO
ELSEIF(f_zeroe(2)) THEN
! neumann boundary
		DO j=0, ny+1
			efx(nx  ,j)=Zero
		END DO
ELSE
! fixed boundary, 1 side derivative
		DO j=0, ny+1
!      efx(nx+1,j,k)= -1.*(3.*phi(nx+1,j,k)-4.*phi(nx,j,k)+phi(nx-1,j,k))
!     #                 /2.*deltxm1
			efx(nx,j)= -1.*(3.*phi(nx,j)-4.*phi(nx-1,j)+phi(nx-2,j))/Two*hx(1)
		END DO
END IF

!  Ey at j=1
IF(f_periodic(3)) THEN
		DO i=0, nx+1
			efy(i,0)=efy(i,ny)
		END DO
ELSEIF(f_zeroe(3)) THEN
		DO i=0, nx+1
			efy(i,1)=Zero
		END DO
ELSE
		DO i=0, nx+1
!         efy(i,0,k)= -1.*(-3.*phi(i,0,k)+4.*phi(i,1,k) -phi(i,2,k))
!     #                 /Two/hx(2)
			efy(i,1)= -1.*(-3.*phi(i,1)+4.*phi(i,2) -phi(i,3))/Two/hx(2)
		END DO
END IF

!  Ey at j=ny
IF(f_periodic(4)) THEN
		DO i=0, nx+1
             efy(i,ny+1)=efy(i,1)
		END DO
ELSEIF(f_zeroe(4)) THEN
		DO i=0, nx+1
             efy(i,ny)=Zero
		END DO
ELSE
		DO i=0, nx+1
!      efy(i,ny+1,k)= -1.*(3.*phi(i,ny+1,k)-4.*phi(i,ny,k)+phi(i,ny-1,k))
!     #                 /Two/hx(2)
			efy(i,ny)= -1.*(3.*phi(i,ny)-4.*phi(i,ny-1)+phi(i,ny-2))/Two/hx(2)

		END DO
END IF


END





