SUBROUTINE SetPhiOuterBC_2D

! Updated:	11/26/2004 04:52 AM
! Purpose:	Set outer boundary phi boundary conditions

! 2003-0326: modified from the original
! now: for Dirichlet boundary conditions,
! i hard wire the outer 2 layers i=0 AND 1, AND i=mx AND mx+1
! to be the given potential. this way, i can always have the boundary
! set from i=1 to mx 
! 990228, this is usefule only when i need to specify phi
! at ny outer boundary

!******************************************************************
USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D

IMPLICIT NONE

INTEGER		i, j, iface

WRITE(6,*)
WRITE(6,*) 'SetPhiOuterBC'

WRITE(6,*) 'phiouter =', (phiouter(iface),iface = 1,4)

! set i = 0 face 
!	IF(outerface(1)) THEN
IF(.NOT.f_zeroe(1).AND..NOT.f_periodic(1)) THEN
		DO  j = 1, ny
			phi(0,j) = phiouter(1)
			i_grid_flag(0,j) = 1
!        solidgrid(0,j,k)=.true.
! 2nd layer for domain BC set at i=1
			phi(1,j) = phiouter(1)
			i_grid_flag(1,j) = 1
		END DO
END IF

! set i = nx+1 face 
!	IF(outerface(2)) THEN
IF(.NOT.f_zeroe(2).AND..NOT.f_periodic(2)) THEN
		DO  j=1, ny
			phi(nx+1,j) = phiouter(2)
			i_grid_flag(nx+1,j) = 1
! 2nd layer for domain BC set at i=nx
			phi(nx,j) = phiouter(2)
			i_grid_flag(nx,j) = 1
		END DO
END IF

! set j = 0 face 
!	IF(outerface(3)) THEN
IF(.NOT.f_zeroe(3).AND..NOT.f_periodic(3)) THEN
		DO  i=0, nx+1
			phi(i,0) = phiouter(3)
			i_grid_flag(i,0) = 1
! 2nd layer for domain BC set at j=1
			phi(i,1) = phiouter(3)
	        i_grid_flag(i,1) = 1
		END DO
END IF

! set j=ny+1 face 
!	IF(outerface(4)) THEN
IF(.NOT.f_zeroe(4).AND..NOT.f_periodic(4)) THEN
		DO  i=0, nx+1
			phi(i,ny+1) = phiouter(4)
			i_grid_flag(i,ny+1) = 1
! 2nd layer for domain BC set at j=ny
			phi(i,ny) = phiouter(4)
			i_grid_flag(i,ny) = 1
		END DO
END IF

END