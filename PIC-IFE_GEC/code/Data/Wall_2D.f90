MODULE Wall_2D

IMPLICIT NONE

! wall geometric location and dimensions	
INTEGER(4), DIMENSION(:), ALLOCATABLE		::	nnode,npo	

REAL(8), DIMENSION(:), ALLOCATABLE			::	alpha
										
! inject angle
REAL(8), DIMENSION(:), POINTER			::	theta,diffY,vx,vy
REAL(8), DIMENSION(:), POINTER			::	all_angle,included_vx,included_vy
! Sputtering energy
REAL(8), DIMENSION(:), ALLOCATABLE			::	k,Y
REAL(8), DIMENSION(:), ALLOCATABLE			::	sputwall,aera
REAL(8)											ionbeita,wallk

! potential
REAL(8), DIMENSION(:), ALLOCATABLE				::	E,Eth

!REAL(8), DIMENSION(:,:), POINTER	::  alphax,Ex,Ethx,kx,xpo
REAL(8), DIMENSION(:,:,:), POINTER	::  alphax,E_par
INTEGER(4), DIMENSION(:,:,:), POINTER			::	nc_par	


INTEGER, DIMENSION(:), ALLOCATABLE		::  n_tot


END MODULE