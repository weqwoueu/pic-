MODULE Domain_2D

USE Boundary_Data_2D
IMPLICIT NONE

! Domain Size and Boundaries
INTEGER											nx, ny, nz
REAL(8)											an_gridt(3),f_left_wall(3),f_right_wall(3)
REAL(8)											hx(2), hxi(2), hz, hzi
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	VertX
REAL(8)											Vert_o(3)
!REAL(8), DIMENSION(:,:), ALLOCATABLE	::	Cell_volume

!INTEGER(4)										k_up, k_down

! Grid Node Type Flag
INTEGER(1), DIMENSION(:,:), ALLOCATABLE	::	i_grid_flag

! Boundary Type
LOGICAL											outerface(6),periodic(6)

! Field
LOGICAL											f_periodic(4), f_zeroe(4)
REAL(8)											phiouter(4)

! Particles
LOGICAL											pabsorb(6),preflect(6),pemit(6)
TYPE(BoundaryType), DIMENSION(:), ALLOCATABLE	::	boundaries
REAL(8), DIMENSION(:,:,:), ALLOCATABLE          ::  collect
REAL(8), DIMENSION(:,:), ALLOCATABLE            ::  collectq

! Coordinate system ------ ab.ZWZ 2021/7/9
INTEGER                                         :: delta_global
INTEGER                                         :: region_type
Real(8)                                         :: dxmin,dxminmin, dxmax,dxmaxmax, dymin,dyminmin, dymax,dymaxmax, dzmin, dzmax
real(8)                                         :: Radius_L,Radius_U,Theta_L,Theta_U   !初始区域为扇形

! NOTE:
!
! Field BC flag
! f_periodic:	periodic for phi
! f_zeroe:		zero e field, zero phi grident
!
! Particle BC flag
! periodic:		currently used only for particles eventually will be used for both field and particles
! pabsorb:		absorb particles
! preflect:		mirror reflection
! pemit:		emit particles, right now for emit ambient


! -------- just for test----------
Integer :: Time1, Time2
!--------------------------------

END MODULE