MODULE Object_2D

USE FluxLine_2D
IMPLICIT NONE

! Objects (i.e. objects representing s/c bus, thruster, solar array, antenna, payloads, ... etc)
INTEGER		MaxObjects, vacuumRegion

!     +                   object_ctr(0:MaxTotObjects,3),
!     +                   object_r(0:MaxTotObjects,3),

! object geometric location and dimensions
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	f_object_l, f_object_r
INTEGER, DIMENSION(:,:), ALLOCATABLE	::	i_object_l, i_object_r

! object shape and axis
INTEGER, DIMENSION(:), ALLOCATABLE		::	shape_object, axis_object

! object size
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	dimension_object
REAL(8), DIMENSION(:), ALLOCATABLE		::	radius_object, radius2_object

! object potential and permittivity
REAL(8), DIMENSION(:), ALLOCATABLE		::	phi_object, eps_object

INTEGER, DIMENSION(:), ALLOCATABLE		::	region_object

!!!! **************** bjw add 2019.9.30 *****************
INTEGER, DIMENSION(:), ALLOCATABLE		::	direction_object
!!!! **************** bjw add 2019.9.30 *****************

! optics currnets
REAL(8), DIMENSION(:), ALLOCATABLE			::	cjcoll
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	cjcoll_s,		&
												currentLeft,	&
												currentRight,	&
												currentInLeft,	&
												currentInRight,	&
												currentHitLeft,	&
												currentHitRight,&
												currentHitSide
INTEGER, DIMENSION(:,:), ALLOCATABLE		::	ncoll_s,	&
												nLeft,		&
												nRight,		&
												nInLeft,	&
												nInRight,	&
												nHitLeft,	&
												nHitRight,	&
												nHitSide

TYPE(FluxLineType), DIMENSION(:), ALLOCATABLE	::	FluxLine

END MODULE