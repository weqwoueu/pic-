MODULE Particle_2D

IMPLICIT NONE

! Particle Stuff
! 1 Particle Array
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	part
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	part_e  !$ ab.ZWZ
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	part_i  !$ ab,ZWZ

! 1.1 Particle Populations
INTEGER(4)									ispe_tot, ispe_load, ispe_inject, ispe_downstm

! 2 Physical Constant
REAL(8), DIMENSION(:), ALLOCATABLE		::	qs, qm, xm

! Particle Charge, q,  this is a variable
REAL(8)										q

! 3 Number and Velocity for Each Species
! 3.0 Number of Loaded, and Injected Particles, and Maximum Particles
!INTEGER(8)									ntot, N_part_tot
INTEGER(8)									 N_part_tot
REAL(8)									ntot
!INTEGER(8), DIMENSION(:), ALLOCATABLE	::	ns
Real(8), DIMENSION(:), ALLOCATABLE	::	ns

! 3.1 Number of Loaded Particles
INTEGER(4), DIMENSION(:), ALLOCATABLE	::	N_part_x,  N_part_y, N_part

REAL(8), DIMENSION(:,:), ALLOCATABLE	::	vd, vt

! 3.2 Collective Momentum and Energy for Each Species
REAL(8), DIMENSION(:), ALLOCATABLE		::	px, py, ek
!REAL(8), DIMENSION(:,:), ALLOCATABLE		::	Ek_one

!---below are problem specific----
! 4 Particle Loading Zone 
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	ploadlength, ploadorigin

! 5 Particles Injection Source
INTEGER(4), DIMENSION(:), ALLOCATABLE	::	ipf, N_inject_x, N_inject_y, N_inject
REAL(8), DIMENSION(:), ALLOCATABLE		::	sdis, Surface, Length, Radius_inject, X_inject, Y_inject
REAL(8), DIMENSION(:), ALLOCATABLE		::	Z_start, R_start

! 6 Adjust Outer
INTEGER(4), DIMENSION(:,:), ALLOCATABLE	::	ntrade_s, nrefl_s, nloss_s
!INTEGER(4), DIMENSION(:), ALLOCATABLE	::	newloss
! 7 Output Removed Particles Option
LOGICAL						OutputRemovedParticles

LOGICAL						OutputGlobalMoment

! DSMC
REAL(8), DIMENSION(:), ALLOCATABLE		::	TMPJ, VMP, SC, Flux

REAL(8), DIMENSION(:,:), ALLOCATABLE	::	colli_point
REAL(8), DIMENSION(:), ALLOCATABLE		::	surface_jump_Q

REAL(8), DIMENSION(:), ALLOCATABLE      ::  t_parameter
REAL(8)                        Te_Sec

!!!! 统计壁面电荷丢失
INTEGER  ::  Nloss,Neloss,Niloss
INTEGER  ::  Nloss_local1,Neloss_local1,Niloss_local1  !!!内壁
INTEGER  ::  Nloss_local2,Neloss_local2,Niloss_local2  !!!外壁

REAL(8), DIMENSION(:), ALLOCATABLE :: YINI    !$ ab.ZWZ
REAL(8), DIMENSION(:), ALLOCATABLE :: Rweight    !$ ab.ZWZ
REAL(8), DIMENSION(2)              :: dy_inject !$ ab.ZWZ for weighting

END MODULE
