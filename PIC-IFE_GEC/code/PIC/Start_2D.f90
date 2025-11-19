SUBROUTINE Start_2D

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name:Start_2D.f90                                        C
!
!  Purpose: PIC input data initialize
!                                                                      C
!  Author: Yuchuan Chu                                Date: 04-May-12  C
!  Comments: modify because normlization                               C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  
!==========LY modification, 2022-7-25==========
Use IFE_Data, only : Field_Size
!==========LY modification, 2022-7-25==========
USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D
USE Particle_2D
USE TimeControl
USE IMPIC_Data_2D !!! bjw add for impic 2019-6-3
IMPLICIT NONE

CHARACTER*80	dummyline
INTEGER			i, j, k, i_part, ispe, jj

IF (irestart==0)	ilap = 0

!	: Total number of simulation particles
WRITE(6,*) 'Total Number of Simulation Particles, N_part_tot = ', N_part_tot
IF(N_part_tot > N_part_max) THEN
   WRITE(6,*) 'Wrong Input!, N_part_max = ',N_part_max
   STOP
ENDIF
!	Output Removed Particles Option
IF(OutputRemovedParticles) WRITE(6,*)'OutputRemovedParticles'
IF(OutputGlobalMoment) WRITE(6,*)'OutputGlobalMoment'


!	ispe_tot: Total number of species
WRITE(6,*) 'Total Number of Populations, ispe_tot = ', ispe_tot
IF(ispe_tot > ispe_max) THEN
   WRITE(6,*) 'Wrong Input!, ispe_max = ',ispe_max
   STOP
ENDIF

! 1 time control
! nt : Number of simulation time steps
! dt : Time step
WRITE(6,*) 'nt, dt = ', nt, dt
! subcyling
! n_updatefld	: Field update each n_updatefld steps
WRITE(6,*) 'n_updatefld = ', n_updatefld

! 2  diagnostics and output
! n_pdiag	: Particle diagnosis each n_pdiag steps
WRITE(6,*) 'n_pdiag    = ', n_pdiag
! n_stride	: Output information of each n_stride particles
WRITE(6,*) 'n_stride   = ', (n_stride(jj), jj=1,ispe_tot)
! n_engydiag: Energy diagnosis each n_engydiag steps
WRITE(6,*) 'n_engydiag = ', n_engydiag
! n_gdiag	: Grid diagnosis each n_gdiag steps
WRITE(6,*) 'n_gdiag = ', n_gdiag
! n_dump	: Dump simulation restart data each n_dump steps
WRITE(6,*) 'n_dump  = ', n_dump

WRITE(6,*) 'index_xyz  = ', index_x, index_y, index_z
WRITE(6,*) 'xyz_center = ', xcenter, ycenter, zcenter


! 3 READ in basic physical parameters

! 3.1  the speed of light no use for ES code
! c	: Simualtion speed of light
WRITE(6,*) 'Simulation Speed of Light, c = ', c

! 3.2  particle quantaties,
! qs	: Charge
! xm	: Mass
! qm	: Charge to Mass ratio
DO jj=1, ispe_tot
	WRITE(6,*) 'Species ', jj
	WRITE(6,*) 'qs = ', qs(jj),'xm = ', xm(jj)
	WRITE(6,*) 'qm = ', qm(jj)
END DO

! 4. default affp and reference density
! affp	: Ratio of 'simulation' desnity to 'real' density
! dens0	: Initial density
WRITE(6,*) 'default affp and dens0'
DO jj=1, ispe_tot
	WRITE(6,*) 'Species ', jj
	WRITE(6,*) 'affp = ', affp(jj),' dens0 = ',dens0(jj)
END DO

! Initial 0 particles in the system
! ntot: Total number of simulation particles
ntot = 0
!ns	: Number of simulation particles of a certain population
ALLOCATE(ns(ispe_tot))
DO jj = 1, ispe_tot
	ns(jj) = 0
END DO

!!! ************************ bjw add for impic 2019-6-3 **********************************************
!! Initialize the particle array
!ALLOCATE(part(N_part_tot,7))
!DO i_part = 1, N_part_tot
!	part(i_part,:) = Zero
!END DO
	
!!! ************************ bjw add for impic 2019-6-3 **********************************************
! Initialize the particle array
!ALLOCATE(part(N_part_tot,10)) !!! 123-xyz,456-vxyz,7-index,8910-axyz£¬ 11-12MCC
!$ ========== ab.ZWZ ======== \\
!IF (delta_global == 1) THEN
!    !ALLOCATE(YINI(N_part_tot))  
!    ALLOCATE(Rweight(N_part_tot))
!    Rweight = 1.
!ENDIF
!$ ========== ab.ZWZ ======== //

!DO i_part = 1, N_part_tot
!	part(i_part,:) = Zero
!END DO

!=========LY modification, 2022-5-19=========
Allocate(Chi(Field_Size, 1))
Chi = 0.0
Allocate(OneAndChi(Field_Size, 1))
OneAndChi = 1.0
Allocate(TransChi(9, Field_Size, 1))
TransChi = 0.0

!Note:The below code used only test non-isotropic finite element convergence.
!TransChi(1,:,:) = 10.0
!TransChi(5,:,:) = 10.0
!TransChi(2,:,:) = 1.0
!TransChi(4,:,:) = 1.0
!=========LY modification, 2022-5-19=========


!ALLOCATE(Chi(0:nx+1,0:ny+1))
!Chi = 0.0
!!ALLOCATE(OneAndChi(0:nx+1,0:ny+1))
!!OneAndChi = 1.0
!ALLOCATE(TransChi(9,0:nx+1,0:ny+1))
!TransChi = 0.0
!TransChi(1,:,:) = 1.0
!TransChi(5,:,:) = 1.0
!TransChi(1,:,:) = 1.0   !$ ab.ZWZ 
!TransChi(2,:,:) = 2.0      
!TransChi(4,:,:) = 2.0   
!TransChi(5,:,:) = 4.0  
!Bfiled_index = .true.

!!! BJW ADD for implicit PIC
!ALLOCATE( A_bar(3,N_part_tot),A_bar_n_minus_1(3,N_part_tot), A_bar_n(3,N_part_tot), A_n_plus_1(3,N_part_tot)) 
!DO i_part = 1, N_part_tot
!	
!    A_bar_n_minus_1(:,i_part) = 0.0
!    A_bar_n(:,i_part) = 0.0
!    A_bar(:,i_part) = 0.0
!    A_n_plus_1(:,i_part) = 0.0
!    
!END DO

!!! ************************ bjw add for impic 2019-6-3 **********************************************

END
