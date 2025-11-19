SUBROUTINE Dump_2D(lap)

! Updated:	11/28/2004 10:30 PM
! Purpose:	Dump the results for restart.

USE PIC_Main_Param_2D
USE Domain_2D
USE Field_2D
USE Particle_2D
Use TimeControl !> ab.ZWZ 2022/3/28

!=========LY modification for Multi-Layer-Mesh, 2022-7-25=========
Use IFE_Data
!=========LY modification for Multi-Layer-Mesh, 2022-7-25=========

IMPLICIT NONE

INTEGER			lap

CHARACTER*30	fname1, fname2
INTEGER			i, j, k, jj, size_part_2

WRITE(*,*) 'Dump at time lap', lap

! the rsult's necessary for restart is
! a) control variables
! things like the domain size, neutral cloud position, etc, 
! may be reconstructe
! b) results: 1) e, b field. 2)particles
!  
!------0. output some control variables
WRITE(fname1, 100) lap
100 FORMAT('./DUMP/var',I7.7,'dump')

OPEN (1, ACTION = 'WRITE', FILE = fname1)
    WRITE(1,*) lap        
    WRITE(1,*) ntot
    WRITE(1,*) (ns(jj),jj=1,ispe_tot)
    WRITE(1,*)  nx,ny
    !> ab.ZWZ 2022/3/28 for RF condition
    Write(1,*) Period_res 
    Write(1,*) iStep_previous 
    Write(1,*) nPeriod
CLOSE(1)

WRITE(6,*) 'control variables for Restart at it=',lap
WRITE(6,*) 'nlast, ntot & ns, and mx_read,my_read,mz_read:'
WRITE(6,*)  lap
WRITE(6,*)  ntot
WRITE(6,*)  (ns(jj),jj=1,ispe_tot)
WRITE(6,*)  nx,ny

!------1. dump phi data----------------------------------------
WRITE(fname1, 101) lap
101 FORMAT('./DUMP/phi',I7.7,'dump')

OPEN (78,FORM='UNFORMATTED',FILE=fname1, STATUS= 'UNKNOWN')
!WRITE(78)((phi(i,j),j=0,ny+1),i=0,nx+1)
WRITE(78)(phi(i,1), i=1,Field_Size)  !LY modification for Multi-Layer-Mesh, 2022-7-25
CLOSE(78)

!------2. dump par data----------------------------------------
WRITE(fname2, 102) lap
102 FORMAT('./DUMP/par',I7.7,'dump')

OPEN (79,FORM='UNFORMATTED',FILE=fname2, STATUS= 'UNKNOWN')
size_part_2 = SIZE(part,2)
!print*,'size_part_2 =',size_part_2 
WRITE(79) ((part(i,j),j=1,size_part_2),i=1,ntot)
CLOSE(79)

!------3. dump rho data----------------------------------------
!WRITE(fname1, 103) lap
!103 FORMAT('rho',I5.5,'dump')
!
!OPEN (80,FORM='UNFORMATTED',FILE=fname1, STATUS= 'UNKNOWN')
!WRITE(80)((rho(i,j),j=0,ny+1),i=0,nx+1)
!CLOSE(80)

!------4. dump surface charge data----------------------------
WRITE(fname1, 104) lap
104 FORMAT('./DUMP/charge',I7.7,'dump')

!OPEN (81,FORM='UNFORMATTED',FILE=fname1, STATUS= 'UNKNOWN')
!!WRITE(81)((collectq(i,j),j=0,ny-1),i=0,nx-1)
!size_part_2 = SIZE(collectq,2)
!WRITE(81)((collectq(i,j),j=1,size_part_2),i=1,6)
!CLOSE(81)


END