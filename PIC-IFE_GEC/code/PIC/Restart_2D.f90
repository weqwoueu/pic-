SUBROUTINE Restart_2D

! Updated:	11/26/2004 05:04 AM
! Purpose:	Restart run

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D
USE Particle_2D
USE TimeControl
USE IFE_Data
USE Wall_2D
! include 'neutral-cex.inc'

IMPLICIT NONE

INTEGER	        lap, lapdummy, i, j, jj, size_part_2,itmp, n_sect
REAL(8)        domain_x_min, domain_x_max, domain_y_min, domain_y_max
        
CHARACTER*30 fname1, fname2, fname3

IF (irestart==1)	THEN
	lap = ilap
!	load dumped field and particles

!	Use this subroutine to restart at a known lap
!	CALL ReStart(lap) 

	WRITE(*,*) 'Restart_2D'

! the result's necessary for restart is
! a) control variables
! things like the domain size, neutral cloud position, etc, 
! may be reconstructe
! b) results: 1) e, b field. 2)particles
!  
!	READ in control variables
	WRITE(fname1, 100) lap
100 FORMAT('./DUMP/var',I7.7,'dump')

	OPEN (1, ACTION = 'READ', FILE = fname1)
        
	    READ(1,*) lapdummy  
	    READ(1,*) ntot
	    READ(1,*) (ns(jj),jj=1,ispe_tot)
      Read(1,*) nx, ny
      !> ab.ZWZ 2022/3/28 for RF condition
      Read(1,*) Period_res
      Read(1,*) iStep_previous
      Read(1,*) nPeriod

	CLOSE(1)

!------0. output some control variables
	WRITE(6,*) 'control variables for Restart at it = ',lap
	WRITE(6,*) 'nlast, ntot & ns, and mx_read,my_read,mz_read:'
	WRITE(6,*)  lap
	WRITE(6,*)  ntot
	WRITE(6,*)  (ns(jj),jj=1,ispe_tot)
	WRITE(6,*)  nx,ny

!------1. READ dumped phi data----------------------------------------
	WRITE(fname1, 101) lap
101 FORMAT('./DUMP/phi',I7.7,'dump')

	OPEN(UNIT=78, FORM='UNFORMATTED', FILE=fname1, STATUS='UNKNOWN')

	!READ(78) ((phi(i,j),j=0,ny+1),i=0,nx+1)
  Read(78) (phi(i,1),i=1,Field_Size)  !LY modification for Multi-Layer-Mesh, 2022-7-25
	CLOSE(UNIT=78)

!------2. READ dumped par data----------------------------------------
	WRITE(fname2, 102) lap
102 FORMAT('./DUMP/par',I7.7,'dump')
	size_part_2 = SIZE(part,2)
	OPEN(UNIT=79, FORM='UNFORMATTED', FILE=fname2, STATUS='UNKNOWN')
!print*,'size_part_2 =',size_part_2 
	READ(79) ((part(i,j),j=1,size_part_2),i=1,ntot)
!	READ(79) (part(i,1:size_part_2),i=1,ntot)
	CLOSE(UNIT=79)

!------3. READ dumped surface charge data----------------------------------------
!	WRITE(fname3, 103) lap
!103 FORMAT('charge',I7.7,'dump')
!	OPEN(UNIT=80, FORM='UNFORMATTED', FILE=fname3, STATUS='UNKNOWN')
!!print*,'size_part_2 =',size_part_2 
!!	READ(80) ((collectq(i,j),j=0,ny-1),i=0,nx-1)
!    size_part_2 = SIZE(collectq,2)
!    READ(80) ((collectq(i,j),j=1,size_part_2),i=1,6)
!!	READ(79) (part(i,1:size_part_2),i=1,ntot)
!	CLOSE(UNIT=80)

!    WRITE(fname3, 104) lap
!104 FORMAT('Sput_IJ_',I6.6,'.dat')    
!    OPEN(123, ACTION = 'READ', FILE = fname3)
!    READ(123,*)
!    READ(123,*)
!    n_sect=SIZE(alphax(:,1,1))
!		DO i=1,n_sect
!			READ(123,*) itmp,alphax(i,1,1), alphax(i,1,2), E_par(i,1,1),    &
!			        E_par(i,1,2), nc_par(i,1,1), nc_par(i,1,2)
!		END DO
!    !601 FORMAT (I6,E15.6,E15.6,E15.6,E15.6,I10,I10)
!    CLOSE(123)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!temp code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!domain_x_min = 0
!!domain_x_max = 100
!!domain_y_min = 19.9999
!!domain_y_max = 100
!!
!!jj = 1
!!
!!DO WHILE (jj <= ntot)
!!    IF(part(jj,1) .LE. domain_x_min .OR. part(jj,1) .GT. domain_x_max .OR.  &
!!        part(jj,2) .LE. domain_y_min .OR. part(jj,2) .GT. domain_y_max) THEN
!!        j=part(jj, 7)
!!        DO i=1,SIZE(part,2)
!!	        part(jj, i) = part(ntot, i)
!!		END DO
!!		ns(j) = ns(j) - 1
!!		ntot=ntot-1
!!		jj = jj - 1
!!    ENDIF
!!    jj = jj + 1
!!ENDDO
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!temp code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END IF

END


