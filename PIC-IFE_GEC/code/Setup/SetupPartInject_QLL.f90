SUBROUTINE SetupPartInject_QLL(Ndisf)

! Updated:	11/28/2004 03:00 AM
! Purpose:	Setup particle injection. ModIFied from SetupPartInject in THRUSTER
!			hard wired for: 
!			1) injection in z direction
!			2) only random flux in the z direction, upstream  
!			   <v> = SQRT(T/m)/SQRT(2pi) , vd is neglected

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Particle_2D
USE Field_2D
USE TimeControl
USE Object_2D

IMPLICIT NONE

INTEGER		jj, Ndisf
!DSMC
REAL(8)		D


double precision randum
external randum
INTEGER		i_part,ntotp,i
REAL(8)		spacing_x, spacing_y, x_pos, y_pos, x_vel, y_vel
INTEGER		i_xpos, i_ypos
REAL(8)		U, G  
REAL(8)		NGRN
real(8)		length_x,length_y

WRITE(6,*) 'SetupPartInject'

IF(.NOT.ALLOCATED(ipf))			ALLOCATE(ipf(ispe_tot))
IF(.NOT.ALlOCATED(N_inject_x))	ALLOCATE(N_inject_x(ispe_tot))
IF(.NOT.ALLOCATED(N_inject_y))	ALLOCATE(N_inject_y(ispe_tot))
IF(.NOT.ALLOCATED(N_inject))	ALLOCATE(N_inject(ispe_tot))
IF(.NOT.ALLOCATED(sdis))		ALLOCATE(sdis(ispe_tot))
IF(.NOT.ALLOCATED(vd))			ALLOCATE(vd(ispe_tot,3))
IF(.NOT.ALLOCATED(vt))			ALLOCATE(vt(ispe_tot,3))
IF(.NOT.ALLOCATED(affp_cell))	ALLOCATE(affp_cell(ispe_tot))
IF(.NOT. ALLOCATED(Length))		ALLOCATE(Length(ispe_tot))
IF(.NOT. ALLOCATED(R_start))	ALLOCATE(R_start(ispe_tot))
IF(.NOT. ALLOCATED(Z_start))	ALLOCATE(Z_start(ispe_tot))
! DSMC
IF(.NOT.ALLOCATED(TMPJ))		ALLOCATE(TMPJ(ispe_tot))
IF(.NOT.ALLOCATED(VMP))			ALLOCATE(VMP(ispe_tot))
IF(.NOT.ALLOCATED(SC))			ALLOCATE(SC(ispe_tot))
IF(.NOT.ALLOCATED(Flux))		ALLOCATE(Flux(ispe_tot))




!IF(.NOT.ALLOCATED(Fndj))		ALLOCATE(Fndj(ispe_tot))
!IF(.NOT.ALLOCATED(Fnum))		ALLOCATE(Fnum(ispe_tot))	
! --- READ in parameters----------------------------
! READ in the loaded speces
WRITE(6,*) 'Number of species injected =', ispe_inject
IF(ispe_inject > ispe_tot) THEN
	WRITE(6,*) 'Wrong input, ispe_tot =', ispe_tot
    STOP
END IF
IF(ispe_inject == 0) THEN
	WRITE(6,*) 'No particle injection, RETURN'
    RETURN
END IF

! Loading injection zone, number, v, for each populations
DO jj = 1, ispe_inject
! ---A info general use------------------------------
! READ particle flag 
    WRITE(6,*) 'ispe',ipf(jj)
    IF(ipf(jj) < 1 .OR. ipf(jj) > ispe_tot) THEN
		WRITE(6,*) 'Wrong input, ispe_tot =', ispe_tot
        STOP
	END IF
! READ density reference
    WRITE(6,*) 'dens0 =', dens0(ipf(jj))
    IF(dens0(ipf(jj)) < 0) THEN
		WRITE(6,*) 'Wrong input, dens0 < 0'
        STOP
	END IF

    WRITE(6,*) 'vd123', vd(ipf(jj),1),vd(ipf(jj),2), vd(ipf(jj),3)
!	WRITE(6,*) 'Length=', Length(ipf(jj))		!,Length(ipf(jj),2)
!	WRITE(6,*) 'Z_start=', Z_start(ipf(jj))		!,Z_start(ipf(jj),2)
!	WRITE(6,*) 'R_start=', R_start(ipf(jj))		!,R_start(ipf(jj),2)


! Particle temperature
	WRITE(6,*) 'TMPJ', TMPJ(ipf(jj))
! Finish READing
END DO




!Get N_inject with the DSMC knowledge
DO jj=1, ispe_inject

	N_inject(ipf(jj))=N_inject_x(ipf(jj))*N_inject_y(ipf(jj))
!	N_inject(ipf(1))=70
!	N_inject(ipf(2))=70

	VMP(ipf(jj)) = SQRT(2.*TMPJ(ipf(jj))/xm(ipf(jj)))
	IF(TMPJ(ipf(jj)) .GT. 0) THEN 
		SC(ipf(jj))  = vd(ipf(jj),2) / VMP(ipf(jj))
		IF(DABS(SC(ipf(jj))) .LT. 10.1)THEN
			CALL ERF(SC(ipf(jj)), D)
			Flux(ipf(jj)) = (DEXP(-SC(ipf(jj)) * SC(ipf(jj))) + SQRT_PI * SC(ipf(jj)) * (1. + D)) / (2. * SQRT_PI)
		ENDIF
	
		IF(SC(ipf(jj)) .GT. 10.) THEN
			Flux(ipf(jj)) = SC(ipf(jj))
		ENDIF	
	
		IF(SC(ipf(jj)) .LT. -10.) THEN
			Flux(ipf(jj)) = 0.
		ENDIF
		WRITE(6,*) 'VMP for species', ipf(jj), '=', VMP(ipf(jj))
		WRITE(6,*) 'SC for species', ipf(jj), '=', SC(ipf(jj))
		WRITE(6,*) 'N_inject for species', ipf(jj), '=', N_inject(ipf(jj))
		WRITE(6,*) 'Flux(ipf(jj)) for species', ipf(jj), '=', Flux(ipf(jj))
	ELSEIF(TMPJ(ipf(jj)) .EQ. 0) THEN
		WRITE(6,*) 'VMP for species', ipf(jj), '=', VMP(ipf(jj))
		WRITE(6,*) 'N_inject for species', ipf(jj), '=', N_inject(ipf(jj))
	ELSE
		WRITE(6,*) 'Temperature Input Error.'
		STOP
	ENDIF
		
ENDDO		
         
! Get the affp here
WRITE(6,*) 'Normalize density using beam exit density'

DO jj = 1, ispe_inject
!          affp(ipf(jj)) = sdis(ipf(jj))*hxi(3)*FLOAT(nx-1)*FLOAT(ny-1)
!     +        /FLOAT(N_inject(ipf(jj)))*dens0(ipf(jj))
!     +                    *(hxi(3)*hxi(1)*hxi(2)) 
! 2001-0915, changed defination of affp, now affp means charge per particle
! affp_cell means charge per particle per cell
	IF(TMPJ(ipf(jj)) .GT. 0) THEN
		sdis(ipf(jj)) = Flux(ipf(jj)) * VMP(ipf(jj)) * dt
	ELSEIF(TMPJ(ipf(jj)) .EQ. 0) THEN
		sdis(ipf(jj)) = abs(vd(ipf(jj), 2)) * dt
	ENDIF
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	affp(ipf(jj)) = an_gridt(1)*sdis(ipf(jj))	&
!					/FLOAT(N_inject(ipf(jj)))*dens0(ipf(jj))
    affp(ipf(jj)) = 3.62395784E-4 !!!!!!qll
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	affp_cell(ipf(jj)) = affp(ipf(jj))*(hxi(1)*hxi(2))
!	affp(ipf(jj)) = sdis(ipf(jj))*an_gridt(1)*an_gridt(2)	&
!					/FLOAT(N_inject(ipf(jj)))*dens0(ipf(jj))
	affp_cell(ipf(jj)) = affp(ipf(jj))*(hxi(1)*hxi(2)*hzi)
!	affp_cell(ipf(jj)) = affp(ipf(jj))*(hxi(1)*hxi(2)*hxi(3))

	WRITE(6,*) 'affp for species',ipf(jj),'=', &
				affp(ipf(jj)),'affp_cell =',affp_cell(ipf(jj))
END DO


!!!!!!!!!!!!!!!!!!!!!!!2014-12-14!!!!!!!!!!!!!!!!!!!!!!
!
!i = 0
!OPEN(UNIT=1, FILE='Phase_0100000.dat')
!
!READ(1,*)
!READ(1,*)
!101 i = i + 1
!READ(1,*,END=100) spacing_x, spacing_y,spacing_x, spacing_y,spacing_x, spacing_y, spacing_x
!GOTO 101
!
!100 CONTINUE
!ntot = i - 1
!PRINT*, 'Total Number of Removed Partical =',ntot		!!计算总的粒子数
!
!REWIND(1)
!
!READ(1,*)
!READ(1,*)
!DO i =1, ntot
!	READ(1,*) part(i,1), part(i,2),part(i,3),part(i,4),part(i,5),part(i,6),part(i,7)
!END DO
!
!CLOSE(1)
!
!ns=0
!
!DO i =1, ntot
!    IF(part(i,7)==1)THEN
!        ns(1)=ns(1)+1
!    ELSE
!        ns(2)=ns(2)+1
!    ENDIF
!END DO
!PRINT*, 'ns(1) =',ns(1),ns(2)
!!!!!!!!!!!!!!!!!!!!!!!2014-12-14!!!!!!!!!!!!!!!!!!!!!!



END
