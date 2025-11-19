SUBROUTINE Init_MCC_2D
  
!! Jinwei Bai 
!! Purpose:		Initializ the process of MCC. 
!! Last Update:	2019-3-10     
    
USE MCC_Data_2D 
USE MCC_Main_Param
USE IFE_Data
USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Particle_2D
USE Field_2D
USE IFE_RHS_PARAM
USE TimeControl
USE Constant_Variable_2D
IMPLICIT NONE


!double precision randum
!external randum

REAL(8)		ranum

INTEGER	   ispe,i,j,i_part,atomNTOT,ioS
INTEGER	   num_ea(2)
REAL(8)	    RO,R1, L,Z1,ytemp
REAL(8)	    t_int(2)
!REAL(8), DIMENSION(:,:), ALLOCATABLE			::	xp,yp      !,phip
REAL(8)	::	xp,yp
REAL(8)	::	p1,p2,p3,p4,p5,p6
INTEGER, DIMENSION(:,:,:), ALLOCATABLE ::     x,y



IF(.NOT. ALLOCATED(bohm))		    ALLOCATE(bohm(0:nx+1,0:ny+1))
IF(.NOT. ALLOCATED(elastic))		ALLOCATE(elastic(0:nx+1,0:ny+1))
IF(.NOT. ALLOCATED(excite))		    ALLOCATE(excite(0:nx+1,0:ny+1))
IF(.NOT. ALLOCATED(ionize))		    ALLOCATE(ionize(0:nx+1,0:ny+1))
IF(.NOT. ALLOCATED(recom))		    ALLOCATE(recom(0:nx+1,0:ny+1))
IF(.NOT. ALLOCATED(cond))		    ALLOCATE(cond(0:nx+1,0:ny+1))
!IF(.NOT. ALLOCATED(capa))		    ALLOCATE(capa(0:nx+1,0:ny+1))
IF(.NOT. ALLOCATED(gden))		    ALLOCATE(gden(0:nx+1,0:ny+1))
IF(.NOT. ALLOCATED(ENER))		    ALLOCATE(ENER(N_part_tot))

bohm=0
elastic=0
excite=0
ionize=0
recom=0
cond=0
!capa=0
gden=0
ENER=0

!CALL Collision_Cross_Section_Table   !!! bjw 2019-3-30

CALL Ele_Energy_Distri_2D

!!Vel_Ref = SQRT(2*kb*T_Ref*11600/Me)           !!!!!2kT/M???
!Vel_Ref = SQRT(kb*T_Ref*11600/Me)
!ALLOCATE(t_parameter(ispe_tot))
!!ALLOCATE(t_parameter(2))
!!t_parameter = 0.5 * xm * Vel_Ref**2 / E
!t_parameter(1) = 0.5 *Me * Vel_Ref**2 / E 
!t_parameter(2) = 0.5 *MI * Vel_Ref**2 / E 
!t_parameter(3) = 0.5 *MI * Vel_Ref**2 / E
!!LMDD=1
!WRITE(*,*)'t_parameter',t_parameter,Me,MI,Vel_Ref,E
!LMDD = SQRT(EPSILON0*kb*11600*T_Ref/(E*E*NERO))
!
!IF(pfactor>0) THEN
!	LMDD=LMDD*sqrt(pfactor)
!ENDIF
!LMDD=nint(LMDD*1.0E3)/1.0E3
!omigp=SQRT(E*E*NERO/(EPSILON0*Me))          !SQRT(E*E*NERO/(EPSILON0*MI))
!!Fnum=NERO*(lmdd**3)
!
!IF(pfactor>1) THEN
!!	lmdd=lmdd*sqrt(pfactor)
!	omigp=omigp / sqrt(pfactor)
!ENDIF
!!B_Ref=Me * omigp / E
!B_Ref=Me * (Vel_Ref/LMDD) / E
!!c=1     !in MKS unit, c is useless
!bmax=0.02
!omigc=E*bmax/Me         !e*bmax/Me
!print *,'omigc=',omigc
!omigc=omigc/omigp

!if(ispe_tot<3) then
!********************load_neutral****************************************************************
!gden=50


!!!! ***************** 给出原子的分布 ***************************
!WRITE(6,*) 'Load_Neutral'
!IF (.NOT.ALLOCATED(x))                  ALLOCATE(x(nx,ny))
!IF (.NOT.ALLOCATED(y))                  ALLOCATE(y(nx,ny))
!
!!den_neu=0.0
!OPEN ( unit=76, ACTION = 'READ', FILE ='DENSITY_NEU_CELL.dat')
!READ(76,*)
!READ(76,*)
!READ(76,*)
!DO j=1,ny
!    DO i=1,nx
!	    READ(76,*)  x(i,j), y(i,j), gden(i,j)
!    END DO
!END DO
!!5556	  FORMAT(3(I5),D18.8)
!CLOSE(UNIT=76)
!gden=gden/NERO !无量纲
!DEALLOCATE(x)
!DEALLOCATE(y)
!!!! ***************** 给出原子的分布 ***************************

!!!! ***************** 原子均布 ***************************
!!!Linear model Poly5:
!!!     f(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
!!!Coefficients (with 95% confidence bounds):
!       p1 = -1.421E7  
!       p2 =  3.125E10 
!       p3 = -2.649E13 
!       p4 =  1.083E16 
!       p5 = -2.168E18 
!       p6 =  1.844E20 
!
!DO j=1,ny
!    DO i=1,nx
!	    gden(i,j) = p1*VertX(2,i,j)**5 + p2*VertX(2,i,j)**4 + &
!                    p3*VertX(2,i,j)**3 + p4*VertX(2,i,j)**2 + p5*VertX(2,i,j) + p6
!    END DO
!END DO


OPEN(270,ACTION='READ',FILE='./INPUT/gden_input.dat')
    READ(270,*)  
    DO j=1,ny
		DO i=1,nx
			READ(270,*) gden(i,j) ! rho_s(i,j,k,3)*NERO
		END DO
	END DO
CLOSE(270)

!DO j=1,ny
!    DO i=1,nx
!        !IF (VertX(1,i,j)<=32.5 ) THEN
!        !    gden(i,j) = 6.0*1.0E21
!        !ELSE
!        !    gden(i,j) = 6.0*1.0E21*EXP(-8.0*(VertX(1,i,j)/70.0-0.5)* &
!        !                             (VertX(1,i,j)/70.0-0.5))
!        !ENDIF
!        
!        IF (VertX(1,i,j)<=20.0 ) THEN
!            gden(i,j) = 1.5*1.0E20
!        ELSE
!            gden(i,j) = 1.5*1.0E20*EXP(-8.0*(VertX(1,i,j)/70.0-2./7.)* &
!                                     (VertX(1,i,j)/70.0-2./7.))
!        ENDIF
!        
!        !gden(i,j) = 3.0*1.0E20
!    END DO
!END DO


OPEN(271, ACTION='WRITE',FILE='./OUTPUT/Gden.dat')
WRITE(271,*) 'TITLE = "Field Plot"'
WRITE(271,*) 'VARIABLES = "x" "y"  "gden"'
WRITE(271,811) nx, ny
	
	DO j=1,ny
		DO i=1,nx
			WRITE(271,821) VertX(1:2,i,j), gden(i,j) ! rho_s(i,j,k,3)*NERO
		END DO
	END DO

811 FORMAT (' ZONE I = ',I6,', J= ',I6)
821 FORMAT (E15.6,' ',E15.6,' ',E15.6)
CLOSE(271)



gden=gden/NERO
!gden=2.0E20/NERO

!WRITE(*,*) MAXVAL(gden)
!PAUSE
!!!! ***************** 原子均布 ***************************

!************************************************************************************************
!gden=50
!elseif(ispe_tot>=3) then
!----------------------------------------------------------------
!IF(irestart/=1)THEN
!atomNTOT=0
!OPEN(123,ACTION='READ',FILE='atomPhase.data')
!Rewind(123)
!read(123,*)
!read(123,*)
!    do while(.true.)
!        read(123,*,ioStat = ioS)
!        if ( ioS /= 0 ) Exit
!        atomNTOT=atomNTOT+1
!    enddo
!Rewind(123)
!close(123)
!
!write(*,*) 'atomNTOT=',atomNTOT
!write(*,*) t_parameter(3)
!ntot=atomNTOT
!ns(3)=atomNTOT
!OPEN(123,ACTION='READ',FILE='atomPhase.data')
!read(123,*)
!read(123,*)
!do i=1,atomNTOT
!part(i,7)=3
!read(123,*)part(i,1),part(i,2),part(i,3),part(i,4),part(i,5),part(i,6)
!!read(123,*)part(i,1),part(i,2),part(i,9),part(i,4),part(i,5)
!part(i,4)=part(i,4)/Vel_Ref
!part(i,5)=part(i,5)/Vel_Ref
!part(i,6)=part(i,6)/Vel_Ref
!part(i,8)=5000
!part(i,9)=part(i,2)
!ENER(i)=t_parameter(3)*(part(i,4)**2+part(i,5)**2+part(i,6)**2)
!CALL PTOG(rho_s(:,:,:,3),i,0)
!enddo
!close(123)
!ENDIF
!-----------------------------------------------------------------


!pause
!endif
!DEALLOCATE(xp,yp)      !,phip

END SUBROUTINE
    