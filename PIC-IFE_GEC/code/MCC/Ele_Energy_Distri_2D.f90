SUBROUTINE Ele_Energy_Distri_2D

!! Jinwei Bai 
!! Purpose:		electrons energy distribution. 
!! Last Update:	2019-3-10   

USE MCC_Data_2D 
USE MCC_Main_Param	
USE Constant_Variable_2D
IMPLICIT NONE
REAL(8)  :: pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10
REAL(8)  :: a1,b1,c1,a2,b2,c2,a3,b3,c3
REAL(8)  :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10
INTEGER  :: N_ela, N_exc, N_ion

REAL(8)          :: V
INTEGER          :: i
REAL(8), EXTERNAL   :: SigmaElastic,SigmaExc,SigmaIz

WRITE(*,*) 'Ele_Energy_Distri_2D'

SV_MAX_E=0

DO i = 1, energy_max   !!! bjw:energy_max???
	V=sqrt(2*E*i/Me)
	SELS(i) = SigmaElastic (REAL(i))*sqrt(1/mfactor)
	SEXC(i) = SigmaExc (REAL(i))*sqrt(1/mfactor)
	SION(i) = SigmaIz (REAL(i))*sqrt(1/mfactor)
	SV_MAX_E = MAX(SV_MAX_E, SELS(i)*V + SEXC(i)*V + SION(i)*V) !碰撞截面和速度乘积的最大值，用于处理空碰撞
END DO

OPEN(1,ACTION = 'WRITE', FILE ='Cross_section_fitting.dat' )
DO i = 1, energy_max
    WRITE(1,*) i, SELS(i), SEXC(i), SION(i)
END DO
CLOSE(1)

END SUBROUTINE

REAL(8) FUNCTION SigmaElastic ( energy )
    USE MCC_Data_2D 
    REAL(8)  :: pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10
    REAL(8)  :: a1,b1,c1,a2,b2,c2,a3,b3,c3
    REAL(8)  :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10
    INTEGER  :: N_ela, N_exc, N_ion
    
	! XSection for Elastic Collisions !!弹性碰撞截面e-Xe
	REAL energy,roote
	roote=sqrt(energy)
    
    !!!!!*********************** e-Xe **************************
	if (energy<=.1592) then
		sigmaElastic=1.699e-19
	else if (energy<=2.8) then
		sigmaElastic=1.0e-17 *(0.07588072747894*energy*energy-	&
			0.34475940259139*energy*rootE+0.58473840309059*energy- &
						0.42726069455393*rootE+0.11430271021684)
	else if (energy<=24.7) then
		sigmaElastic =1.0e-17 *(-0.00199145459640*energy*energy+&
				0.02974653588357*energy*rootE-0.16550787909579*energy+ &
						0.40171310068942*rootE-0.31727871240879)
	else if (energy<=50) then
		sigmaElastic =1.0e-17 *(-0.00217736834537*energy*rootE+	&
			0.04302155076778*energy-0.28567311384223*rootE+0.65180228051047)
	else if (energy<=500) then
		sigmaElastic =1.0e-18 *(-0.00002249610521*energy*rootE+ &
				0.00109930275788*energy-0.02071463195923*rootE+ &
					 0.22876772390428)
	else
		sigmaElastic=6.400000000000000e-20
    end if
    !!!!!*********************** e-Xe **************************
    
    !!!!!*********************** e-He **************************
    !N_ela = SIZE(section_ela,2)
    !DO i = 1, N_ela
    !    IF (energy<=section_ela(1,1)) THEN
    !        sigmaElastic=section_ela(2,1)
    !    ELSEIF ( (energy>section_ela(1,i)) .AND. (energy<=section_ela(1,i+1)) )THEN
    !        sigmaElastic= 0.5*(section_ela(2,i) + section_ela(2,i+1))
    !    ELSEIF(energy>section_ela(1,N_ela)) THEN
    !        sigmaElastic=section_ela(2,N_ela)
    !    END IF
    !END DO
    
    !a1 =       2.182  
    !b1 =       2.647  
    !c1 =       4.125 
    !a2 =       2.169 
    !b2 =       7.976 
    !c2 =        10.9  
    !a3 =  3.114D8  
    !b3 =       -1431  
    !c3 =       333.7  
    !
    !p1 =  3.998e-037 
    !p2 = -1.386e-033  
    !p3 =  1.921e-030 
    !p4 = -1.359e-027 
    !p5 =  5.186e-025 
    !p6 =  -1.03e-022 
    !p7 =  8.926e-021 
    !
    !IF (energy<=100) THEN
    !    sigmaElastic = 1.0D-20 *( a1*exp(-((energy-b1)/c1)**2) + a2*exp(-((energy-b2)/c2)**2) + &
    !                              a3*exp(-((energy-b3)/c3)**2))  
    !ELSE
    !    sigmaElastic = 1.0D0 *(p1*energy**6 + p2*energy**5 + p3*energy**4 + p4*energy**3 + &
    !                           p5*energy**2 + p6*energy + p7)
    !END IF
    !!!!!*********************** e-He **************************
    
END FUNCTION SigmaElastic

REAL(8) FUNCTION SigmaExc (energy)
    USE MCC_Data_2D 
    REAL(8)  :: pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10
    REAL(8)  :: a1,b1,c1,a2,b2,c2,a3,b3,c3
    REAL(8)  :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10
    INTEGER  :: N_ela, N_exc, N_ion
    
	!XSection for Excitation Collisions !!激发碰撞截面
	real energy,rootE
	rootE=sqrt(energy)
    
    !!!!!*********************** e-Xe **************************
	if (energy<=8.32) then
		sigmaExc=0
	else if (energy<=11) then
		sigmaExc =1.0e-16 *(0.00194724369808*energy*energy-	&
			0.02261576374741*energy*rootE+0.09807793114366*energy- &
					0.18808539260191*rootE+0.13446494003922)
	else if (energy<=25) then
		sigmaExc =1.0e-17 *(0.00069390658261*energy*energy-	&
			0.01241570210985*energy*rootE+0.08109737428153*energy- &
					0.22730324307635*rootE+0.23122639784590)
	else if (energy<=500) then
		sigmaExc =1.0e-18 *(0.00000121267639*energy*energy-	&
			0.00008169557347*energy*rootE+0.00207211887803*energy- &
					0.02409700583197*rootE+0.11701534311188)
	else 
		sigmaExc=3.950000000000000e-21
    end if
    !!!!!*********************** e-Xe **************************
    
    !!!!!*********************** e-He **************************
    !N_exc = SIZE(section_exc,2)
    !DO i = 1, N_exc
    !    IF (energy<=section_exc(1,1)) THEN
    !        sigmaExc=section_exc(2,1)
    !    ELSEIF ( (energy>section_exc(1,i)) .AND. (energy<=section_exc(1,i+1)) )THEN
    !        sigmaExc= 0.5*(section_exc(2,i) + section_exc(2,i+1))
    !    ELSEIF(energy>section_exc(1,N_exc)) THEN
    !        sigmaExc=section_exc(2,N_exc)
    !    END IF
    !END DO
    
    !p1 = -2.877e-032 
    !p2 =  1.048e-030  
    !p3 =  8.592e-028  
    !p4 = -1.298e-025  
    !p5 =  8.231e-024  
    !p6 = -2.711e-022 
    !p7 =   4.62e-021 
    !p8 = -3.184e-020 
    !
    !pp1 =  1.785e-045  
    !pp2 = -8.688e-042  
    !pp3 =  1.802e-038  
    !pp4 = -2.079e-035  
    !pp5 =  1.463e-032  
    !pp6 = -6.496e-030 
    !pp7 =  1.815e-027 
    !pp8 = -3.048e-025 
    !pp9 =  2.468e-023 
    !pp10 =  1.048e-021 
    !
    !IF (energy<=60) THEN
    !    sigmaExc = 1.0D0 *(p1*energy**7 + p2*energy**6 + p3*energy**5 + p4*energy**4 + p5*energy**3 + &
    !                       p6*energy**2 + p7*energy + p8)  
    !ELSE
    !    sigmaExc = 1.0D0 *(pp1*energy**9 + pp2*energy**8 + pp3*energy**7 + pp4*energy**6 + pp5*energy**5 + &
    !                       pp6*energy**4 + pp7*energy**3 + pp8*energy**2 + pp9*energy + pp10)
    !END IF
    !
    !IF (sigmaExc < 0.0) THEN
    !    sigmaExc = 1.0D-25
    !END IF
    !!!!!*********************** e-He **************************
      
END FUNCTION SigmaExc

REAL(8) FUNCTION SigmaIz (energy)
    USE MCC_Data_2D 
    REAL(8)  :: pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10
    REAL(8)  :: a1,b1,c1,a2,b2,c2,a3,b3,c3
    REAL(8)  :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10
    INTEGER  :: N_ela, N_exc, N_ion
    
	!XSection for Ionization Collisions  !!!电离碰撞截面
	REAL energy,rootE
	rootE=sqrt(energy)
    
    !!!!!*********************** e-Xe **************************
	if (energy<=12.1) then
		sigmaIz=0
	else if (energy<=20) then
		sigmaIz =1.0e-17 *(0.00135612832973*energy*energy-&
				0.02258559839486*rootE*energy+ 0.14035004086532	&
				*energy-0.38335664819867*rootE+0.38736677629904)
	else if (energy<=44) then
		sigmaIz=1.0e-18 *(-0.00061869954583*energy*energy+	&
				0.01448501832638*energy*rootE-0.13321973517308	&
				*energy+0.57375481836921*rootE-0.92720818547058)
	else if (energy<=360) then
		sigmaIz =1.0e-19 *(-0.00001627288393*energy*energy+	&
				0.00103294012446*energy*rootE-0.02400846159171*	&
					energy+0.21746827014037*rootE-0.18814292010734)
	else
		sigmaIz=2.440000000000000e-20	
    end if
    !!!!!*********************** e-Xe **************************
    
    
    !!!!!*********************** e-He **************************
    !N_ion = SIZE(section_ion,2)
    !DO i = 1, N_ion
    !    IF (energy<=section_ion(1,1)) THEN
    !        sigmaIz=section_ion(2,1)
    !    ELSEIF ( (energy>section_ion(1,i)) .AND. (energy<=section_ion(1,i+1)) )THEN
    !        sigmaIz= 0.5*(section_ion(2,i) + section_ion(2,i+1))
    !    ELSEIF(energy>section_ion(1,N_ion)) THEN
    !        sigmaIz=section_ion(2,N_ion)
    !    END IF
    !END DO
    
    !p1 =  2.773e-032  
    !p2 = -1.219e-029 
    !p3 =  2.078e-027 
    !p4 =  -1.66e-025 
    !p5 =  5.142e-024 
    !p6 =  7.837e-023 
    !p7 = -3.234e-021  
    !
    !pp1 =  2.145e-045
    !pp2 = -1.136e-041
    !pp3 =  2.601e-038 
    !pp4 =  -3.37e-035 
    !pp5 =  2.716e-032
    !pp6 = -1.407e-029
    !pp7 =  4.667e-027
    !pp8 = -9.494e-025
    !pp9 =  1.016e-022 
    !pp10 = -6.045e-022
    !
    !IF (energy<=120) THEN
    !    sigmaIz = 1.0D0 *(p1*energy**6 + p2*energy**5 + p3*energy**4 + p4*energy**3 + &
    !                      p5*energy**2 + p6*energy + p7)  
    !ELSE
    !    sigmaIz = 1.0D0 *(pp1*energy**9 + pp2*energy**8 + pp3*energy**7 + pp4*energy**6 + &
    !                      pp5*energy**5 + pp6*energy**4 + pp7*energy**3 + pp8*energy**2 + pp9*energy + pp10)
    !END IF
    !IF (sigmaIz < 0.0) THEN
    !    sigmaIz = 1.0D-25
    !END IF
    !!!!!*********************** e-He **************************
     
END FUNCTION SigmaIz