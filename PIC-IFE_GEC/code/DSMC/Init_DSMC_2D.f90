SUBROUTINE	Init_DSMC_2D

USE IFE_Data
USE IFE_MAIN_PARAM
USE Constant_Variable_2D
USE DSMC_Data_2D
USE Particle_2D
USE Domain_2D
IMPLICIT NONE

!double precision      randum
!external               randum

REAL(8)                   ranum
INTEGER	                nxy,i,j,inum_a,mc


nxy=nx*ny

!IF(.NOT. ALLOCATED(CS))		    ALLOCATE(CS(7,nxy))
IF(.NOT. ALLOCATED(CCG))		ALLOCATE(CCG(2,nxy))
IF(.NOT. ALLOCATED(CC))		    ALLOCATE(CC(nxy))
IF(.NOT. ALLOCATED(IR))		    ALLOCATE(IR(N_part_tot))
!CS=0
CCG=0
CC=0
IR=0

    SP(1)=2.18E-25 !SP(1) is the molecular mass

    SP(2)=4.939E-10 !SP(2) is the molecular diameter

    SPM(2)=116000                !单位：开 k
    SPM(3)=0.5
    SPM(4)=1
    SPM(1)=PI*SP(2)**2
    CALL GAM(2.5-SPM(3),SPM(5))
!SPM(5)=GAM(2.5-SPM(3))
!SPM(1) 碰撞截面
!SPM(2) the reference temperature 
!SPM(3) the viscosity temperature power law 
!SPM(4) the reciprocal of the VSS scattering parameter !VSS散射参数的倒数
!SPM(5) the gamma function of (5/2-viscosity-temperature power law)
	DO i=1,nx-1
		DO j=1,ny-1
			mc=(i-1)*(ny-1)+j
			!CC(mc)=Cell_volume(i,j)*LMDD**3           !real volume of the cell
            CC(mc)=hx(1)*hx(2)*L_ref**2           !real volume of the cell ： bjw
        	!CALL random_number(ranum)
            CALL DRandom(ranum)
			CCG(2,mc)=ranum
			!CCG(1,mc)=SPM(1)*SPM(2)*SQRT(11600*tmpj(3)*T_Ref/SPM(2))  !!! tmpj(3)--表示原子温度
            CCG(1,mc)=SPM(1)*SPM(2)*SQRT(11600*tmpj(2)*T_Ref/SPM(2))  !!! bjw：无原子，temp0(3)没用
!the maximum value of the (rel. speed)*(cross-section) is set to a reasonable, 
!but low, initial value and will be increased as necessary
		END DO
	END DO



END SUBROUTINE
