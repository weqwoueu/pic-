SUBROUTINE  Collision_Frequency_2D(ENERGY,i_part,nju1,nju2,nju3,njut,nju_bohm)

USE MCC_Data_2D
USE MCC_Main_Param
USE Domain_2D
USE Field_2D
USE Constant_Variable_2D
USE Particle_2D
USE IFE_Data
USE IFE_INTERFACE, ONLY: GTOP
IMPLICIT NONE

INTEGER		ispe, i_part,i,j,k

REAL(8)		ENERGY,SE1,SE2,SE3
INTEGER     FLAG_ENERGY
REAL(8)		nju1,nju2,nju3,njut,nju_bohm
REAL(8)		V,bxpre,bypre,bzpre,btpre,gdenpre
REAL(8), DIMENSION(:,:), ALLOCATABLE    ::	para	
!REAL(8), DIMENSION(:,:), ALLOCATABLE    ::  bfx,bfy  !!! bjw 2019

ispe = part(i_part,7)
ALLOCATE(para(0:nx+1,0:ny+1))
DO j = 0, ny+1
	DO i = 0, nx+1
			para(i,j)=gden(i,j)*NERO            !*1E18
    END DO
END DO
      

CALL GTOP(i_part,gdenpre,para)  ! 插值得到电子位置处的原子密度


IF (ENERGY > ENERGY_MAX) ENERGY = ENERGY_MAX -1
FLAG_ENERGY =ANINT( ENERGY)
IF (FLAG_ENERGY==0) FLAG_ENERGY=1	
	
SE1 = SELS(FLAG_ENERGY)
SE2 = SEXC(FLAG_ENERGY)
SE3 = SION(FLAG_ENERGY)

!V=SQRT(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)*Vel_Ref

V=SQRT(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)*v_Ref
nju1=gdenpre* V * Se1
nju2=gdenpre* V * Se2
nju3=gdenpre* V * Se3

!ALLOCATE(bfx(0:nx+1,0:ny+1),bfy(0:nx+1,0:ny+1))  !! bjw 2019
!bfx = Bfield(1,:,:)
!bfy = Bfield(2,:,:)

CALL GTOP(i_part,bxpre,bfx)
CALL GTOP(i_part,bypre,bfy)
CALL GTOP(i_part,bzpre,bfz) !$ mb.ZWZ 2021/9/11
!btpre=B_Ref*SQRT(bxpre**2+bypre**2)       !4**bmax
!btpre=Bfield_ref*SQRT(bxpre**2+bypre**2)       !4**bmax
btpre=Bfield_ref*SQRT(bxpre**2+bypre**2+bzpre**2)       !$ mb.ZWZ 2021/9/11
!    if(ipret<=50) then
    nju_bohm=co_bohm*E*btpre/Me*SQRT(1/mfactor)
!    else
!       nju_bohm=co_bohm1*e*btpre/m(ispe)*sqrt(1/mfactor)
!    end if
!if(gdenpre<2.e18/sfactor) nju_bohm=nju_bohm*64
!if(ipret>=60) nju_bohm=nju_bohm*4
njut=nju1+nju2+nju3+nju_bohm
!WRITE(*,*) 'Collision_Frequency_2D:',nju1,nju2,nju3,nju_bohm
!STOP

DEALLOCATE(para)    !$ ab.ZWZ 2021/9/11

END 
