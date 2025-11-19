SUBROUTINE MCC_2D(it)
    
!! Jinwei Bai 
!! Purpose:		simulate the process of MCC. 
!! Last Update:	2019-3-10  

USE MCC_Data_2D
USE MCC_Main_Param
USE IFE_Data

USE Particle_2D
USE Field_2D
USE Domain_2D

USE Constant_Variable_2D
USE PIC_MAIN_PARAM_2D
USE TimeControl
IMPLICIT NONE

! niucollide is the collission frequency,njubohm is the bohm frequency
! Niucollide为碰撞频率，Njubohm为波姆频率
! double precision          randum
! external                   randum

INTEGER		it
REAL(8)    :: RANUM,RANUM1,RANUM2
REAL(8)	   :: nju_bohm,njut,nju1,nju2,nju3
REAL(8)	   :: gden_max,nju_collide_max,nju_bohm_max,nju_max,p1,p2,p3,p4,para

INTEGER    :: i,j,k,ns1 
REAL(8)	   :: ENERGY,P, V
INTEGER	   :: i_part,isp,num_old,ionized_particle     !num_old:total amount of electron before ionization
CHARACTER*40	:: fname, filename
!INTEGER, DIMENSION(2,(nx-1)*(ny-1)*(nz-1))		::	IC		!LC 2014
INTEGER, DIMENSION(2,(nx-1)*(ny-1))		::	IC		!LC 2014
!REAL(8), DIMENSION(:,:,:), ALLOCATABLE      ::  VertX  ! bjw 2019-3-11
REAL(8)             :: nsold, bmax

WRITE(6,*)
WRITE(6,*) 'MCC_2D '
WRITE(6,*) 'Before MCC_2D:'
WRITE(6,*) 'ntot         =', ntot
WRITE(6,*) 'ns           =', (ns(j),j=1,ispe_tot)

ns=0   !!!!! qll

DO i_part=1,ntot
    IF(part(i_part,7)== 1)THEN
        ns(1)=ns(1)+1
    ELSEIF(part(i_part,7)== 2)THEN
        ns(2)=ns(2)+1
    ELSEIF(part(i_part,7)== 3)THEN
        ns(3)=ns(3)+1       
    END IF
END DO
ntot=0

DO i=1, ispe_tot
    ntot = ntot+ns(i)
ENDDO

 
WRITE(6,*)
WRITE(6,*) 'MCC_2D '
WRITE(6,*) 'Before MCC_2D:'  
WRITE(6,*) 'ntot         =', ntot
WRITE(6,*) 'ns           =', (ns(j),j=1,ispe_tot)

nsold = 0

!	if(ispe_tot>=3) then
!	    gden(:,:,:)=rho_s(:,:,:,3) !如果模拟了原子的运动，则将原子密度给gden，否则原子密度用表达式
!	end if


	gden_max=MAXVAL(gden)
    !bmax = MAXVAL(bfz) 
    bmax = MAXVAL(bt)  !$ mb.ZWZ 2021/9/11
	!write(*,*) gden_max,bmax

	nju_collide_max=gden_max*NERO *  SV_MAX_E !最大电子和原子碰撞频率1E18
    
	!nju_bohm_max=co_bohm1*E*bmax/Me*SQRT(1/mfactor) !最大波姆碰撞频率*B_Ref
	nju_bohm_max=co_bohm1*E*Bfield_ref/Me*SQRT(1/mfactor) !最大波姆碰撞频率*B_Ref
    
	nju_max=nju_collide_max+nju_bohm_max
!	P=nju_max*dt/omigp	 !电子的最大碰撞概率
	!P=nju_max*dt*LMDD/Vel_Ref	 !电子的最大碰撞概率
    
    P=nju_max*dt*L_ref/v_ref	 !电子的最大碰撞概率 bjw (take taylor expand's first two term)
    
    !P = 1. - EXP(-nju_max*dt*L_ref/v_ref) !$ the maximum collision probability of electron ( expotential form )
    
	ionized_particle=0 !每次计算mcc之前电离产生的等离子首先清零
    !print *,SV_MAX_E,nju_collide_max,nju_bohm_max,nju_max,P
   
    !IF(it>atomstep)THEN
        IF(ipf(ispe_inject)>=3) CALL INDEXM(3, IC, nx-1, ny-1) !确定原子在每个网格内的分布，以便于考虑原子的电离过程	!LC 2014
    !ENDIF
	num_old=ntot

	DO i_part = 1,num_old    
	!!!!!DO i_part = ns(1)+ns(2),num_old  !!!!qll 
	    isp = INT(part(i_part,7))

	    IF(isp==1) THEN
		   
            CALL DRandom(RANUM)

		    IF (RANUM < P) THEN    !小于电子的最大碰撞概率
            !用于计算各种碰撞率


               !if(free_path_diag==1 .AND. j<=free_path_num) call ave_free_path(insp,i_part)

                ENERGY = t_parameter(isp)*(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)
            
		        CALL Collision_Frequency_2D(ENERGY,i_part,nju1,nju2,nju3,njut,nju_bohm)
		        !print *,ENERGY,nju1,nju2,nju3,njut,nju_bohm
		        !STOP
		        njut=nju1+nju2+nju3+nju_bohm
		
                CALL DRandom(RANUM1)
            
     !           IF (RANUM1 <= (nju1+nju2+nju_bohm)/nju_max) THEN
     !               energy=1
    !				v=SQRT(energy/t_parameter(isp))
    !                CALL Anisotropicv(v,part(i_part,4),part(i_part,5),part(i_part,6))	  	 
			    IF (RANUM1 <= nju_bohm / nju_max) THEN  !计算波姆碰撞率
                    V = SQRT(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)
				    CALL Anisotropicv(v,part(i_part,4),part(i_part,5),part(i_part,6))
                    CALL PTOG(bohm,i_part)
    !				ENER(i_part)=t_parameter(isp)*(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)			
    !                call oned_av(insp,inum_a,i_part,onedf,6,imax1,2)
                    CALL Anomalous_Difusion(i_part)

			    ELSE IF(RANUM1>nju_bohm/nju_max .AND. RANUM1<=(nju1+nju_bohm)/nju_max) then 		!计算弹性碰撞率	!ELASTIC COLLISION
                    V = SQRT(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)
                    CALL Anisotropicv(v,part(i_part,4),part(i_part,5),part(i_part,6))
                    CALL PTOG(elastic,i_part)
    !				ENER(i_part)=t_parameter(isp)*(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)
    !                call oned_av(insp,inum_a,i_part,onedf,6,imax1,3)
                    CALL Anomalous_Difusion(i_part)
                
			    ELSE IF (RANUM1 > (nju1+nju_bohm)/nju_max .AND. RANUM1 <= (nju1+nju2+nju_bohm)/nju_max) THEN !		!计算激发碰撞率	
				    ENERGY=ENERGY-8.32
				    v=SQRT(ENERGY/t_parameter(isp))
                    CALL Anisotropicv(v,part(i_part,4),part(i_part,5),part(i_part,6))
                    CALL PTOG(excite,i_part)
    !                call oned_av(insp,inum_a,i_part,onedf,6,imax1,4)			
				    ENER(i_part)=ENERGY !t_parameter(isp)*(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)
                
 		        ELSE IF (RANUM1 > (nju1+nju2+nju_bohm)/nju_max .AND. RANUM1<njut/nju_max) THEN !!计算电离碰撞率IONIZATION COLLISION
                    !WRITE(*,*) '11111111111111111',ENERGY
        !            ENERGY=ENERGY-12.13
        !            CALL DRandom(RANUM2)
	       !         ENERGY=ENERGY-abs(8.7*tan(RANUM2*atan(ENERGY/(2*8.7))))
				    !v=SQRT(ENERGY/t_parameter(isp))
        !            CALL Anisotropicv(v,part(i_part,4),part(i_part,5),part(i_part,6))
				    CALL Ionization_2D(i_part,ENERGY,ionized_particle, IC )              
                    CALL PTOG(ionize,i_part)	                                !$ 获取电离密度分布
    !               call oned_av(insp,inum_a,i_part,onedf,6,imax1,5)
                    ENER(i_part)=ENERGY !t_parameter(isp)*(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)	
                
	            END IF
    !			ENER(i_part)=t_parameter(isp)*(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)
		    END IF
		END IF
	END DO
!$OMP END PARALLEL DO
   !DO isp=1,ispe_tot
   !    ns(isp)=ns(isp)+ionized_particle !将新产生的电子数和离子数加入数组
   !END DO


!       j=1
! 100   continue   
!       IF(part(j,8)<0)THEN
!       CALL Remove_2D(j,ntot)
!       j=j-1
!       ENDIF
!       j=j+1
!       if(j<=ntot)goto 100

WRITE(*,*) 'ionized_particle', ionized_particle
print *,'+++++++++',P,ENERGY,nju1,nju2,nju3,njut,nju_bohm


ns=0
DO i_part=1,ntot
    IF(part(i_part,7)== 1)THEN
        ns(1)=ns(1)+1
    ELSEIF(part(i_part,7)== 2)THEN
        ns(2)=ns(2)+1
    ELSEIF(part(i_part,7)== 3)THEN
        ns(3)=ns(3)+1       
    END IF
END DO
ntot=0

DO i=1, ispe_tot
    ntot = ntot+ns(i)
ENDDO


WRITE(6,*)
WRITE(6,*) 'After MCC_2D:'
WRITE(6,*) 'ntot         =', ntot
WRITE(6,*) 'ns           =', (ns(j),j=1,ispe_tot)
WRITE(6,*) 'nsold           =',nsold

!!!----------------------------------------------------LC 2015.1.23
IF(MOD(it,2000)==0)THEN

WRITE(fname,999) it
999 FORMAT(I6.6)

filename = './OUTPUT/MCC/MCC_IJK_'//TRIM(fname)//'.dat'

PRINT*, 'Writing to file: ', TRIM(filename)

OPEN(123, ACTION = 'WRITE', FILE = TRIM(filename))
WRITE(123,*) 'TITLE = "Field Plot"'
WRITE(123,*) 'VARIABLES = "x" "y" "bohm" "elastic" "excite" "cond" "ionize"'
WRITE(123,456) nx, ny

	DO j=1,ny
		DO i=1,nx
			WRITE(123,567) VertX(1:2,i,j), bohm(i,j), elastic(i,j), &
			            excite(i,j), cond(i,j), ionize(i,j)
		END DO
	END DO

456 FORMAT (' ZONE I = ',I6,', J= ',I6)
567 FORMAT (F8.2,' ',F8.2,' ',F8.2,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
CLOSE(123)
ENDIF
!!!----------------------------------------------------

END SUBROUTINE MCC_2D