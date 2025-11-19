SUBROUTINE	Anomalous_Difusion(i_part)

USE Particle_2D
USE Constant_Variable_2D
USE MCC_Main_Param
USE MCC_Data_2D
USE Field_2D
!USE PIC_MAIN_PARAM_2D
USE IFE_Data
USE IFE_INTERFACE, ONLY: GTOP    
IMPLICIT NONE

!    double precision           randum
!    external                   randum
REAL(8)					    RANUM,RANUM1
         
INTEGER    i,j,i_part
REAL(8)    xp,yp,dx,dy,xold
REAL(8)    bz,br,btp,r_larmor,phi1,phi2,phi3,t_new,v_new,a1,a2,a3,a4,v1
REAL(8)    bx, by   !$ ab.ZWZ 2021/9/11
    
    

	!CALL GTOP(i_part,bz,bfz)
	!CALL GTOP(i_part,br,bfy)

!$ ====== mb.ZWZ 2021/9/11 ==== \\
    CALL GTOP(i_part,bx,bfx)    
    CALL GTOP(i_part,by,bfy)
    CALL GTOP(i_part,bz,bfz)
!$ ====== mb.ZWZ 2021/9/11 ==== //

    CALL GTOP(i_part,phi2,Phi)
	    
	!btp=B_Ref*sqrt(bz**2+br**2)        !*bmax
    !btp=Bfield_ref*sqrt(bz**2+br**2)        !*bmax
    btp=Bfield_ref*sqrt(bx**2+by**2+bz**2)        !$ mb.ZWZ 2021/9/11

	v1=SQRT(part(i_part,4)**2+part(i_part,5)**2+part(i_part,6)**2)*v_ref

    CALL DRandom(RANUM)
	phi1=2.0*pi*ranum
    CALL DRandom(RANUM1)
!    r_larmor=0.3
    r_larmor=me*v1/(e*btp)/L_ref
    IF(r_larmor>0.3)  r_larmor=0.3
    
    !WRITE(*,*) 'btp,r_larmor',btp,r_larmor
    !STOP
    
    !part(i_part,2)=part(i_part,2)+r_larmor*cos(phi1)
    !xold=part(1,i_part)
    !part(i_part,1)=part(i_part,1)-abs(r_larmor*sin(phi1)) !向阳极牵移,牵移距离0到0.1之间
    
    part(i_part,1)=part(i_part,1)+r_larmor*cos(phi1)
    !xold=part(1,i_part)
    xold=part(i_part,1)     !$ mb.ZWZ 2021/7/10
    part(i_part,2)=part(i_part,2)-abs(r_larmor*sin(phi1)) !向阳极牵移,牵移距离0到0.1之间
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! think more abou here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !IF(part(1,i_part)<=boundaries(2)%Locations(2,1))THEN
    !   
    !    if(part(2,i_part)>boundaries(1)%Locations(2,2) .OR.     &
    !       part(2,i_part)<boundaries(1)%Locations(1,2))then
    !      part(2,i_part)=part(2,i_part)-r_larmor*cos(phi1)
    !    endif
    !    if(part(1,i_part)<boundaries(1)%Locations(1,1) .or. xold>boundaries(2)%Locations(2,1))then
    !     part(1,i_part)=part(1,i_part)+abs(r_larmor*sin(phi1))
    !    endif
    !    
    !ELSE 
    !    if(part(2,i_part)>boundaries(5)%Locations(2,2) .OR.     &
    !       part(2,i_part)<boundaries(1)%Locations(1,2))then
    !      part(2,i_part)=part(2,i_part)-r_larmor*cos(phi1)
    !    endif
    ! ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! think more abou here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

	CALL GTOP(i_part,phi3,Phi)
	
    t_new=ENER(i_part)+(phi3-phi2)*tmpj(1)  !! bjw
    !t_new=ENER(i_part)+(phi3-phi2)*temp0(1)
    
    if(t_new<500 .AND. t_new>0) then
        ENER(i_part)=t_new
        v_new=sqrt(t_new/t_parameter(1)) 
        CALL anisotropicv(v_new,part(i_part,4),part(i_part,5),part(i_part,6))
    else
        v_new=sqrt(ENER(i_part)/t_parameter(1)) 
        CALL anisotropicv(v_new,part(i_part,4),part(i_part,5),part(i_part,6))
    end if
    
    CALL PTOG(COND,i_part)

END SUBROUTINE
