SUBROUTINE  Ionization_2D(i_part, energy, ionized_particle, IC)

! Updated:	11/29/2004 05:54 PM
! Purpose:	Check particles exiting through the domain outer boundaries.

USE MCC_Data_2D
USE MCC_Main_Param

USE IFE_Data
USE DSMC_Data_2D

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Particle_2D

USE Constant_Variable_2D

IMPLICIT NONE

!double precision randum
!external randum

REAL(8) 	ranum

REAL(8)		energy, rengy, ta
INTEGER		i_part,i,iionization,ia

INTEGER		new_electron,new_ion, ionized_particle
INTEGER		ipret, jpret, mc, sum_fnuma, ifnuma

REAL(8)		ve0, v, vnew, fnuma_res

INTEGER(1), DIMENSION(:), ALLOCATABLE	::	pc
INTEGER, DIMENSION(2,(nx-1)*(ny-1))		::	IC

!WRITE(*,*) '11111111111111'
!STOP

    
!原电子和新电子的能量变化		
!	print *,'energy=',energy
	energy=energy-12.1
!	print *,'energy=',energy
    CALL DRandom(RANUM)
	rengy=abs(8.7*tan(RANUM*atan(energy/(2*8.7)))) !new electron's energy
	energy=abs(energy-rengy) !old electron's energy
    ve0=SQRT(energy/t_parameter(1))
!	call newvel(energy,vr(insp,j),vtheta(j),vz(insp,j),0)
    CALL Anisotropicv(ve0,part(i_part,4),part(i_part,5),part(i_part,6))

!新产生的电子的速度和位置
    ionized_particle=ionized_particle+1
    ntot=ntot+1 
    ns(1)=ns(1)+1
    new_electron=ntot                   !+2*ionized_particle-1
	vnew=sqrt(rengy/t_parameter(1))
	
!	print *,'energy,rengy,ve0,vnew=',energy,rengy,ve0,vnew

!WRITE(*,*) i_part, part(i_part,1:3)

	CALL Anisotropicv(vnew,part(new_electron,4),part(new_electron,5),part(new_electron,6))
	ENER(new_electron)=rengy
	part(new_electron,1)=part(i_part,1)
	part(new_electron,2)=part(i_part,2)
    part(new_electron,3)=part(i_part,3)
    part(new_electron,7)=1
    part(new_electron,11)=1    !!不考虑原子运动时无效
	part(new_electron,12)=part(i_part,12)

!新产生的离子速度和位置
	IF(ispe_tot>=2) THEN !add a new ion
	    ntot=ntot+1
	    ns(2)=ns(2)+1
	    new_ion=ntot            !+2*ionized_particle
	    part(new_ion,1)=part(i_part,1)
	    part(new_ion,2)=part(i_part,2)
        part(new_ion,3)=part(i_part,3)
        part(new_ion,7)=2
        part(new_ion,11)=1     !!不考虑原子运动时无效
		part(new_ion,12)=part(i_part,2)
        
!        if(mod(iloop,iont)/=0) then !当没计算离子运动时，离子密度也没有计算，此时电离造成离子密度的变化
!  			para=1.
!			call ptogt(rot,nsp,imax1,imax2,new_ion,insp,para,ipre,jpre,p1,p2,p3,p4)	     
!        end if
		ta=0.03 !如果没有考虑原子，则碰撞后的离子的平均速度为0.03eV

!如果考虑了原子，则原子的权重减一,如果一个原子的权重不足在0到1之间，则该原子权重变为0，并接着寻找下一个原子
!		if(ispe_tot>=3 .AND. iatom==1) then !minus atom
!!			inspa=3		
!			ipret = INT((part(i_part,1) - Vert_o(1)) / hx(1))
!			jpret = INT((part(i_part,2) - Vert_o(2)) / hx(2))
!			mc=ipret*(ny-1)+jpret+1
!			allocate(pc(ic(2,mc))) !将权重大于0的原子序号统一保存在pc中
!            sum_fnuma=0
!            ifnuma=0
!            do i=1,ic(2,mc)
!                ia=Ic(1,mc)+i
!                if(part(ir(ia),11)>0) then        !if(fnuma(ir(ia))>0) then
!                    ifnuma=ifnuma+1
!                    pc(ifnuma)=ir(ia)
!                    sum_fnuma=sum_fnuma+part(pc(ifnuma),11)
!                end if
!            end do
!            if(sum_fnuma>1) then
!                iionization=1
!                fnuma_res=1
!                do while(iionization==1)
!                    iionization=0
!                    call random_number(ranum)
!                    if(ifnuma>0) then                    
!                        ia=ranum*ifnuma+1
!                        ta=ENER(pc(ia)) !如果模拟了原子则产生的离子继承原子的温度
!                        if(part(pc(ia),11)>=fnuma_res) then
!			                part(pc(ia),11)=part(pc(ia),11)-fnuma_res
!			                fnuma0=part(pc(ia),11)
!			                part(pc(ia),11)=-1
!!			                !密度相应地减小
!!                            para=-1.
!                            call PTOG(rho_s(:,:,3),pc(ia),1)
!!			                call ptogt(rot,nsp,imax1,imax2,pc(ia),inspa,para,ipre,jpre,p1,p2,p3,p4)		            
!			                part(pc(ia),11)=fnuma0
!			            else if(part(pc(ia),11)<fnuma_res .AND. part(pc(ia),11)>0) then
!!			                para=-1.
!                            part(pc(ia),11)=-1
!                            call PTOG(rho_s(:,:,3),pc(ia),1)
!!			                call ptogt(rot,nsp,imax1,imax2,pc(ia),inspa,para,ipre,jpre,p1,p2,p3,p4)	    
!			                fnuma_res=fnuma_res-part(pc(ia),11) !剩余的需要减掉的原子权重
!			                part(pc(ia),11)=0
!			                iionization=1
!			                pc(ia)=pc(ifnuma)
!			                ifnuma=ifnuma-1
!			                !密度相应地减小			            
!			            end if
!!                        r_ini(2,new_ion)=r_ini(inspa,pc(ia))
! !                       r_ini(1,new_electron)=r_ini(inspa,pc(ia))
!			        end if              
!                end do
!!            else
! !               ionized_particle=ionized_particle-1 !如果网格内总的原子的权重小于1，则此次电离无效			    			    
!            end if
!            if(ta<=0) ta=0.03 !如果计算得到的原子温度等于零，则赋予原子一个温度，避免出现计算错误
!            deallocate(pc)
!	    end if
	    !离子继承原子的动能
        if(ta>10) then
            write(*,*) "ta>10"
        !    pause
        end if

        
	    V = sqrt(1.5*ta/t_parameter(2)) !由于只计算了原子轴向和径向两个方向的动能，三个方向的动能需要乘1.5
	    call anisotropicv(v,part(new_ion,4),part(new_ion,5),part(new_ion,6))
        !CALL Loadv(part(new_ion,4), 1.5*ta, 2)   
        !CALL Loadv(part(new_ion,5), 1.5*ta, 2)
        !CALL Loadv(part(new_ion,6), 1.5*ta, 2)
!v1=v1+sqrt(2*e*10/mi)/veth
!		vz(2,new_ion)=v1
        !part(new_ion,4)=part(new_ion,4)+vaz/Vel_Ref !加上原子的定向速度
        part(new_ion,4)=part(new_ion,4)+vaz/V_Ref !加上原子的定向速度
        !part(new_ion,5)=part(new_ion,5)+vaz/V_Ref !加上原子的定向速度
        !print*, part(new_ion,4)
	    ENER(new_ion)=t_parameter(2)*(part(new_ion,4)**2+part(new_ion,5)**2+part(new_ion,6)**2)
!	    if(trajectory_diag==1 .AND. particle_diag<ipc_num) then
!	        particle_diag=particle_diag+1
!	        ir1(particle_diag)=new_ion
!	    end if

	END IF

    IF(ispe_tot>=3 .AND. iatom==1 .and. sum_fnuma<=1) then
        ionized_particle=ionized_particle-1 !如果网格内总的原子的权重小于1，则此次电离无效
        CALL Remove_2D(new_electron,ntot)
        CALL Remove_2D(new_ion,ntot)
    ENDIF
                
                
END