SUBROUTINE	Quasi_2D(delta,it)

USE Domain_2D
USE Particle_2D
USE IMPIC_Data_2D
USE Constant_Variable_2D
USE TimeControl
USE PIC_MAIN_PARAM_2D
USE Field_2D

USE MDL_DEBUG   !$ ab.ZWZ for checking atom distribution 
IMPLICIT NONE
INTEGER,INTENT(IN) :: delta !$ ab.ZWZ
INTEGER       :: ii, jj, i_part, ntotp
!REAL(8)       :: ranum
double precision ::  ranum
REAL(8)       :: v_r, FS1, FS2, QA, UU, UN, GG, V_x, V_y, V_z
REAL(8)		  :: ns_local(2)
INTEGER       :: add_number(2)
REAL(8)       ::  b_amb(2,2)

INTEGER, DIMENSION(2,(nx-1)*(ny-1))		::	ic1, IC
INTEGER									::	n_cell_in
INTEGER									::	insp,i,nnx,nny
INTEGER									::	ipre, jpre
INTEGER, DIMENSION(nx-1,ny-1)			::	net_at_cell
REAL(8), DIMENSION(nx-1,ny-1)			::	pro_at_cell
INTEGER									::	net_at_boundary = 0
INTEGER									::	mc, num_ion, num_elec, num, ne_old


INTEGER     :: it
CHARACTER*40	fname
INTEGER, DIMENSION(nx-1)    :: nab

print*,'before quasi:'
print*,'ns(1)=',ns(1)
print*,'ns(2)=',ns(2)

!!! ******************************* 方法一（bjw 2019-5-7）: 各个网格准中性条件 ***********************************************
ne_old = ns(1)
!n_cell_in = 10 / hx(2)  !!! 准中性区域的网格厚度
n_cell_in = 10 / hx(1)  !!! 准中性区域的网格厚度

nnx = nx - 1
nny = ny - 1

DO insp = 1, 2
	CALL INDEXM(insp, IC, nnx, nny)  !!! nnx, nny网格数
	ic1(insp,:) = IC(2,:)       !$ 每个网格的粒子数
ENDDO

DO ipre=nnx - n_cell_in, nnx  !!!x方向入射(x最大位置处)
	net_at_cell = 0
	net_at_boundary = 0
	pro_at_cell = 0
	DO jpre = 1, nny
		mc = (ipre - 1) * nny + jpre
		net_at_cell(ipre,jpre) = ic1(2,mc) - ic1(1,mc)
		net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
	ENDDO

	IF(net_at_boundary > 0) THEN
		net_at_boundary = 0
		DO jpre = 1, nny
			IF(net_at_cell(ipre,jpre) > 0) THEN
				pro_at_cell(ipre,jpre) = REAL(net_at_cell(ipre,jpre))
				net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
			ELSE
				pro_at_cell(ipre,jpre) = 0.
			ENDIF
		ENDDO

		!pro_at_cell(:,jpre) = pro_at_cell(:,jpre) / REAL(net_at_boundary)
        pro_at_cell(ipre,:) = pro_at_cell(ipre,:) / REAL(net_at_boundary)

		DO jpre = 2, nny
			pro_at_cell(ipre,jpre) = pro_at_cell(ipre,jpre) + pro_at_cell(ipre,jpre-1)
		ENDDO

		DO i = 1, net_at_boundary
			num_elec = ns(1)
			num_ion = ns(2)
			num = num_elec + num_ion
            
            CALL DRandom(ranum)
			jpre = 1
			DO WHILE(ranum > pro_at_cell(ipre, jpre) .AND. jpre < nnx)
				jpre = jpre + 1
			ENDDO

			CALL DRandom(ranum)
			part(num+1,1) = 0. + (ipre - ranum) * hx(1)			
			CALL DRandom(ranum)
            !$ ========================= mb.ZWZ ================================= \\
            IF(delta == 0) THEN
                part(num+1,2)= f_left_wall(2) + (jpre - ranum) * hx(2)                      !$ y位置
                part(num+1,3) = 0.
            ELSEIF(delta == 1) THEN
                part(num+1,2)= f_left_wall(2) + SQRT(hx(2)**2*((jpre -1)*(jpre -1) &
                                +((jpre)+(jpre -1)) &
                                *((jpre)-(jpre -1))*ranum))
                CALL DRandom(ranum)
                part(num+1,3) = 2.*PI*ranum
            ENDIF
            !$ ========================= mb.ZWZ ================================= //

2000        CALL Loadv(V_x, tmpj(ipf(1)), 1)   !!! 只补电子，所以用 temp0(ipf(1))
		    CALL Loadv(V_y, tmpj(ipf(1)), 1)
            CALL Loadv(V_z, tmpj(ipf(1)), 1)
			IF(V_x >= 0) THEN
                !WRITE(*,*) '++++++++++++'
				goto 2000
			ENDIF
			part(num+1,4) = V_x
			part(num+1,5) = V_y
			part(num+1,6) = V_z

            part(num+1,8) = 0.
			part(num+1,9) = 0.
			part(num+1,10) = 0.
            
			part(num+1,7) = 1

			ns(1) = ns(1) + 1
		ENDDO
    ENDIF
ENDDO
WRITE(6,*) '# Quasi: to inject on one side',ns(1)-ne_old

ntot = ns(1) + ns(2)


!DO insp = 1, 2
!	CALL INDEXM(insp, IC, nnx, nny)  !!! nnx, nny网格数
!	ic1(insp,:) = IC(2,:)       !$ 每个网格的粒子数
!ENDDO
!net_at_cell = 0
!nab = 0
!DO ipre=1, 1 + n_cell_in  !!!x方向入射(x最小位置处)
!    net_at_boundary = 0
!    DO jpre = 1, nny
!		mc = (ipre - 1) * nny + jpre
!        net_at_cell(ipre,jpre) = ic1(2,mc) - ic1(1,mc) 
!        net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
!    ENDDO
!    nab(ipre) = net_at_boundary
!ENDDO
!DO ipre=nnx - n_cell_in, nnx  !!!x方向入射(x最大位置处)
!    net_at_boundary = 0
!    DO jpre = 1, nny
!		mc = (ipre - 1) * nny + jpre
!        net_at_cell(ipre,jpre) = ic1(2,mc) - ic1(1,mc)
!    ENDDO
!    nab(ipre) = net_at_boundary
!ENDDO
!
!IF(MOD(it,1000) == 0 ) THEN
!    WRITE(fname,"(I6.6)") it
!    OPEN(500,ACTION='WRITE', FILE='quasi'//TRIM(fname)//'.dat',POSITION='APPEND', STATUS='REPLACE')
!    WRITE(500,*) 'TITLE = "Field Plot"'
!    WRITE(500,"(A100)")'VARIABLES = "x" "y" "net_at_cell" "net_at_boundary"'
!    WRITE(500,"('ZONE I = ',I6,', J= ', I6)")nnx,nny
!       DO jpre=1,nny
!		    DO ipre=1,nnx
!                WRITE(500,"(2(E15.6),2(I15.6))")VertX(1:2,ipre,jpre),net_at_cell(ipre,jpre),nab(ipre)           
!		    END DO
!	    END DO
!    CLOSE(500)
!ENDIF


!DO jpre=nny - n_cell_in, nny  !!!y方向入射
!	net_at_cell = 0
!	net_at_boundary = 0
!	pro_at_cell = 0
!	DO ipre = 1, nnx
!		mc = (ipre - 1) * nny + jpre
!		net_at_cell(ipre,jpre) = ic1(2,mc) - ic1(1,mc)
!		net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
!	ENDDO
!
!	IF(net_at_boundary > 0) THEN
!		net_at_boundary = 0
!		DO ipre = 1, nnx
!			IF(net_at_cell(ipre,jpre) > 0) THEN
!				pro_at_cell(ipre,jpre) = REAL(net_at_cell(ipre,jpre))
!				net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
!			ELSE
!				pro_at_cell(ipre,jpre) = 0.
!			ENDIF
!		ENDDO
!
!		pro_at_cell(:,jpre) = pro_at_cell(:,jpre) / REAL(net_at_boundary)
!
!		DO ipre = 2, nnx
!			pro_at_cell(ipre,jpre) = pro_at_cell(ipre,jpre) + pro_at_cell(ipre-1,jpre)
!		ENDDO
!
!		DO i = 1, net_at_boundary
!			num_elec = ns(1)
!			num_ion = ns(2)
!			num = num_elec + num_ion
!            
!            CALL DRandom(ranum)
!			ipre = 1
!			DO WHILE(ranum > pro_at_cell(ipre, jpre) .AND. ipre < nny)
!				ipre = ipre + 1
!			ENDDO
!
!			CALL DRandom(ranum)
!			part(num+1,1) = 0. + (ipre - ranum) * hx(1)			
!			CALL DRandom(ranum)
!			part(num+1,2) = 0. + (jpre - ranum) * hx(2)
!			part(num+1,3) = 0.
!
!			CALL Loadv(V_x, tmpj(ipf(1)), 1)   !!! 只补电子，所以用 temp0(ipf(1))
!100         CALL Loadv(V_y, tmpj(ipf(1)), 1)
!            CALL Loadv(V_z, tmpj(ipf(1)), 1)
!			!IF(V_y >= 0) THEN
!			!	goto 100
!			!ENDIF
!			part(num+1,4) = V_x
!			part(num+1,5) = V_y
!			part(num+1,6) = V_z
!
!            part(num+1,8) = 0.
!			part(num+1,9) = 0.
!			part(num+1,10) = 0.
!            
!			part(num+1,7) = 1
!
!			ns(1) = ns(1) + 1
!		ENDDO
!	ENDIF
!ENDDO
!
!WRITE(6,*) '# Quasi: to inject on one side',ns(1)-ne_old
!
!ntot = ns(1) + ns(2)

!!! ******************************* 方法一（bjw 2019-5-7）: 各个网格准中性条件 ***********************************************



!!! ******************************* 方法二（bjw 2019-5-7）: 全区域准中性条件 *************************************************

!!!! 各种粒子在入射方向上的平均速度（速度为正）
!!b_amb(1,2)=vt(1,2)/SQRT(2.*PI)*DEXP(-(vd(1,2))**2/(2*vt(1,2)**2))+vd(1,2)/2*(1+ERF(vd(1,2)/(SQRT(2.)*vt(1,2))))
!!b_amb(2,2)=vt(2,2)/SQRT(2.*PI)*DEXP(-(vd(2,2))**2/(2*vt(2,2)**2))+vd(2,2)/2*(1+ERF(vd(2,2)/(SQRT(2.)*vt(2,2))))
!!WRITE(*,*) b_amb(2,1),b_amb(2,2)
!!!!! 各种粒子在入射方向上的平均速度（速度为负）
!!b_amb(1,2)=-vt(1,2)/SQRT(2.*PI)*DEXP(-(vd(1,2))**2/(2*vt(1,2)**2))+vd(1,2)/2*(1-ERF(vd(1,2)/(SQRT(2.)*vt(1,2))))
!!b_amb(2,2)=-vt(2,2)/SQRT(2.*PI)*DEXP(-(vd(2,2))**2/(2*vt(2,2)**2))+vd(2,2)/2*(1-ERF(vd(2,2)/(SQRT(2.)*vt(2,2))))
!
!DO jj=1,ispe_inject 
!    
!    IF (vd(ipf(jj),2) > 0.0) THEN
!        !!! 各种粒子在入射方向上的平均速度（速度为正）
!        b_amb(ipf(jj),2)=vt(ipf(jj),2)/SQRT(2.*PI)*DEXP(-(vd(ipf(jj),2))**2/(2*vt(ipf(jj),2)**2))+ &
!                         vd(ipf(jj),2)/2*(1+ERF(vd(ipf(jj),2)/(SQRT(2.)*vt(ipf(jj),2))))
!    ELSE
!        !!!! 各种粒子在入射方向上的平均速度（速度为负）
!        b_amb(ipf(jj),2)=-vt(ipf(jj),2)/SQRT(2.*PI)*DEXP(-(vd(ipf(jj),2))**2/(2*vt(ipf(jj),2)**2))+ &
!                         vd(ipf(jj),2)/2*(1-ERF(vd(ipf(jj),2)/(SQRT(2.)*vt(ipf(jj),2))))
!    END IF
!
!END DO
!!************************统计入口处离子比电子多的个数*********************************************************************
!DO jj=1,2
!	ns_local(jj)=0
!	add_number(jj)=0
!END DO
!
!DO i_part=1, ntot
!    
!	!IF(part(2,i_part)>=domain%ymin .and. part(2,i_part)<=domain%ymin+10.)THEN    !!! 此处7.0个德拜长度为准中性区域
!    IF(part(i_part,2)>=f_right_wall(2)-10. .and. part(i_part,2)<=f_right_wall(2))THEN    !!! 此处40.0个德拜长度为准中性区域
!	
!		ii=part(i_part,7)
!		ns_local(ii)=ns_local(ii)+1
!
!	END IF
!
!END DO
!
!print*,'ns_local',ns_local(1),ns_local(2)
!
!DO jj =1, 2
!    ns_local(jj)=affp_bjw(ipf(jj))*ns_local(jj)
!END DO
!
!IF(ns_local(1)-ns_local(2)>=0) THEN	
!	add_number(1)=0 !ns_local(1)+ns_local(3)-ns_local(2)-N_ina(2,6)
!ELSE
!    CALL DRandom(ranum)
!	add_number(1)=int(-(ns_local(1)-ns_local(2))/affp_bjw(1)+ranum)
!END IF
!
!!!!*******************************注入电子*******************************************
!DO jj=1,1
!    
!    WRITE(6,*) '# to inject on one side',add_number(jj)
!    ntotp = ntot
!	ntot = ntotp + add_number(jj)
!	ns(jj) = ns(jj) +add_number(jj)	
!	WRITE(6,*) 'injection at inlet, ntot=(eletron)', ntot, 'ns=', ns(jj)
!
!	DO i_part=ntotp+1, ntot
!            
!        !CALL Loadv(V_x, tmpj(ipf(jj)), 1)
!        !!CALL Loadv(V_y, temp0(ipf(jj)), 1)
!        !
!        !!!! 此过程可参考G.A.Bird书中的公式（12.5）
!        !v_r=abs(vd(ipf(jj),2))/(sqrt(2.)*vt(ipf(jj),2))
!        !FS1=v_r+SQRT(v_r*v_r+2.)
!        !FS2=0.5*(1.+v_r*(2.*v_r-FS1))
!        !    
!        !QA=3.			 
!        !IF (v_r.LT.-3.) QA=ABS(v_r)+1.
!        !103  CALL DRandom(ranum)	
!        !UU=-QA+2.*QA*ranum
!        !UN=UU+v_r
!        !!!!*--UN is a potential inward velocity component
!        !IF (UN.LT.0.) GO TO 103
!        !GG=(2.*UN/FS1)*DEXP(FS2-UU*UU)
!        !CALL DRandom(ranum)
!        !IF (GG.LT.ranum) GO TO 103
!        !!!!*--the inward normalised vel. component has been selected (eqn (12.5))
!        !!V_y=UN*(sqrt(2.)*vt(ipf(jj),2))
!        !V_y=-UN*(sqrt(2.)*vt(ipf(jj),2))   !!!! bjw
!          
!        
!        CALL Loadv(V_x, tmpj(ipf(1)), 1)   !!! 只补电子，所以用 temp0(ipf(1))
!100     CALL Loadv(V_y, tmpj(ipf(1)), 1)
!        IF(V_y >= 0) THEN
!            goto 100
!        ENDIF
!
!        part(i_part,4) = V_x
!        part(i_part,5) = V_y
!		part(i_part,6) = 0.0
!
!		part(i_part,8) = 0.0	
!        part(i_part,9) = 0.0
!        part(i_part,10) = 0.0
!
!        part(i_part,7) = ipf(jj)
!
!	    CALL DRandom(ranum)
!	    part(i_part,1)= (f_right_wall(1)-f_left_wall(1)) * ranum
!        CALL DRandom(ranum)
!	    part(i_part,2)= R_start(ipf(jj)) + 0.5 * b_amb(ipf(jj),2) * dt * ranum
!	    part(i_part,3)= 50.0
!             
!    END DO
!END DO

!!! ******************************* 方法二（bjw 2019-5-7）: 全区域准中性条件 *************************************************



END