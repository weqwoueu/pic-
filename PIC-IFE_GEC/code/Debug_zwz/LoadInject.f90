SUBROUTINE LoadInject(it, delta)
    
USE Domain_2D
USE Particle_2D
USE IMPIC_Data_2D
USE Constant_Variable_2D
USE TimeControl
USE PIC_MAIN_PARAM_2D
USE Field_2D
USE MDL_DEBUG   !$ ab.ZWZ for checking atom distribution 
IMPLICIT NONE
INTEGER		it, delta

REAL(8)       ::  b_amb(2,2)
INTEGER       ::  jj, ntotp, i_part, i
INTEGER       ::  add_N(ispe_inject)
!REAL(8)       ::  ranum
double precision ::  ranum, R_ranum, Theta_ranum
REAL(8)       ::  v_r, FS1, FS2, QA, UU, UN, GG, V_x, V_y, V_z
CHARACTER*25	:: fname, sname, filename
REAL(8)       ::  r1,r2,x, y, r
INTEGER ::count1 
INTEGER ::count2 


DO jj=1,ispe_inject 
    
    IF (vd(ipf(jj),1) > 0.0) THEN
        !!! 各种粒子在入射方向上的平均速度（速度为正）
        b_amb(ipf(jj),1)=vt(ipf(jj),1)/SQRT(2.*PI)*DEXP(-(vd(ipf(jj),1))**2/(2*vt(ipf(jj),1)**2))+ &
                         vd(ipf(jj),1)/2*(1+ERF(vd(ipf(jj),1)/(SQRT(2.)*vt(ipf(jj),1))))
    ELSE
        !!!! 各种粒子在入射方向上的平均速度（速度为负）
        b_amb(ipf(jj),1)=-vt(ipf(jj),1)/SQRT(2.*PI)*DEXP(-(vd(ipf(jj),1))**2/(2*vt(ipf(jj),1)**2))+ &
                         vd(ipf(jj),1)/2*(1-ERF(vd(ipf(jj),1)/(SQRT(2.)*vt(ipf(jj),1))))
    END IF
    
    CALL DRandom(ranum)
    !$ ============= mb.ZWZ =================== \\
    IF (delta == 0) THEN
        add_N(ipf(jj))=int(dens0(ipf(jj))*n_ref*ABS(b_amb(ipf(jj),1))*v_ref*dt*time_ref* &
                          (f_right_wall(2)-f_left_wall(2))*L_ref/affp_bjw(ipf(jj))+ranum) 
    ELSEIF (delta == 1) THEN
        add_N(ipf(jj))=int(dens0(ipf(jj))*n_ref*ABS(b_amb(ipf(jj),1))*v_ref*dt*time_ref* &
                          PI*(f_right_wall(2)**2-f_left_wall(2)**2)*L_ref**2/affp_bjw(ipf(jj))+ranum)  
        !add_N(ipf(jj))=int(dens0(ipf(jj))*n_ref*ABS(b_amb(ipf(jj),1))*v_ref*dt*time_ref* &
        !                  (f_right_wall(2)-f_left_wall(2))*L_ref/affp_bjw(ipf(jj))+ranum) 
        print*,add_N(ipf(jj))
        add_N(ipf(jj)) = 25
        !pause
    ENDIF
    !$ ============= mb.ZWZ =================== //
END DO
!add_N(1) == 100

DO jj = 1, 2
    dy_inject(ipf(jj)) = (f_right_wall(2)-f_left_wall(2))/(add_N(ipf(jj))) !$ ab.ZWZ for adding weighting
ENDDO
print*,dy_inject

print*,'before injection:'
print*,'ns(1)=',ns(1)
print*,'ns(2)=',ns(2)
count1 = 0
count2 = 0


    !!!! 电子
    DO jj =1,1
        
        WRITE(6,*) '# to inject on one side -- electron',add_N(ipf(jj))
	    ntotp=ntot
	    ntot = ntotp + add_N(ipf(jj))
	    ns(ipf(jj)) = ns(ipf(jj)) + add_N(ipf(jj))
	    !WRITE(6,*) 'injection at inlet, ntot=(eletron)', ntot, 'ns=', ns(jj)
	    WRITE(6,*) 'injection at inlet, ntot=', ntot, 'ns(eletron)=', ns(jj)

	    DO i_part=ntotp+1, ntot	     
	     
                part(i_part,:) = 0.
                
	            CALL DRandom(R_ranum)
                CALL DRandom(Theta_ranum)

                !part(i_part,1)= f_left_wall(1)+100+R_ranum*(f_right_wall(1)-100-f_left_wall(1)-100)           !$ x位置
                part(i_part,1) = 1
                IF(delta == 0) THEN
                    part(i_part,2)= Theta_ranum*(f_right_wall(2)-f_left_wall(2))                        !$ y位置
                ELSEIF(delta == 1) THEN
                    !part(i_part,2)= SQRT(f_left_wall(2)*f_left_wall(2) &
                    !                + (f_right_wall(2)+f_left_wall(2)) &
                    !                *(f_right_wall(2)-f_left_wall(2))*Theta_ranum)
                    
                    part(i_part,2) = (i_part-ntotp-0.5)*dy_inject(ipf(jj))
                    !print*,i_part
                    !print*,ntotp
                    !print*,part(i_part,2)
                    !pause
                    Rweight(i_part) = 2.*PI*part(i_part,2)*dy_inject(ipf(jj))
                    
                    
                    IF (part(i_part,2) < hx(2)/2) THEN
                        count1 = count1 + 1
                    ELSE IF(part(i_part,2) < hx(2)*1.5) THEN
                        count2 = count2 + 1
                    ENDIF
                    
                ENDIF
                part(i_part,7) = ipf(jj)
                

        END DO
    END DO

    WRITE(6,*) 'updated ntot=', ntot, 'updated ns(electron)=',ns(1)
    print*,' count1=',count1
    print*,' count2=',count2
    !pause
    
    
END