SUBROUTINE AdjustBoundary_2D(time, dt, delta, N_Objects, objects)
    
!!! bjw add 2018-05-23 调整运动到区域外的粒子

USE Particle_2D
USE Domain_2D
USE Object_Data_2D
USE IMPIC_Data_2D
USE PIC_MAIN_PARAM_2D !$ ab.ZWZ
Use ModuleMCCInterface,ONLY:ControlFlowGlobal, ParticleGlobal
IMPLICIT NONE
INTEGER    :: i_part,ispe
REAL(8)    :: x,y,z, xtemp, ytemp

REAL(8)    :: time, dt
INTEGER    :: delta
TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
INTEGER    :: N_Objects


REAL(8)     :: xp, yp, zp
REAL(8)     :: delt, V_n, V_t

Integer     :: isp,i,j

!Do isp = 1, ispe_tot
Do isp=0,ControlFlowGlobal%Ns
    !do i=ParticleGlobal(isp)%NPar,1,-1
    i = ParticleGlobal(isp)%NPar
    do while(i > 0)
        x = ParticleGlobal(isp)%PO(i)%X
        y = ParticleGlobal(isp)%PO(i)%Y
        ! asorb boundary
        If (x < dxmin .or. x > dxmax) then
            Call ParticleGlobal(isp)%DelOne(i) 
            goto 100
        End If
        
        ! asorb objects
        If (n_objects /=0) Then
            Do j = 1,N_Objects
                If (objects(j)%Shape == 3) Then
                    If (x > objects(j)%Locations(1,1) .And. x < objects(j)%Locations(2,1) .AND. &
                        y > objects(j)%Locations(1,2) .AND. y < objects(j)%Locations(2,2) ) THEN  
                        Call ParticleGlobal(isp)%DelOne(i)
                        goto 100
                    End If
                Elseif (objects(j)%Shape == 1) Then
                    If ( (x - objects(j)%Locations(1,1))**2 + (y - objects(j)%Locations(1,2))**2 < objects(j)%Dimensions(1)**2) Then
                        Call ParticleGlobal(isp)%DelOne(i)
                        goto 100
                    End If
                End If
            End Do
        End If
        ! reflex boundary
        If (y < dymin) then
            If (delta == 0) Then
                ParticleGlobal(isp)%PO(i)%Y = ParticleGlobal(isp)%PO(i)%Y - ParticleGlobal(isp)%PO(i)%Vy*dt
                ParticleGlobal(isp)%PO(i)%Vy = -ParticleGlobal(isp)%PO(i)%Vy
                i = i + 1
                goto 100
            Else If (delta == 1) Then
                Print*,'nearly impossible cross!'
                Print*,ParticleGlobal(isp)%PO(i)
                Stop
            End If
        End If
        If (y > dymax) then
            If (delta == 0) Then
                ParticleGlobal(isp)%PO(i)%Y = ParticleGlobal(isp)%PO(i)%Y - ParticleGlobal(isp)%PO(i)%Vy*dt
                ParticleGlobal(isp)%PO(i)%Vy = -ParticleGlobal(isp)%PO(i)%Vy
                i = i + 1
                goto 100                
            Else If (delta == 1) Then
                !delt = y - dymax
                !ParticleGlobal(isp)%PO(i)%Y = dymax - delt
                !
                !V_n =  ParticleGlobal(isp)%PO(i)%Vy*DCos(ParticleGlobal(isp)%PO(i)%Z) + ParticleGlobal(isp)%PO(i)%Vz*Dsin(ParticleGlobal(isp)%PO(i)%Z)
                !V_t = -ParticleGlobal(isp)%PO(i)%Vy*DSin(ParticleGlobal(isp)%PO(i)%Z) + ParticleGlobal(isp)%PO(i)%Vz*DCos(ParticleGlobal(isp)%PO(i)%Z)
                !V_n = -V_n
                !
                !ParticleGlobal(isp)%PO(i)%Vy = V_n*DCos(ParticleGlobal(isp)%PO(i)%Z) - V_t*Dsin(ParticleGlobal(isp)%PO(i)%Z)
                !ParticleGlobal(isp)%PO(i)%Vz = V_n*DSin(ParticleGlobal(isp)%PO(i)%Z) + V_t*DCos(ParticleGlobal(isp)%PO(i)%Z)
                !
                !i = i + 1
                !goto 100
                
                Call ParticleGlobal(isp)%DelOne(i) 
                goto 100
            End If
        End If
100     i = i - 1
    end do
End Do

!print*,'before AdjustBoundary:'
!print*,'ns(1)=',ns(1)
!print*,'ns(2)=',ns(2)

!i_part = 1
!DO WHILE(i_part <= ntot)
!    !write(*,*) i_part, part(i_part,:)
!    !!! 吸收边界
!    x = part(i_part,1)
!    y = part(i_part,2)
!    
!    IF (x < f_left_wall(1) .OR. x > f_right_wall(1)) THEN
!           
!        ispe = part(i_part,7)
!		ns(ispe) = ns(ispe) - 1
!        part(i_part,:) = part(ntot,:)
!        part(ntot,:) = 0.0
!            
!        ntot = ntot - 1
!        i_part = i_part - 1
!        GOTO 100   
!    END IF
!
!    !IF (SQRT(x*x+y*y) < r1 .OR. SQRT(x*x+y*y) > r2) THEN
!    IF (n_objects /= 0) THEN
!        
!        IF (x > objects(1)%Locations(1,1) .AND. x < objects(1)%Locations(2,1) .AND. &
!            y > objects(1)%Locations(1,2) .AND. y < objects(1)%Locations(2,2) ) THEN  !! 物体1
!           
!            ispe = part(i_part,7)
!		    ns(ispe) = ns(ispe) - 1
!            part(i_part,:) = part(ntot,:)
!            part(ntot,:) = 0.0
!            
!            
!            ntot = ntot - 1
!            i_part = i_part - 1
!            GOTO 100   
!        END IF
!           
!        IF (x > objects(2)%Locations(1,1) .AND. x < objects(2)%Locations(2,1) .AND. &
!            y > objects(2)%Locations(1,2) .AND. y < objects(2)%Locations(2,2) ) THEN  !! 物体2
!           
!            ispe = part(i_part,7)
!		    ns(ispe) = ns(ispe) - 1
!            part(i_part,:) = part(ntot,:)
!            part(ntot,:) = 0.0
!            
!            ntot = ntot - 1
!            i_part = i_part - 1
!            GOTO 100   
!        END IF
!    END IF
!    
!    !!!!! 反射边界
!    x = part(i_part,1)
!    y = part(i_part,2)
!    IF (y < f_left_wall(2)) THEN
!        
!        IF (delta == 0) THEN
!        
!            part(i_part,2) = part(i_part,2)-part(i_part,5)*dt
!	        part(i_part,5) = -part(i_part,5)
!            i_part = i_part - 1
!            
!        ELSEIF (delta == 1) THEN
!            
!            print*,'impossible cross'
!            print*,i_part
!            print*,part(i_part,:)
!            pause
!        
!            !*****transform positions to cartesian coordinates*****
!		    xp=part(i_part,2)*DCOS(part(i_part,3))
!		    yp=part(i_part,2)*DSIN(part(i_part,3))
!		    zp=part(i_part,1)
!            !*****update cartesian positions*****
!		    xp=xp-dt*part(i_part,5)
!		    IF(xp == 0) xp=hx(2)*1.0E-5
!		    yp=yp-dt*part(i_part,6)
!            !*****update polar positions
!		    part(i_part,2) = DSQRT(xp*xp+yp*yp)
!            part(i_part,3) = DATAN(yp/xp)
!            !*****place particle in proper quadrant*****
!		    IF (xp <= 0.0) THEN
!			    part(i_part,3)=part(i_part,3)+PI
!		    ENDIF
!        
!            part(i_part,5) = -part(i_part,5)
!            part(i_part,6) = -part(i_part,6)
!        
!            i_part = i_part - 1
!        
!        ENDIF
!        
!    ELSEIF (y > f_right_wall(2)) THEN
!        
!      !  IF (delta == 0) THEN
!      !  
!      !      part(i_part,2) = part(i_part,2)-part(i_part,5)*dt
!      !      part(i_part,5) = -part(i_part,5)
!      !      i_part = i_part - 1
!      !      
!      !  ELSEIF (delta == 1) THEN
!      !      
!      !!      ispe = part(i_part,7)
!		    !!ns(ispe) = ns(ispe) - 1
!      !!      part(i_part,:) = part(ntot,:)
!      !!      part(ntot,:) = 0.0
!      !!      
!      !!      ntot = ntot - 1
!      !!      i_part = i_part - 1
!      !!      GOTO 100 
!      !      
!      !      delt = part(i_part,2) - f_right_wall(2)
!      !      part(i_part,2) = f_right_wall(2) - delt
!      !      
!      !      V_n =  part(i_part,5)*DCOS(part(i_part,3)) + part(i_part,6)*DSIN(part(i_part,3))
!      !      V_t = -part(i_part,5)*DSIN(part(i_part,3)) + part(i_part,6)*DCOS(part(i_part,3))
!      !      V_n = -V_n
!      !      
!      !      part(i_part, 5) = V_n*DCOS(part(i_part, 3)) - V_t*DSIN(part(i_part, 3))
!      !      part(i_part, 6) = V_n*DSIN(part(i_part, 3)) + V_t*DCOS(part(i_part, 3))
!      !      
!      !      i_part = i_part - 1
!      !      
!      !  ENDIF
!        
!        ispe = part(i_part,7)
!		ns(ispe) = ns(ispe) - 1
!        part(i_part,:) = part(ntot,:)
!        part(ntot,:) = 0.0
!        ntot = ntot - 1
!        i_part = i_part - 1
!        GOTO 100 
!    END IF
!    
!     !IF (y < f_left_wall(2)) THEN    
!        
!  !          ispe = part(i_part,7)
!  !          
!  !          !Nloss_local1 = Nloss_local1 + 1
!  !          !IF (ispe==1) THEN
!  !          !    Neloss_local1 = Neloss_local1 + 1
!  !          !ELSEIF (ispe==2) THEN
!  !          !    Niloss_local1 = Niloss_local1 + 1
!  !          !ENDIF
!  !          
!		!    ns(ispe) = ns(ispe) - 1
!  !          part(i_part,:) = part(ntot,:)
!  !          part(ntot,:) = 0.0
!  !          
!  !          A_bar_n_minus_1(:,i_part) = A_bar_n_minus_1(:,ntot)    !!!! bjw 2019-4-19 修复此BUG
!  !          A_bar_n(:,i_part) = A_bar_n(:,ntot)
!  !          A_bar(:,i_part) = A_bar(:,ntot)
!  !          A_n_plus_1(:,i_part) = A_n_plus_1(:,ntot)
!  !          A_bar_n_minus_1(:,ntot) = 0.0
!  !          A_bar_n(:,ntot) = 0.0
!  !          A_bar(:,ntot) = 0.0
!  !          A_n_plus_1(:,ntot) = 0.0
!  !          
!  !          ntot = ntot - 1
!  !          i_part = i_part - 1
!  !          GOTO 100   
!  !  END IF
!  !  
!  !  IF (y > f_right_wall(2)) THEN     
!  !      
!  !      ispe = part(i_part,7)
!  !      
!  !      !Nloss_local2 = Nloss_local2 + 1
!  !      !IF (ispe==1) THEN
!  !      !    Neloss_local2 = Neloss_local2 + 1
!  !      !ELSEIF (ispe==2) THEN
!  !      !    Niloss_local2 = Niloss_local2 + 1
!  !      !ENDIF
!  !          
!		!ns(ispe) = ns(ispe) - 1
!  !      part(i_part,:) = part(ntot,:)
!  !      part(ntot,:) = 0.0
!  !          
!  !      A_bar_n_minus_1(:,i_part) = A_bar_n_minus_1(:,ntot)    !!!! bjw 2019-4-19 修复此BUG
!  !      A_bar_n(:,i_part) = A_bar_n(:,ntot)
!  !      A_bar(:,i_part) = A_bar(:,ntot)
!  !      A_n_plus_1(:,i_part) = A_n_plus_1(:,ntot)
!  !      A_bar_n_minus_1(:,ntot) = 0.0
!  !      A_bar_n(:,ntot) = 0.0
!  !      A_bar(:,ntot) = 0.0
!  !      A_n_plus_1(:,ntot) = 0.0
!  !          
!  !      ntot = ntot - 1
!  !      i_part = i_part - 1
!  !      GOTO 100   
!  !      
!  !  END IF
!    
!    !!!! 周期边界（左右边界）
!    !y = part(i_part,2)
!    !IF (x < f_left_wall(1)) THEN
!    !    
!	   ! part(i_part,1) = part(i_part,1) + (f_right_wall(1) - f_left_wall(1))
!    !    i_part = i_part - 1
!    !    
!    !ELSEIF (x > f_right_wall(1)) THEN
!    !    
!    !    part(i_part,1) = part(i_part,1) - (f_right_wall(1) - f_left_wall(1))
!    !    i_part = i_part - 1
!	   !
!    !END IF
!    !
!    
!
!    !!!! 吸收边界
!    !x = part(i_part,1)
!    !y = part(i_part,2)
!    !
!    !IF (SQRT(x*x+y*y) < r1 .OR. SQRT(x*x+y*y) > r2) THEN
!    !       
!    !        ispe = part(i_part,7)
!		  !  ns(ispe) = ns(ispe) - 1
!    !        part(i_part,:) = part(ntot,:)
!    !        part(ntot,:) = 0.0
!    !        
!    !        A_bar_n_minus_1(:,i_part) = A_bar_n_minus_1(:,ntot)    !!!! bjw 2019-4-19 修复此BUG
!    !        A_bar_n(:,i_part) = A_bar_n(:,ntot)
!    !        A_bar(:,i_part) = A_bar(:,ntot)
!    !        A_n_plus_1(:,i_part) = A_n_plus_1(:,ntot)
!    !        A_bar_n_minus_1(:,ntot) = 0.0
!    !        A_bar_n(:,ntot) = 0.0
!    !        A_bar(:,ntot) = 0.0
!    !        A_n_plus_1(:,ntot) = 0.0
!    !        
!    !        ntot = ntot - 1
!    !        i_part = i_part - 1
!    !        GOTO 100   
!    !END IF
!    
!  !  z = part(i_part,3)
!  !  IF (z < 0.0) THEN
!  !      
!	 !   ispe = part(i_part,7)
!		!ns(ispe) = ns(ispe) - 1
!  !      part(i_part,:) = part(ntot,:)
!  !      part(ntot,:) = 0.0
!  !          
!  !      A_bar_n_minus_1(:,i_part) = A_bar_n_minus_1(:,ntot)    !!!! bjw 2019-4-19 修复此BUG
!  !      A_bar_n(:,i_part) = A_bar_n(:,ntot)
!  !      A_bar(:,i_part) = A_bar(:,ntot)
!  !      A_n_plus_1(:,i_part) = A_n_plus_1(:,ntot)
!  !      A_bar_n_minus_1(:,ntot) = 0.0
!  !      A_bar_n(:,ntot) = 0.0
!  !      A_bar(:,ntot) = 0.0
!  !      A_n_plus_1(:,ntot) = 0.0
!  !          
!  !      ntot = ntot - 1
!  !      i_part = i_part - 1
!  !      GOTO 100   
!  !      
!  !  ELSEIF (z > 100.0*1.5) THEN
!  !      
!  !      ispe = part(i_part,7)
!		!ns(ispe) = ns(ispe) - 1
!  !      part(i_part,:) = part(ntot,:)
!  !      part(ntot,:) = 0.0
!  !          
!  !      A_bar_n_minus_1(:,i_part) = A_bar_n_minus_1(:,ntot)    !!!! bjw 2019-4-19 修复此BUG
!  !      A_bar_n(:,i_part) = A_bar_n(:,ntot)
!  !      A_bar(:,i_part) = A_bar(:,ntot)
!  !      A_n_plus_1(:,i_part) = A_n_plus_1(:,ntot)
!  !      A_bar_n_minus_1(:,ntot) = 0.0
!  !      A_bar_n(:,ntot) = 0.0
!  !      A_bar(:,ntot) = 0.0
!  !      A_n_plus_1(:,ntot) = 0.0
!  !          
!  !      ntot = ntot - 1
!  !      i_part = i_part - 1
!  !      GOTO 100   
!	 !  
!  !  END IF
!    
!    
!    
!    
!    !!!! 周期边界
!    !x = part(i_part,1)
!    !y = part(i_part,2)
!    !IF (x < 0.0) THEN
!    !    
!    !    xtemp = part(i_part,1)
!    !    ytemp = part(i_part,2)
!	   ! part(i_part,1) = ytemp
!    !    part(i_part,2) = -xtemp
!    !    i_part = i_part - 1
!    !    
!    !ELSEIF (y < 0.0) THEN
!    !    
!    !    xtemp = part(i_part,1)
!    !    ytemp = part(i_part,2)
!	   ! part(i_part,1) = -ytemp
!    !    part(i_part,2) = xtemp
!    !    i_part = i_part - 1
!	   !
!    !END IF
!    
!    
!    
!    
!    !!!! 吸收边界（上下边界）
!    !y = part(i_part,2)
!    !
!    !IF (n_objects /= 0) THEN
!    !
!    !    IF (y < boundaries(1)%Locations(1,2) .OR. y > f_right_wall(2)) THEN
!    !       
!    !        ispe = part(i_part,7)
!		  !  ns(ispe) = ns(ispe) - 1
!    !        part(i_part,:) = part(ntot,:)
!    !        part(ntot,:) = 0.0
!    !        
!    !        A_bar_n_minus_1(:,i_part) = A_bar_n_minus_1(:,ntot)    !!!! bjw 2019-4-19 修复此BUG
!    !        A_bar_n(:,i_part) = A_bar_n(:,ntot)
!    !        A_bar(:,i_part) = A_bar(:,ntot)
!    !        A_n_plus_1(:,i_part) = A_n_plus_1(:,ntot)
!    !        A_bar_n_minus_1(:,ntot) = 0.0
!    !        A_bar_n(:,ntot) = 0.0
!    !        A_bar(:,ntot) = 0.0
!    !        A_n_plus_1(:,ntot) = 0.0
!    !        
!    !        ntot = ntot - 1
!    !        i_part = i_part - 1
!    !        GOTO 100   
!    !    END IF
!    !    
!    !ELSE
!    !    
!    !    IF (y < f_left_wall(2) .OR. y > f_right_wall(2)) THEN
!    !        ispe = part(i_part,7)
!		  !  ns(ispe) = ns(ispe) - 1
!    !        part(i_part,:) = part(ntot,:)
!    !        part(ntot,:) = 0.0
!    !        
!    !        A_bar_n_minus_1(:,i_part) = A_bar_n_minus_1(:,ntot)    !!!! bjw 2019-4-19 修复此BUG
!    !        A_bar_n(:,i_part) = A_bar_n(:,ntot)
!    !        A_bar(:,i_part) = A_bar(:,ntot)
!    !        A_n_plus_1(:,i_part) = A_n_plus_1(:,ntot)
!    !        A_bar_n_minus_1(:,ntot) = 0.0
!    !        A_bar_n(:,ntot) = 0.0
!    !        A_bar(:,ntot) = 0.0
!    !        A_n_plus_1(:,ntot) = 0.0
!    !        
!    !        ntot = ntot - 1
!    !        i_part = i_part - 1
!    !        GOTO 100
!    !    END IF
!    !
!    !END IF
!    !
!    !!!! 周期边界（左右边界）
!    !x = part(i_part,1)
!    !IF (x < f_left_wall(1)) THEN
!    !    
!	   ! part(i_part,1) = part(i_part,1) + (f_right_wall(1) - f_left_wall(1))
!    !    i_part = i_part - 1
!    !    
!    !ELSEIF (x > f_right_wall(1)) THEN
!    !    
!    !    part(i_part,1) = part(i_part,1) - (f_right_wall(1) - f_left_wall(1))
!    !    i_part = i_part - 1
!	   !
!    !END IF
!    
!    
!    !!!! 周期边界（前后边界-2D3V）:设置300个德拜长度
!    !z = part(i_part,3)
!    !IF (z < 0.0) THEN
!    !    
!	   ! part(i_part,3) = part(i_part,3) + 100.0*10.0
!    !    ispe = part(i_part,7)
!    !    IF (ispe == 1) THEN !! 电子
!    !        CALL Loadv(V_z, tmpj(ipf(ispe)), ispe)
!    !        part(i_part,6) = 0.0 + V_z
!    !    ELSEIF (ispe == 2) THEN !! 离子
!    !        CALL Loadv(V_z, tmpj(ipf(ispe)), ispe)
!    !        part(i_part,6) = vd(ispe,3) + V_z
!    !    END IF
!    !    
!    !    i_part = i_part - 1
!    !    
!    !ELSEIF (z > 100.0*10.0) THEN
!    !    
!    !    part(i_part,3) = part(i_part,3) - 100.0*10.0
!    !    ispe = part(i_part,7)
!    !    IF (ispe == 1) THEN !! 电子
!    !        CALL Loadv(V_z, tmpj(ipf(ispe)), ispe)
!    !        part(i_part,6) = 0.0 + V_z
!    !    ELSEIF (ispe == 2) THEN !! 离子
!    !        CALL Loadv(V_z, tmpj(ipf(ispe)), ispe)
!    !        part(i_part,6) = vd(ispe,3) + V_z
!    !    END IF
!    !    
!    !    i_part = i_part - 1
!	   !
!    !END IF
!    
!    
!    
!100 continue
!    i_part = i_part + 1
!
!END DO

!WRITE(*,*) Neloss_local1,Niloss_local1
!WRITE(*,*) Neloss_local2,Niloss_local2
!print*,'after AdjustBoundary:'
!print*,'ns(1)=',ns(1)
!print*,'ns(2)=',ns(2)
    
END