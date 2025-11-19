SUBROUTINE PrePush
  
USE Domain_2D
USE Particle_2D
USE Field_2D
USE TimeControl
USE IMPIC_Data_2D
USE PIC_MAIN_PARAM_2D !$ ab.ZWZ
Use ModuleMCCInterface,ONLY:ControlFlowGlobal, ParticleGlobal   !$ ab.ZWZ for using JW's particle data structure
IMPLICIT NONE
INTEGER   :: i, j, i_part, ispe
REAL(8)   :: SVelocity(1:3)=0.d0
REAL(8)   :: xp, dx, yp, dy, f, g

!$ ========= ab.ZWZ 2021/7/11 ====== \\
Integer   :: isp
REAL(8)   :: xcellmdx, ycellmdy
REAL(8)   :: R, R1, R2, den
REAL(8)   :: P1, P2, P3, P4

REAL(8)   :: zefield, refield, tefield
REAL(8)	  :: xefield, yefield
REAL(8)   :: zbfield, rbfield, tbfield
REAL(8)	  :: xbfield, ybfield

!$ cylindrical moving
REAL(8)   :: X, Y, Z, Theta
!$ ========= ab.ZWZ 2021/7/11 ====== //

Do isp=0,ControlFlowGlobal%Ns
    Do i_part = 1, ParticleGlobal(isp)%NPar
        IF (Bfiled_index) THEN   !!!! 有磁场
            xp = (ParticleGlobal(isp)%PO(i_part)%X - Vert_o(1))*hxi(1)
	        i = xp        !!! 粒子在第i个单元里
	        dx = xp-i
    
            yp = (ParticleGlobal(isp)%PO(i_part)%Y - Vert_o(2))*hxi(2)
	        j = yp        !!! 粒子在第j个单元里
	        dy = yp-j
        
            xcellmdx = 1. -dx
            ycellmdy = 1. -dy
        
            IF(delta_global == 0) THEN
                P1 = xcellmdx*ycellmdy
                P2 = dx*ycellmdy
                P3 = xcellmdx*dy
                P4 = dx*dy
            ELSEIF(delta_global == 1)THEN
                R1=dymin + float(j - 1)*hx(2)
	            R2=dymin + float(j)*hx(2)
                R = ParticleGlobal(isp)%PO(i_part)%Y
	            den=R2*R2-R1*R1
            
                P1 = xcellmdx*(R2*R2-R*R)/den
                P2 = dx*(R2*R2-R*R)/den
                P3 = xcellmdx*(R*R-R1*R1)/den
                P4 = dx*(R*R-R1*R1)/den          
            ENDIF
        
            IF (delta_global == 0) THEN
        
                xbfield=bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4		
                ybfield=bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4		
                zbfield=bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4
            
                Omega(1) = xbfield * qm(isp+1) * 0.5 * dt
                Omega(2) = ybfield * qm(isp+1) * 0.5 * dt
                Omega(3) = zbfield * qm(isp+1) * 0.5 * dt
        
            ELSE IF (delta_global == 1) THEN
                !$ axisymmetric e field
                zbfield = bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4  
                rbfield = bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4   
                tbfield = bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4  
                !$ Convert axisymmetric efield to Cartesian efield
                Theta = ParticleGlobal(isp)%PO(i_part)%Z
                xbfield = rbfield*DCOS(Theta) - tbfield*DSIN(Theta)
                ybfield = rbfield*DSIN(Theta) + tbfield*DCOS(Theta)
            
                Omega(1) = zbfield * qm(isp+1) * 0.5 * dt
                Omega(2) = xbfield * qm(isp+1) * 0.5 * dt
                Omega(3) = ybfield * qm(isp+1) * 0.5 * dt
            ENDIF         
            
            SVelocity(1) = ParticleGlobal(isp)%PO(i_part)%Vx + 0.5 * ParticleGlobal(isp)%PO(i_part)%Ax * dt + &
                           ParticleGlobal(isp)%PO(i_part)%Vy * Omega(3) - ParticleGlobal(isp)%PO(i_part)%Vz * Omega(2)
            SVelocity(2) = ParticleGlobal(isp)%PO(i_part)%Vy + 0.5 * ParticleGlobal(isp)%PO(i_part)%Ay * dt + &
                           ParticleGlobal(isp)%PO(i_part)%Vz * Omega(1) - ParticleGlobal(isp)%PO(i_part)%Vx * Omega(3)
            SVelocity(3) = ParticleGlobal(isp)%PO(i_part)%Vz + 0.5 * ParticleGlobal(isp)%PO(i_part)%Az * dt + &
                           ParticleGlobal(isp)%PO(i_part)%Vx * Omega(2) - ParticleGlobal(isp)%PO(i_part)%Vy * Omega(1)
        
            OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
            TransB(1,1)=1.d0+Omega(1)*Omega(1)
            TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
            TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
            TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
            TransB(2,2)=1.d0+Omega(2)*Omega(2)
            TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
            TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
            TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
            TransB(3,3)=1.d0+Omega(3)*Omega(3)
            TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT
           
            ParticleGlobal(isp)%PO(i_part)%Vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
            ParticleGlobal(isp)%PO(i_part)%Vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
            ParticleGlobal(isp)%PO(i_part)%Vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)
        
        ELSE  !!!! 无磁场
    
            ParticleGlobal(isp)%PO(i_part)%Vx = ParticleGlobal(isp)%PO(i_part)%Vx + 0.5*dt*ParticleGlobal(isp)%PO(i_part)%Ax
            ParticleGlobal(isp)%PO(i_part)%Vy = ParticleGlobal(isp)%PO(i_part)%Vy + 0.5*dt*ParticleGlobal(isp)%PO(i_part)%Ay
            ParticleGlobal(isp)%PO(i_part)%Vz = ParticleGlobal(isp)%PO(i_part)%Vz + 0.5*dt*ParticleGlobal(isp)%PO(i_part)%Az
    
        END IF
    
        IF (delta_global == 0) THEN
            ParticleGlobal(isp)%PO(i_part)%X = ParticleGlobal(isp)%PO(i_part)%X + dt*ParticleGlobal(isp)%PO(i_part)%Vx
            ParticleGlobal(isp)%PO(i_part)%Y = ParticleGlobal(isp)%PO(i_part)%Y + dt*ParticleGlobal(isp)%PO(i_part)%Vy
    !Debug, Z is meaningless for 2-D planar problem
            ParticleGlobal(isp)%PO(i_part)%Z = ParticleGlobal(isp)%PO(i_part)%Z + dt*ParticleGlobal(isp)%PO(i_part)%Vz
        ELSEIF (delta_global == 1) THEN
            
    !------ transform positions to cartesian coordinates ------
		    X=ParticleGlobal(isp)%PO(i_part)%Y * DCOS(ParticleGlobal(isp)%PO(i_part)%Z)
		    Y=ParticleGlobal(isp)%PO(i_part)%Y * DSIN(ParticleGlobal(isp)%PO(i_part)%Z)
		    Z=ParticleGlobal(isp)%PO(i_part)%X
    !------ update cartesian positions ---------
		    X=X+dt * ParticleGlobal(isp)%PO(i_part)%Vy
            IF(X == 0.) THEN
                X=hx(2)*1.0E-5
            ENDIF
		    Y=Y+dt * ParticleGlobal(isp)%PO(i_part)%Vz
            
            ParticleGlobal(isp)%PO(i_part)%X = ParticleGlobal(isp)%PO(i_part)%X + dt*ParticleGlobal(isp)%PO(i_part)%Vx
    !---------- update polar positions
		    ParticleGlobal(isp)%PO(i_part)%Y = DSQRT(X*X+Y*Y)
            ParticleGlobal(isp)%PO(i_part)%Z = DATAN(Y/X)
    !--------- place particle in proper quadrant --------
		    IF (X <= 0.0) THEN
			    ParticleGlobal(isp)%PO(i_part)%Z=ParticleGlobal(isp)%PO(i_part)%Z+PI
		    ENDIF
        ENDIF
    END DO
END DO


!DO i_part = 1, ntot
!    
!    IF (Bfiled_index) THEN   !!!! 有磁场
!        
!        ispe  =  INT(part(i_part,7))
!	
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        xp = (part(i_part,1) - Vert_o(1)) * hxi(1)
!	    i = xp        !!! 粒子在第i个单元里
!	    dx = xp-i
!    
!        yp = (part(i_part,2) - Vert_o(2)) * hxi(2)
!	    j = yp        !!! 粒子在第j个单元里
!	    dy = yp-j
!        
!        !$ ===================== mb.ZWZ 2021/7/11 ======================= \\
!        xcellmdx = 1. -dx
!        ycellmdy = 1. -dy
!        
!        IF(delta_global == 0) THEN
!            P1 = xcellmdx*ycellmdy
!            P2 = dx*ycellmdy
!            P3 = xcellmdx*dy
!            P4 = dx*dy
!        ELSEIF(delta_global == 1)THEN
!            R1=f_left_wall(2) + float(j - 1)*hx(2)
!	        R2=f_left_wall(2) + float(j)*hx(2)
!            R = part(i_part,2)
!	        den=R2*R2-R1*R1
!            
!            P1 = xcellmdx*(R2*R2-R*R)/den
!            P2 = dx*(R2*R2-R*R)/den
!            P3 = xcellmdx*(R*R-R1*R1)/den
!            P4 = dx*(R*R-R1*R1)/den          
!        ENDIF
!        
!        IF (delta_global == 0) THEN
!        
!            xbfield=bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4		
!            ybfield=bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4		
!            zbfield=bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4
!            
!            Omega(1) = xbfield * qm(ispe) * 0.5 * dt
!            Omega(2) = ybfield * qm(ispe) * 0.5 * dt
!            Omega(3) = zbfield * qm(ispe) * 0.5 * dt
!        
!        ELSE IF (delta_global == 1) THEN
!            !$ axisymmetric e field
!            zbfield = bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4  
!            rbfield = bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4   
!            tbfield = bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4  
!            !$ Convert axisymmetric efield to Cartesian efield
!            xbfield = rbfield*DCOS(part(i_part,3)) - tbfield*DSIN(part(i_part,3))
!            ybfield = rbfield*DSIN(part(i_part,3)) + tbfield*DCOS(part(i_part,3))
!            
!            Omega(1) = zbfield * qm(ispe) * 0.5 * dt
!            Omega(2) = xbfield * qm(ispe) * 0.5 * dt
!            Omega(3) = ybfield * qm(ispe) * 0.5 * dt
!        ENDIF         
!        
!        !$ ===================== mb.ZWZ 2021/7/11 ======================= //
!        
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !f = bfx(i,j)+dx*(bfx(i+1,j)-bfx(i,j))
!        !g = bfx(i,j+1)+dx*(bfx(i+1,j+1)-bfx(i,j+1))
!        !Omega(1) = (f + dy*(g-f)) * qm(ispe) * 0.5 * dt
!        !
!        !f = bfy(i,j)+dy*(bfy(i,j+1)-bfy(i,j))
!        !g = bfy(i+1,j)+dy*(bfy(i+1,j+1)-bfy(i+1,j))
!        !Omega(2) = (f + dx*(g-f)) * qm(ispe) * 0.5 * dt
!        !
!        !f = bfz(i,j)+dx*(bfz(i+1,j)-bfz(i,j))
!        !g = bfz(i,j+1)+dx*(bfz(i+1,j+1)-bfz(i,j+1))
!        !Omega(3) = (f + dy*(g-f)) * qm(ispe) * 0.5 * dt 
!    
!        !SVelocity(1) = part(i_part,4) + 0.5 * A_bar_n(1,i_part) * dt + &
!        !               part(i_part,5)*Omega(3) - part(i_part,6)*Omega(2)
!        !SVelocity(2) = part(i_part,5) + 0.5 * A_bar_n(2,i_part) * dt + &
!        !               part(i_part,6)*Omega(1) - part(i_part,4)*Omega(3)
!        !SVelocity(3) = part(i_part,6) + 0.5 * A_bar_n(3,i_part) * dt + &
!        !               part(i_part,4)*Omega(2) - part(i_part,5)*Omega(1)
!        SVelocity(1) = part(i_part,4) + 0.5 * part(i_part,8) * dt + &
!                       part(i_part,5)*Omega(3) - part(i_part,6)*Omega(2)
!        SVelocity(2) = part(i_part,5) + 0.5 * part(i_part,9) * dt + &
!                       part(i_part,6)*Omega(1) - part(i_part,4)*Omega(3)
!        SVelocity(3) = part(i_part,6) + 0.5 * part(i_part,10) * dt + &
!                       part(i_part,4)*Omega(2) - part(i_part,5)*Omega(1)
!        
!        OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
!        TransB(1,1)=1.d0+Omega(1)*Omega(1)
!        TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
!        TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
!        TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
!        TransB(2,2)=1.d0+Omega(2)*Omega(2)
!        TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
!        TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
!        TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
!        TransB(3,3)=1.d0+Omega(3)*Omega(3)
!        TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT
!           
!        part(i_part,4) = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
!        part(i_part,5) = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
!        part(i_part,6) = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)
!        
!        !part(i_part,1) = part(i_part,1) + dt*part(i_part,4)
!        !part(i_part,2) = part(i_part,2) + dt*part(i_part,5)
!        !part(i_part,3) = part(i_part,3) + dt*part(i_part,6)
!        
!    ELSE  !!!! 无磁场
!    
!        !!! BJW:A_bar is a 2D matrix, A_bar(:,1) is x direction and A_bar(:,2) is y direction.
!        part(i_part,4) = part(i_part,4) + 0.5*dt*part(i_part,8)  
!        part(i_part,5) = part(i_part,5) + 0.5*dt*part(i_part,9)
!        part(i_part,6) = part(i_part,6) + 0.5*dt*part(i_part,10)
!    
!        !part(i_part,1) = part(i_part,1) + dt*part(i_part,4)
!        !part(i_part,2) = part(i_part,2) + dt*part(i_part,5)
!        !part(i_part,3) = part(i_part,3) + dt*part(i_part,6)   
!        
!    END IF
!    
!    IF (delta_global == 0) THEN
!        part(i_part,1:2) = part(i_part,1:2) + dt*part(i_part,4:5)
!!Debug, Z is meaningless for 2-D planar problem
!		part(i_part,3) = part(i_part,3) + dt*part(i_part,6)
!    ELSEIF (delta_global == 1) THEN
!    
!!------ transform positions to cartesian coordinates ------
!		X=part(i_part,2)*DCOS(part(i_part,3))
!		Y=part(i_part,2)*DSIN(part(i_part,3))
!		Z=part(i_part,1)
!!------ update cartesian positions ---------
!		X=X+dt*part(i_part,5)
!        IF(X == 0.) THEN
!            X=hx(2)*1.0E-5
!        ENDIF
!		Y=Y+dt*part(i_part,6)
!		part(i_part,1) = part(i_part,1) + dt*part(i_part,4)
!!---------- update polar positions
!		part(i_part,2) = DSQRT(X*X+Y*Y)
!        part(i_part,3) = DATAN(Y/X)
!!--------- place particle in proper quadrant --------
!		IF (X <= 0.0) THEN
!			part(i_part,3)=part(i_part,3)+PI
!		ENDIF
!    ENDIF
!    
!END DO

END