SUBROUTINE PostPush
    
USE Domain_2D
USE Particle_2D
USE Field_2D
USE TimeControl
USE IMPIC_Data_2D
USE PIC_MAIN_PARAM_2D !$ ab.ZWZ
Use ModuleMCCInterface,ONLY:ControlFlowGlobal, ParticleGlobal   !$ ab.ZWZ for using JW's particle data structure
IMPLICIT NONE
INTEGER   :: i_part, i, j, ispe
REAL(8)   :: SVelocity(1:3)=0.d0
REAL(8)   :: xp, dx, yp, dy, f, g, A_x, A_y, A_z  
REAL(8)   :: deta_vx, deta_vy, deta_vz

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
    DO i_part = 1, ParticleGlobal(isp)%NPar
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
            xefield=efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4		
            yefield=efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4		
            zefield=efz(i,j)*P1 + efz(i+1,j)*P2 + efz(i,j+1)*P3 + efz(i+1,j+1)*P4
        
            A_x = xefield * qm(isp+1)    
            A_y = yefield * qm(isp+1)
            A_z = zefield * qm(isp+1)
        ELSEIF (delta_global == 1) THEN
            !$ axisymmetric e field
            zefield = efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4  
            refield = efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4   
            tefield = efz(i,j)*P1 + efz(i+1,j)*P2 + efz(i,j+1)*P3 + efz(i+1,j+1)*P4  
            !$ Convert axisymmetric efield to Cartesian efield
            Theta = ParticleGlobal(isp)%PO(i_part)%Z
            xefield = refield*DCOS(Theta) - tefield*DSIN(Theta)
            yefield = refield*DSIN(Theta) + tefield*DCOS(Theta)
            
            A_x = zefield * qm(isp+1)    
            A_y = xefield * qm(isp+1)
            A_z = yefield * qm(isp+1)
        ENDIF

        A_n_plus_1(1) = A_x
        A_n_plus_1(2) = A_y
        A_n_plus_1(3) = A_z
        
        IF (Bfiled_index) THEN   !!!! 有磁场
        
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
        
            SVelocity(1) = A_n_plus_1(1) * 0.5 * dt
            SVelocity(2) = A_n_plus_1(2) * 0.5 * dt
            SVelocity(3) = A_n_plus_1(3) * 0.5 * dt
        
            deta_vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
            deta_vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
            deta_vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)
        
            ParticleGlobal(isp)%PO(i_part)%Vx = ParticleGlobal(isp)%PO(i_part)%Vx + deta_vx
            ParticleGlobal(isp)%PO(i_part)%Vy = ParticleGlobal(isp)%PO(i_part)%Vy + deta_vy
            ParticleGlobal(isp)%PO(i_part)%Vz = ParticleGlobal(isp)%PO(i_part)%Vz + deta_vz
        
            IF(delta_global == 0)THEN
                ParticleGlobal(isp)%PO(i_part)%X = ParticleGlobal(isp)%PO(i_part)%X + dt*deta_vx
                ParticleGlobal(isp)%PO(i_part)%Y = ParticleGlobal(isp)%PO(i_part)%Y + dt*deta_vy
                ParticleGlobal(isp)%PO(i_part)%Z = ParticleGlobal(isp)%PO(i_part)%Z + dt*deta_vz
            ELSEIF(delta_global == 1)THEN
                X=ParticleGlobal(isp)%PO(i_part)%Y*DCOS(ParticleGlobal(isp)%PO(i_part)%Z)
		        Y=ParticleGlobal(isp)%PO(i_part)%Y*DSIN(ParticleGlobal(isp)%PO(i_part)%Z)
		        Z=ParticleGlobal(isp)%PO(i_part)%X
            
                Z = Z + dt*deta_vx
                X = X + dt*deta_vy 
                IF(X == 0.) THEN
                    X=hx(2)*1.0E-5
                ENDIF
                Y = Y + dt*deta_vz
            
                ParticleGlobal(isp)%PO(i_part)%X = Z
                ParticleGlobal(isp)%PO(i_part)%Y = DSQRT(X*X+Y*Y)
                ParticleGlobal(isp)%PO(i_part)%Z = DATAN(Y/X)
                IF (X <= 0.0) THEN
                    ParticleGlobal(isp)%PO(i_part)%Z=ParticleGlobal(isp)%PO(i_part)%Z+PI
                ENDIF
            ENDIF
        
        ELSE  !!!! 无磁场
 
            ParticleGlobal(isp)%PO(i_part)%Vx = ParticleGlobal(isp)%PO(i_part)%Vx + 0.5*dt*A_n_plus_1(1)  
            ParticleGlobal(isp)%PO(i_part)%Vy = ParticleGlobal(isp)%PO(i_part)%Vy + 0.5*dt*A_n_plus_1(2)
            ParticleGlobal(isp)%PO(i_part)%Vz = ParticleGlobal(isp)%PO(i_part)%Vz + 0.5*dt*A_n_plus_1(3)
        
            IF(delta_global == 0)THEN
                ParticleGlobal(isp)%PO(i_part)%X = ParticleGlobal(isp)%PO(i_part)%X + 0.5*dt*dt*A_n_plus_1(1)  
                ParticleGlobal(isp)%PO(i_part)%Y = ParticleGlobal(isp)%PO(i_part)%Y + 0.5*dt*dt*A_n_plus_1(2)
                ParticleGlobal(isp)%PO(i_part)%Z = ParticleGlobal(isp)%PO(i_part)%Z + 0.5*dt*dt*A_n_plus_1(3)
            ELSEIF(delta_global == 1)THEN
                X=ParticleGlobal(isp)%PO(i_part)%Y*DCOS(ParticleGlobal(isp)%PO(i_part)%Z)
		        Y=ParticleGlobal(isp)%PO(i_part)%Y*DSIN(ParticleGlobal(isp)%PO(i_part)%Z)
		        Z=ParticleGlobal(isp)%PO(i_part)%X
            
                Z = Z + 0.5*dt*dt*A_n_plus_1(1)
                X = X + 0.5*dt*dt*A_n_plus_1(2) 
                IF(X == 0.) THEN
                    X=hx(2)*1.0E-5
                ENDIF
                Y = Y + 0.5*dt*dt*A_n_plus_1(3) 
            
                ParticleGlobal(isp)%PO(i_part)%X = Z
                ParticleGlobal(isp)%PO(i_part)%Y = DSQRT(X*X+Y*Y)
                ParticleGlobal(isp)%PO(i_part)%Z = DATAN(Y/X)
                IF (X <= 0.0) THEN
                    ParticleGlobal(isp)%PO(i_part)%Z=ParticleGlobal(isp)%PO(i_part)%Z+PI
                ENDIF
            ENDIF
              
        END IF
    
        A_bar_n_minus_1(1) = ParticleGlobal(isp)%PO(i_part)%Ax   
        A_bar_n_minus_1(2) = ParticleGlobal(isp)%PO(i_part)%Ay   
        A_bar_n_minus_1(3) = ParticleGlobal(isp)%PO(i_part)%Az   
    
        A_bar_n(1) = 0.5*(A_bar_n_minus_1(1) + A_n_plus_1(1)) 
        A_bar_n(2) = 0.5*(A_bar_n_minus_1(2) + A_n_plus_1(2))
        A_bar_n(3) = 0.5*(A_bar_n_minus_1(3) + A_n_plus_1(3))
        
        ParticleGlobal(isp)%PO(i_part)%Ax = A_bar_n(1) 
        ParticleGlobal(isp)%PO(i_part)%Ay = A_bar_n(2)
        ParticleGlobal(isp)%PO(i_part)%Az = A_bar_n(3)
        
    END DO
END DO

!DO i_part = 1, ntot
!    
!    ispe  =  INT(part(i_part,7))
!	
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    xp = (part(i_part,1) - Vert_o(1)) * hxi(1)
!	i = xp        !!! 粒子在第i个单元里
!	dx = xp-i
!    
!    yp = (part(i_part,2) - Vert_o(2)) * hxi(2)
!	j = yp        !!! 粒子在第j个单元里
!	dy = yp-j
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!    !$ ===================== mb.ZWZ 2021/7/11 ======================= \\
!    xcellmdx = 1. -dx
!    ycellmdy = 1. -dy
!        
!    IF(delta_global == 0) THEN
!        P1 = xcellmdx*ycellmdy
!        P2 = dx*ycellmdy
!        P3 = xcellmdx*dy
!        P4 = dx*dy
!    ELSEIF(delta_global == 1)THEN
!        R1=f_left_wall(2) + float(j - 1)*hx(2)
!	    R2=f_left_wall(2) + float(j)*hx(2)
!        R = part(i_part,2)
!	    den=R2*R2-R1*R1
!            
!        P1 = xcellmdx*(R2*R2-R*R)/den
!        P2 = dx*(R2*R2-R*R)/den
!        P3 = xcellmdx*(R*R-R1*R1)/den
!        P4 = dx*(R*R-R1*R1)/den          
!    ENDIF
!    
!    IF (delta_global == 0) THEN
!        xefield=efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4		
!        yefield=efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4		
!        zefield=efz(i,j)*P1 + efz(i+1,j)*P2 + efz(i,j+1)*P3 + efz(i+1,j+1)*P4
!        
!        A_x = xefield * qm(ispe)    
!        A_y = yefield * qm(ispe)
!        A_z = zefield * qm(ispe)
!    ELSEIF (delta_global == 1) THEN
!        !$ axisymmetric e field
!        zefield = efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4  
!        refield = efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4   
!        tefield = efz(i,j)*P1 + efz(i+1,j)*P2 + efz(i,j+1)*P3 + efz(i+1,j+1)*P4  
!        !$ Convert axisymmetric efield to Cartesian efield
!        xefield = refield*DCOS(part(i_part,3)) - tefield*DSIN(part(i_part,3))
!        yefield = refield*DSIN(part(i_part,3)) + tefield*DCOS(part(i_part,3))
!        
!        A_x = zefield * qm(ispe)    
!        A_y = xefield * qm(ispe)
!        A_z = yefield * qm(ispe)
!    ENDIF
!    !$ ===================== mb.ZWZ 2021/7/11 ======================= //
!    
!    !f = efx(i,j)+dx*(efx(i+1,j)-efx(i,j))
!    !g = efx(i,j+1)+dx*(efx(i+1,j+1)-efx(i,j+1))
!    !A_x = (f + dy*(g-f)) * qm(ispe)
!    !
!    !f = efy(i,j)+dy*(efy(i,j+1)-efy(i,j))
!    !g = efy(i+1,j)+dy*(efy(i+1,j+1)-efy(i+1,j))
!    !A_y = (f + dx*(g-f)) * qm(ispe)
!    !
!    !f = efz(i,j)+dy*(efz(i,j+1)-efz(i,j))
!    !g = efz(i+1,j)+dy*(efz(i+1,j+1)-efz(i+1,j))
!    !A_z = (f + dx*(g-f)) * qm(ispe)
!    
!    !A_n_plus_1(1,i_part) = A_x
!    !A_n_plus_1(2,i_part) = A_y
!    !A_n_plus_1(3,i_part) = A_z
!    A_n_plus_1(1) = A_x
!    A_n_plus_1(2) = A_y
!    A_n_plus_1(3) = A_z
!        
!    IF (Bfiled_index) THEN   !!!! 有磁场
!        
!        !$ ===================== mb.ZWZ 2021/7/11 ======================= \\      
!        IF (delta_global == 0) THEN
!        
!            xbfield=bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4		
!            ybfield=bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4		
!            zbfield=bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4
!        
!            Omega(1) = xbfield * qm(ispe) * 0.5 * dt
!            Omega(2) = ybfield * qm(ispe) * 0.5 * dt
!            Omega(3) = zbfield * qm(ispe) * 0.5 * dt 
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
!        !$ ===================== mb.ZWZ 2021/7/11 ======================= //
!        
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
!        !SVelocity(1) = A_n_plus_1(1,i_part) * 0.5 * dt
!        !SVelocity(2) = A_n_plus_1(2,i_part) * 0.5 * dt
!        !SVelocity(3) = A_n_plus_1(3,i_part) * 0.5 * dt
!        SVelocity(1) = A_n_plus_1(1) * 0.5 * dt
!        SVelocity(2) = A_n_plus_1(2) * 0.5 * dt
!        SVelocity(3) = A_n_plus_1(3) * 0.5 * dt
!        
!        deta_vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
!        deta_vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
!        deta_vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)
!        
!        part(i_part,4) = part(i_part,4) + deta_vx
!        part(i_part,5) = part(i_part,5) + deta_vy
!        part(i_part,6) = part(i_part,6) + deta_vz
!        
!        !part(i_part,1) = part(i_part,1) + dt*deta_vx
!        !part(i_part,2) = part(i_part,2) + dt*deta_vy
!        !part(i_part,3) = part(i_part,3) + dt*deta_vz
!        
!        IF(delta_global == 0)THEN
!            part(i_part,1) = part(i_part,1) + dt*deta_vx
!            part(i_part,2) = part(i_part,2) + dt*deta_vy
!            !Debug, Z is meaningless for 2-D planar problem
!            part(i_part,3) = part(i_part,3) + dt*deta_vz
!        ELSEIF(delta_global == 1)THEN
!            X=part(i_part,2)*DCOS(part(i_part,3))
!		    Y=part(i_part,2)*DSIN(part(i_part,3))
!		    Z=part(i_part,1)
!            
!            Z = Z + dt*deta_vx
!            X = X + dt*deta_vy 
!            IF(X == 0.) THEN
!                X=hx(2)*1.0E-5
!            ENDIF
!            Y = Y + dt*deta_vz
!            
!            part(i_part,1) = Z
!            part(i_part,2) = DSQRT(X*X+Y*Y)
!            part(i_part,3) = DATAN(Y/X)
!            IF (X <= 0.0) THEN
!                part(i_part,3)=part(i_part,3)+PI
!            ENDIF
!        ENDIF
!        
!    ELSE  !!!! 无磁场
! 
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST 1 (Wei Jiang)
!        !part(i_part,1) = part(i_part,1) + 0.5*dt*dt*A_n_plus_1(1,i_part)  ! BJW:A_x,A_y,A_z are the accelerations for different directions.
!        !part(i_part,2) = part(i_part,2) + 0.5*dt*dt*A_n_plus_1(2,i_part)
!        !!Debug, Z is meaningless for 2-D planar problem
!        !part(i_part,3) = part(i_part,3) + 0.5*dt*dt*A_n_plus_1(3,i_part)
!        !
!        !part(i_part,4) = part(i_part,4) + 0.5*dt*A_n_plus_1(1,i_part)  
!        !part(i_part,5) = part(i_part,5) + 0.5*dt*A_n_plus_1(2,i_part)
!        !!Debug, Z is meaningless for 2-D planar problem
!        !part(i_part,6) = part(i_part,6) + 0.5*dt*A_n_plus_1(3,i_part)
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TEST 1 (Wei Jiang)
!        !part(i_part,1) = part(i_part,1) + 0.5*dt*dt*A_n_plus_1(1)  ! BJW:A_x,A_y,A_z are the accelerations for different directions.
!        !part(i_part,2) = part(i_part,2) + 0.5*dt*dt*A_n_plus_1(2)
!        !!Debug, Z is meaningless for 2-D planar problem
!        !part(i_part,3) = part(i_part,3) + 0.5*dt*dt*A_n_plus_1(3)
!    
!        part(i_part,4) = part(i_part,4) + 0.5*dt*A_n_plus_1(1)  
!        part(i_part,5) = part(i_part,5) + 0.5*dt*A_n_plus_1(2)
!        !Debug, Z is meaningless for 2-D planar problem
!        part(i_part,6) = part(i_part,6) + 0.5*dt*A_n_plus_1(3)
!        
!        IF(delta_global == 0)THEN
!            part(i_part,1) = part(i_part,1) + 0.5*dt*dt*A_n_plus_1(1)  
!            part(i_part,2) = part(i_part,2) + 0.5*dt*dt*A_n_plus_1(2)
!            !Debug, Z is meaningless for 2-D planar problem
!            part(i_part,3) = part(i_part,3) + 0.5*dt*dt*A_n_plus_1(3)
!        ELSEIF(delta_global == 1)THEN
!            X=part(i_part,2)*DCOS(part(i_part,3))
!		    Y=part(i_part,2)*DSIN(part(i_part,3))
!		    Z=part(i_part,1)
!            
!            Z = Z + 0.5*dt*dt*A_n_plus_1(1)
!            X = X + 0.5*dt*dt*A_n_plus_1(2) 
!            IF(X == 0.) THEN
!                X=hx(2)*1.0E-5
!            ENDIF
!            Y = Y + 0.5*dt*dt*A_n_plus_1(3) 
!            
!            part(i_part,1) = Z
!            part(i_part,2) = DSQRT(X*X+Y*Y)
!            part(i_part,3) = DATAN(Y/X)
!            IF (X <= 0.0) THEN
!                part(i_part,3)=part(i_part,3)+PI
!            ENDIF
!        ENDIF
!              
!    END IF
!    
!    !A_bar_n_minus_1(i_part,:) = A_bar_n(i_part,:) 
!    !A_bar_n_minus_1(:,i_part) = part(i_part,8:10)   !!! bjw add 2018-3-5
!    A_bar_n_minus_1(:) = part(i_part,8:10)   !$ mb.ZWZ
!    
!    !A_bar_n(1,i_part) = 0.5*(A_bar_n_minus_1(1,i_part) + A_n_plus_1(1,i_part)) ! BJW:A_bar is a 2D matrix, A_bar(:,1) is x direction and A_bar(:,2) is y direction.
!    !A_bar_n(2,i_part) = 0.5*(A_bar_n_minus_1(2,i_part) + A_n_plus_1(2,i_part))
!    !A_bar_n(3,i_part) = 0.5*(A_bar_n_minus_1(3,i_part) + A_n_plus_1(3,i_part))
!    A_bar_n(1) = 0.5*(A_bar_n_minus_1(1) + A_n_plus_1(1)) ! BJW:A_bar is a 2D matrix, A_bar(:,1) is x direction and A_bar(:,2) is y direction.
!    A_bar_n(2) = 0.5*(A_bar_n_minus_1(2) + A_n_plus_1(2))
!    A_bar_n(3) = 0.5*(A_bar_n_minus_1(3) + A_n_plus_1(3))
!        
!    !!!****************** bjw add for test bug 2018-3-5 *****************
!    !part(i_part,8) = A_bar_n(1,i_part) 
!    !part(i_part,9) = A_bar_n(2,i_part)
!    !!Debug, Z is meaningless for 2-D planar problem
!    !part(i_part,10) = A_bar_n(3,i_part)
!    part(i_part,8) = A_bar_n(1) 
!    part(i_part,9) = A_bar_n(2)
!    !Debug, Z is meaningless for 2-D planar problem
!    part(i_part,10) = A_bar_n(3)
!    !!!******************************************************************
!
!        
!END DO


END


