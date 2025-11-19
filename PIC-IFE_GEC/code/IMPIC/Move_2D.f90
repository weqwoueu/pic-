SUBROUTINE Move_2D(it, dt, delta)
    
!!! bjw add 2018-5-22 推进粒子

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Particle_2D
USE Field_2D
!USE IMPIC_Data_2D
USE IMPIC_Data_2D, ONLY: Bfiled_index, Omega
IMPLICIT NONE
REAL(8)		dt
INTEGER		delta, it

INTEGER   :: i, j, i_part, ispe
REAL(8)   :: xp, dx, deta_vx, yp, dy, deta_vy, deta_vz
REAL(8)   :: Efield(3), f, g
REAL(8)   :: SVelocity(1:3)=0.d0

INTEGER         :: index

!$ ========= ab.ZWZ 2021/7/13 ====== \\
REAL(8)   :: xcellmdx, ycellmdy
REAL(8)   :: R, R1, R2, den
REAL(8)   :: P1, P2, P3, P4
REAL(8)	  :: Bfield(3)
REAL(8)   :: X, Y, Z
REAL(8)   :: zefield, refield, tefield
REAL(8)   :: xefield, yefield
REAL(8)   :: zbfield, rbfield, tbfield
REAL(8)   :: xbfield, ybfield
!$ ========= ab.ZWZ 2021/7/13 ====== //




DO i_part = 1, ntot
    
    ispe  =  INT(part(i_part,7))
    
    xp = (part(i_part,1) - Vert_o(1))*hxi(1)
	i = xp        !!! 粒子在第i个单元里
	dx = xp-i
    
    yp = (part(i_part,2) - Vert_o(2))*hxi(2)
	j = yp        !!! 粒子在第i个单元里
	dy = yp-j

    !$ ===================== mb.ZWZ 2021/7/11 ======================= \\
    xcellmdx = 1. -dx
    ycellmdy = 1. -dy
        
    IF(delta_global == 0) THEN
        P1 = xcellmdx *ycellmdy
        P2 = dx       *ycellmdy
        P3 = xcellmdx *dy
        P4 = dx       *dy
    ELSEIF(delta_global == 1)THEN
        R1=f_left_wall(2) + float(j - 1)*hx(2)
	    R2=f_left_wall(2) + float(j)*hx(2)
        R = part(i_part,2)
	    den=R2*R2-R1*R1
            
        P1 = xcellmdx *(R2*R2-R*R)/den
        P2 = dx       *(R2*R2-R*R)/den
        P3 = xcellmdx *(R*R-R1*R1)/den
        P4 = dx       *(R*R-R1*R1)/den 
    ENDIF
    !$ ===================== mb.ZWZ 2021/7/11 ======================= //
    
    IF (Bfiled_index) THEN   !!!! 有磁场
        
        !!!!!!!!!!!!!!!!!!!!!!!! 算法一：同姜巍 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !f = efx(i,j)+dx*(efx(i+1,j)-efx(i,j))
        !g = efx(i,j+1)+dx*(efx(i+1,j+1)-efx(i,j+1))
        !Efield(1) = f + dy * (g - f)
        !deta_vx = qm(ispe)*Efield(1)*dt      !!! 1/2*q*E/m 
        !
        !f = efy(i,j)+dy*(efy(i,j+1)-efy(i,j))
        !g = efy(i+1,j)+dy*(efy(i+1,j+1)-efy(i+1,j))
        !Efield(2) = f + dx * (g - f)
        !deta_vy = qm(ispe)*Efield(2)*dt      !!! 1/2*q*E/m 
        !
        !f = efz(i,j)+dy*(efz(i,j+1)-efz(i,j))
        !g = efz(i+1,j)+dy*(efz(i+1,j+1)-efz(i+1,j))
        !Efield(3) = f + dx * (g - f)
        !deta_vz = qm(ispe)*Efield(3)*dt      !!! 1/2*q*E/m
        !
        !f = bfx(i,j)+dx*(bfx(i+1,j)-bfx(i,j))
        !g = bfx(i,j+1)+dx*(bfx(i+1,j+1)-bfx(i,j+1))
        !Omega(1) = (f + dy*(g-f)) * qm(ispe) * 0.5 * dt
        !
        !f = bfy(i,j)+dy*(bfy(i,j+1)-bfy(i,j))
        !g = bfy(i+1,j)+dy*(bfy(i+1,j+1)-bfy(i+1,j))
        !Omega(2) = (f + dx*(g-f)) * qm(ispe) * 0.5 * dt
        !
        !f = bfz(i,j)+dx*(bfz(i+1,j)-bfz(i,j))
        !g = bfz(i,j+1)+dx*(bfz(i+1,j+1)-bfz(i,j+1))
        !Omega(3) = (f + dy*(g-f)) * qm(ispe) * 0.5 * dt 
        !
        !OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
        !TransB(1,1)=1.d0+Omega(1)*Omega(1)
        !TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
        !TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
        !TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
        !TransB(2,2)=1.d0+Omega(2)*Omega(2)
        !TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
        !TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
        !TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
        !TransB(3,3)=1.d0+Omega(3)*Omega(3)
        !TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT
        !
        !SVelocity(1)=(TransB(1,1)*deta_vx+TransB(1,2)*deta_vy+TransB(1,3)*deta_vz)*0.5d0
        !SVelocity(2)=(TransB(2,1)*deta_vx+TransB(2,2)*deta_vy+TransB(2,3)*deta_vz)*0.5d0
        !SVelocity(3)=(TransB(3,1)*deta_vx+TransB(3,2)*deta_vy+TransB(3,3)*deta_vz)*0.5d0
        !
        !part(i_part,4) = part(i_part,4) + SVelocity(1)
        !part(i_part,5) = part(i_part,5) + SVelocity(2)
        !part(i_part,6) = part(i_part,6) + SVelocity(3)
        !
        !part(i_part,8) = 0.5*(part(i_part,8) + deta_vx)
        !part(i_part,9) = 0.5*(part(i_part,9) + deta_vy)
        !part(i_part,10) = 0.5*(part(i_part,10) + deta_vz)
        !
        !SVelocity(1) = part(i_part,4) + 0.5 * part(i_part,8) * dt + &
        !               part(i_part,5)*Omega(3) - part(i_part,6)*Omega(2)
        !SVelocity(2) = part(i_part,5) + 0.5 * part(i_part,9) * dt + &
        !               part(i_part,6)*Omega(1) - part(i_part,4)*Omega(3)
        !SVelocity(3) = part(i_part,6) + 0.5 * part(i_part,10) * dt + &
        !               part(i_part,4)*Omega(2) - part(i_part,5)*Omega(1)
        !
        !part(i_part,4) = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
        !part(i_part,5) = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
        !part(i_part,6) = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)
        !
        !part(i_part,1) = part(i_part,1) + dt*part(i_part,4)
        !part(i_part,2) = part(i_part,2) + dt*part(i_part,5)
        !part(i_part,3) = part(i_part,3) + dt*part(i_part,6)
        !!!!!!!!!!!!!!!!!!!!!!!! 算法一：同姜巍 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !!!!!!!!!!!!!!!!!!!!!!!! 算法二：同原实验室程序 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !f = efx(i,j)+dx*(efx(i+1,j)-efx(i,j))
        !g = efx(i,j+1)+dx*(efx(i+1,j+1)-efx(i,j+1))
        !Efield(1) = f + dy * (g - f)
        !deta_vx = 0.5*qm(ispe)*Efield(1)*dt      !!! 1/2*q*E/m 
        !
        !f = efy(i,j)+dy*(efy(i,j+1)-efy(i,j))
        !g = efy(i+1,j)+dy*(efy(i+1,j+1)-efy(i+1,j))
        !Efield(2) = f + dx * (g - f)
        !deta_vy = 0.5*qm(ispe)*Efield(2)*dt      !!! 1/2*q*E/m 
        !
        !f = efz(i,j)+dy*(efz(i,j+1)-efz(i,j))
        !g = efz(i+1,j)+dy*(efz(i+1,j+1)-efz(i+1,j))
        !Efield(3) = f + dx * (g - f)
        !deta_vz = 0.5*qm(ispe)*Efield(3)*dt      !!! 1/2*q*E/m
        !
        !f = bfx(i,j)+dx*(bfx(i+1,j)-bfx(i,j))
        !g = bfx(i,j+1)+dx*(bfx(i+1,j+1)-bfx(i,j+1))
        !Omega(1) = (f + dy*(g-f)) * qm(ispe) * 0.5 * dt
        !
        !f = bfy(i,j)+dy*(bfy(i,j+1)-bfy(i,j))
        !g = bfy(i+1,j)+dy*(bfy(i+1,j+1)-bfy(i+1,j))
        !Omega(2) = (f + dx*(g-f)) * qm(ispe) * 0.5 * dt
        !
        !f = bfz(i,j)+dx*(bfz(i+1,j)-bfz(i,j))
        !g = bfz(i,j+1)+dx*(bfz(i+1,j+1)-bfz(i,j+1))
        !Omega(3) = (f + dy*(g-f)) * qm(ispe) * 0.5 * dt 
        
        !$ ========================== mb.ZWZ 2021/7/13 ======================== \\
        IF (delta == 0) THEN
            Efield(1) = efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4   
            Efield(2) = efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4   
            Efield(3) = efz(i,j)*P1 + efz(i+1,j)*P2 + efz(i,j+1)*P3 + efz(i+1,j+1)*P4   
        
            deta_vx = 0.5*qm(ispe)*Efield(1)*dt
            deta_vy = 0.5*qm(ispe)*Efield(2)*dt
            deta_vz = 0.5*qm(ispe)*Efield(3)*dt
        
            Bfield(1) = bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4   
            Bfield(2) = bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4   
            Bfield(3) = bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4
        
            Omega(1) = Bfield(1) * qm(ispe) * 0.5 * dt
            Omega(2) = Bfield(2) * qm(ispe) * 0.5 * dt
            Omega(3) = Bfield(3) * qm(ispe) * 0.5 * dt
        
        ELSE IF (delta == 1) THEN
            !$ axisymmetric e field
            zefield = efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4  
            refield = efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4   
            tefield = efz(i,j)*P1 + efz(i+1,j)*P2 + efz(i,j+1)*P3 + efz(i+1,j+1)*P4  
            !$ Convert axisymmetric efield to Cartesian efield
            xefield = refield*DCOS(part(i_part,3)) - tefield*DSIN(part(i_part,3))
            yefield = refield*DSIN(part(i_part,3)) + tefield*DCOS(part(i_part,3))
            
            deta_vx =  0.5*qm(ispe)*zefield*dt
            deta_vy =  0.5*qm(ispe)*xefield*dt
            deta_vz =  0.5*qm(ispe)*yefield*dt
            
            !$ axisymmetric e field
            zbfield = bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4  
            rbfield = bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4   
            tbfield = bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4  
            !$ Convert axisymmetric efield to Cartesian efield
            xbfield = rbfield*DCOS(part(i_part,3)) - tbfield*DSIN(part(i_part,3))
            ybfield = rbfield*DSIN(part(i_part,3)) + tbfield*DCOS(part(i_part,3))
            
            Omega(1) = zbfield * qm(ispe) * 0.5 * dt
            Omega(2) = xbfield * qm(ispe) * 0.5 * dt
            Omega(3) = ybfield * qm(ispe) * 0.5 * dt
                     
        ENDIF
        !$ ========================== mb.ZWZ 2021/7/13 ======================== //
        
        part(i_part,4) = part(i_part,4) + deta_vx
        part(i_part,5) = part(i_part,5) + deta_vy
        part(i_part,6) = part(i_part,6) + deta_vz
        
        f = 2.0/(1.0 + Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3))
        
        SVelocity(1) = (part(i_part,4) + part(i_part,5)*Omega(3) - part(i_part,6)*Omega(2)) * f
        SVelocity(2) = (part(i_part,5) + part(i_part,6)*Omega(1) - part(i_part,4)*Omega(3)) * f
        SVelocity(3) = (part(i_part,6) + part(i_part,4)*Omega(2) - part(i_part,5)*Omega(1)) * f
        
        part(i_part,4) = part(i_part,4) + SVelocity(2)*Omega(3) - SVelocity(3)*Omega(2) + deta_vx
        part(i_part,5) = part(i_part,5) + SVelocity(3)*Omega(1) - SVelocity(1)*Omega(3) + deta_vy
        part(i_part,6) = part(i_part,6) + SVelocity(1)*Omega(2) - SVelocity(2)*Omega(1) + deta_vz
        
        !part(i_part,1) = part(i_part,1) + dt*part(i_part,4)
        !part(i_part,2) = part(i_part,2) + dt*part(i_part,5)
        !part(i_part,3) = part(i_part,3) + dt*part(i_part,6)
        !!!!!!!!!!!!!!!!!!!!!!!! 算法二：同原实验室程序 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    ELSE !!!! 无磁场
        
        !f = efx(i,j)+dx*(efx(i+1,j)-efx(i,j))
        !g = efx(i,j+1)+dx*(efx(i+1,j+1)-efx(i,j+1))
        !Efield(1) = f + dy * (g - f)
        !deta_vx = qm(ispe)*Efield(1)*dt      !!! 1/2*q*E/m 
        !part(i_part,4) = part(i_part,4) + deta_vx
        !
        !f = efy(i,j)+dy*(efy(i,j+1)-efy(i,j))
        !g = efy(i+1,j)+dy*(efy(i+1,j+1)-efy(i+1,j))
        !Efield(2) = f + dx * (g - f)
        !deta_vy = qm(ispe)*Efield(2)*dt      !!! 1/2*q*E/m 
        !part(i_part,5) = part(i_part,5) + deta_vy
        
        !$ ========================== mb.ZWZ 2021/7/13 ======================== \\
        IF(delta == 0) THEN
            Efield(1) = efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4   
            Efield(2) = efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4     
        
            deta_vx = qm(ispe)*Efield(1)*dt
            deta_vy = qm(ispe)*Efield(2)*dt
        
            part(i_part,4) = part(i_part,4) + deta_vx
            part(i_part,5) = part(i_part,5) + deta_vy
        
        !$ ========================== mb.ZWZ 2021/7/13 ======================== //
        
        !part(i_part,1:3) = part(i_part,1:3) + part(i_part,4:6)*dt
        
        ELSEIF(delta == 1) THEN
            
            !$ axisymmetric e field
            zefield = efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4  
            refield = efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4   
            tefield = efz(i,j)*P1 + efz(i+1,j)*P2 + efz(i,j+1)*P3 + efz(i+1,j+1)*P4   
            !$ Convert axisymmetric efield to Cartesian efield
            xefield = refield*DCOS(part(i_part,3)) - tefield*DSIN(part(i_part,3))
            yefield = refield*DSIN(part(i_part,3)) + tefield*DCOS(part(i_part,3))
        
            deta_vx =  qm(ispe)*zefield*dt
            deta_vy =  qm(ispe)*xefield*dt
            deta_vz =  qm(ispe)*yefield*dt
            
            part(i_part,4) = part(i_part,4) + deta_vx
            part(i_part,5) = part(i_part,5) + deta_vy
            part(i_part,6) = part(i_part,6) + deta_vz
            
            !IF( ABS(part(i_part,5)) > 10)THEN
            !    print*,'speed is too fast'
            !    print*,i_part
            !    print*,part(i_part,:)
            !    pause
            !ENDIF
        ENDIF
    END IF
    
    
    IF (delta == 0) THEN
        part(i_part,1:2) = part(i_part,1:2) + dt*part(i_part,4:5)
!Debug, Z is meaningless for 2-D planar problem
		part(i_part,3) = part(i_part,3) + dt*part(i_part,6)
    ELSEIF (delta == 1) THEN
    
!*****transform positions to cartesian coordinates*****
		X=part(i_part,2)*DCOS(part(i_part,3))
		Y=part(i_part,2)*DSIN(part(i_part,3))
		Z=part(i_part,1)
!*****update cartesian positions*****
		X=X+dt*part(i_part,5)
		IF(X == 0) X=hx(2)*1.0E-5
		Y=Y+dt*part(i_part,6)
		part(i_part,1) = part(i_part,1) + dt*part(i_part,4)
!*****update polar positions
		part(i_part,2) = DSQRT(X*X+Y*Y)
        part(i_part,3) = DATAN(Y/X)
!*****place particle in proper quadrant*****
		IF (X <= 0.0) THEN
			part(i_part,3)=part(i_part,3)+PI
		ENDIF
    ENDIF
END DO
    

    
END
