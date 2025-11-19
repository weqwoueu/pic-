SUBROUTINE InjectBeams_2D(it, delta, PB)
    
USE Domain_2D
USE Particle_2D
USE IMPIC_Data_2D
USE Constant_Variable_2D
USE TimeControl
USE PIC_MAIN_PARAM_2D, only : zero
USE Field_2D
Use ModuleParticleBundle ! wsy add
IMPLICIT NONE
! ------------ input------------
INTEGER		it, delta
Type(ParticleBundle), intent(inout) :: PB(:)
!-------------------------------


!------------- Temporary variables ----------------
Type(ParticleOneIndex), Allocatable :: InjectedParticles(:)
REAL(8)       ::  b_amb(2,2)
INTEGER       ::  jj, ntotp, i_part, i, AddPar, ParticleIndex
INTEGER       ::  add_N(ispe_inject), n_element_old
!REAL(8)       ::  ranum
double precision ::  ranum, R_ranum, Theta_ranum, Z_ranum
REAL(8)       ::  v_r, FS1, FS2, QA, UU, UN, GG, V_x, V_y, V_z, XL, XU, YL, YU, ZL, ZU
CHARACTER*25	:: fname, sname, filename
REAL(8)       ::  r1,r2,x, y

REAL(8) :: Ic, part_x, part_y
REAL(8) :: cathode_left_location, cathode_right_location, cathode_radius_location
! --------------------------------------------------

!------test----
REAL(8)       :: Avg_Vx, Avg_Vy
!--------
!Avg_Vx = 0
!Avg_Vy = 0
n_element_old = 0

DO jj=1,ispe_inject 
    
    IF (vd(ipf(jj),1) > 0.0) THEN
        !!! Average velocity of particles in the direction of incidence£¨v > 0£©
        b_amb(ipf(jj),1)=vt(ipf(jj),1)/SQRT(2.*PI)*DEXP(-(vd(ipf(jj),1))**2/(2*vt(ipf(jj),1)**2))+ &
                         vd(ipf(jj),1)/2*(1+ERF(vd(ipf(jj),1)/(SQRT(2.)*vt(ipf(jj),1))))
    ELSE
        !!!! Average velocity of particles in the direction of incidence£¨v < 0£©
        b_amb(ipf(jj),1)=-vt(ipf(jj),1)/SQRT(2.*PI)*DEXP(-(vd(ipf(jj),1))**2/(2*vt(ipf(jj),1)**2))+ &
                         vd(ipf(jj),1)/2*(1-ERF(vd(ipf(jj),1)/(SQRT(2.)*vt(ipf(jj),1))))
    END IF
    
    CALL DRandom(ranum)

    IF (delta == 0) THEN
        
        add_N(ipf(jj))=int(dens0(ipf(jj))*n_ref*ABS(b_amb(ipf(jj),1))*v_ref*dt*time_ref* &
                          (f_right_wall(2)-f_left_wall(2))*L_ref/affp_bjw(ipf(jj))+ranum) 
        !add_N(ipf(jj))=PB(jj)%nLoss(1)   
    ELSEIF (delta == 1) THEN
        add_N(ipf(jj))=int(dens0(ipf(jj))*n_ref*ABS(b_amb(ipf(jj),1))*v_ref*dt*time_ref* &
                          PI*(f_right_wall(2)-f_left_wall(2))**2*L_ref**2/affp_bjw(ipf(jj))+ranum)
    ENDIF
    
END DO

Allocate(InjectedParticles(sum(add_N)))
!add_N(1) = 0
AddPar = 0

 Do jj = 1, ispe_inject
    
    If (ipf(jj) == 1) Then ! electron
        XL = Min(R_start(ipf(jj)), R_start(ipf(jj)) + 0.5 * b_amb(ipf(jj),1) * dt)
        XU = Max(R_start(ipf(jj)), R_start(ipf(jj)) + 0.5 * b_amb(ipf(jj),1) * dt)
        YL = f_left_wall(2)
        YU = f_right_wall(2)
    Else ! ion
        XL = Min(R_start(ipf(jj)), R_start(ipf(jj)) + 0.5 * b_amb(ipf(jj),1) * dt)
        XU = Max(R_start(ipf(jj)), R_start(ipf(jj)) + 0.5 * b_amb(ipf(jj),1) * dt)
        YL = f_left_wall(2)
        YU = f_right_wall(2)
    End If
    
    
    Do i_part = 1, add_N(ipf(jj))
    
        AddPar = AddPar + 1
        
        !---------------- Velocity Sampling Method 1 --------
            !! In G.A.Bird Book (12.5£©
	        v_r=abs(vd(ipf(jj),1))/(sqrt(2.)*vt(ipf(jj),1))
		    FS1=v_r+SQRT(v_r*v_r+2.)
		    FS2=0.5*(1.+v_r*(2.*v_r-FS1))
            QA=3.			 
            IF (v_r.LT.-3.) QA=ABS(v_r)+1.
            103  CALL DRandom(ranum)	
		    UU=-QA+2.*QA*ranum
		    UN=UU+v_r
            !!!*--UN is a potential inward velocity component
		    IF (UN.LT.0.) GO TO 103
		    GG=(2.*UN/FS1)*DEXP(FS2-UU*UU)
            CALL DRandom(ranum)
		    IF (GG.LT.ranum) GO TO 103
            !!!*--the inward normalised vel. component has been selected (eqn (12.5))

            V_x=UN*(sqrt(2.)*vt(ipf(jj),1))   !!!! wsy
        ! ---------------------------------------------------
        
        
        
        If (ipf(jj) == 1) Then ! electron
            
        
            !---------------- Velocity Sampling Method 2 --------
            !! In G.A.Bird Book Appendix C(P425£©
            CALL Loadv(V_y, tmpj(ipf(jj)), jj)
            CALL Loadv(V_z, tmpj(ipf(jj)), jj)
            ! ---------------------------------------------------
        Else    ! ion
            
            CALL Loadv(V_y, tmpj(ipf(jj)), jj)
            CALL Loadv(V_z, tmpj(ipf(jj)), jj)
            
            !V_x = 0
            V_y = vd(ipf(jj),2) + V_y
            V_z = vd(ipf(jj),3) + V_z
        End If
        
        
        !Avg_Vx = Avg_Vx + V_x
        !Avg_Vy = Avg_Vy + V_y
        Call InjectedParticles(AddPar)%IndexInit(ipf(jj))
        Call InjectedParticles(AddPar)%ParticleOne%PosInit(XL, XU, YL, YU)
        Call InjectedParticles(AddPar)%ParticleOne%VelInpInit(V_x, V_y, V_z)
        Call InjectedParticles(AddPar)%ParticleOne%AccInpInit( zero, zero, zero)
        
        part_x = InjectedParticles(AddPar)%ParticleOne%X
        part_y = InjectedParticles(AddPar)%ParticleOne%Y
        
        Call InjectPositioning(part_x, part_y, n_element_old)
        
        Call InjectedParticles(AddPar)%ParticleOne%ParLocate(n_element_old)
         
    End Do
    !Avg_Vx = Avg_Vx / add_N(ipf(1))
    !Avg_Vy = Avg_Vy / add_N(ipf(1))
End Do


Do jj = 1, Sum(add_N)
    
    ParticleIndex = InjectedParticles(jj)%Index
    Call PB(ParticleIndex)%AddOne(InjectedParticles(jj)%ParticleOne)

End Do


Deallocate(InjectedParticles)


END