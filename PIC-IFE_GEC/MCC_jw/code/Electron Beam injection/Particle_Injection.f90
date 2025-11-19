Module Particle_Injection
    use ModuleParticleBundle
    use ModuleControlFlow
    use Constants
    implicit none
        integer(4),parameter :: Inject_Particle_type = 0
        integer(4),parameter :: Inject_Direction_Mode = 1001   !1001--X,1002--Y,1114--RandoM
        integer(4),parameter :: Inject_Source_Mode = 1         !1_single,2-Double
        real(8),parameter :: Position = 0.1                     !Relative position of the first or last grid
        real(8),parameter :: Derta_P = 1.d-3                   
        real(8),parameter :: Particle_energy = 10               !eV
        real(8),parameter :: Emission_Current = 1.0            ! A 
        
    contains  
    subroutine Particle_Beam_InjectionS(CF,PB)
       Implicit none
       TYPE(ParticleBundle),intent(inout) ::  PB
       Type(ControlFlow),intent(in):: CF
       type(ParticleOne) :: PO(2)
       integer :: i,Number,Emission_Particle_Count = 0
       real(8) :: V,Xtemp,XFactor,VFactor,number_temp 
       XFactor = 1.0 / PB%XFactor
       VFactor = 1.0 / PB%VFactor
       
       number_temp = (CF%Dt * Emission_Current) / (ElectronCharge * PB%weight * CF%Dx)
       Number = int(number_temp )
       call RANDOM_NUMBER(R)
       if (R < (number_temp - Number)) then
           Number = Number + 1
       end if

       Xtemp = Position * CF%Dx  
       !Xtemp = Position * ZLength 
       V = sqrt(((Particle_energy * ElectronCharge) + (Particle_energy * ElectronCharge)) / ElectronMass)

       PO(1)%Ax=0.d0
       PO(1)%Ay=0.d0
       PO(1)%Az=0.d0
       do i=1,Number
           call RANDOM_NUMBER(R)
           PO(1)%X = Xtemp +  CF%Dx * Derta_P * R
           Select case(Inject_Direction_Mode)
               case(1001)
                   PO(1)%Vx = V
                   PO(1)%Vy = 0
                   PO(1)%Vy = 0
               case(1010)
                   PO(1)%Vx = 0
                   PO(1)%Vy = V
                   PO(1)%Vy = 0
               case(1111)
                   PO(1)%Vx = V
                   PO(1)%Vx = V
                   PO(1)%Vx = V
               case(1114)
                   Call VelocityRandomInitializationParticleOne(PO(1),V)
               case(1115)
           end select
           call PositionRescaleParticleOne(PO(1),XFactor)
           call VelocityRescaleParticleOne(PO(1),VFactor)
           PB%NPar = PB%NPar + 1
           PB%PO(PB%NPar) = PO(1)
           Emission_Particle_Count =  Emission_Particle_Count + 1
           if (Inject_Source_Mode == 2) then
               PO(2) = PO(1)
               call RANDOM_NUMBER(R)
               PO(2)%X = ZLength - (Xtemp + CF%Dx * Derta_P * R)
               PO(2)%Vx = (-PO(2)%Vx)
               PO(2)%Vy = (-PO(2)%Vy)
               PO(2)%Vy = (-PO(2)%Vy)
               call PositionRescaleParticleOne(PO(2),XFactor)
               PB%NPar = PB%NPar + 1
               PB%PO(PB%NPar) = PO(2)
               Emission_Particle_Count =  Emission_Particle_Count + 1
           end if
       end do
       return
    END  subroutine Particle_Beam_Injections
    
    end Module Particle_Injection