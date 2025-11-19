MOdule ModuleTrackingParticle
    Use ModuleControlFlow
    Use ModuleParticleBundle
    Implicit none
    
    type Tracking_Particle_Information
        integer(4) :: Particle_Type(2) = (/0,1/)
        integer(4) :: Particle_Number(2) = (/1,10/)
        Real(8) :: X(2)  = (/1.d-5,5.d-3/)
        Real(8) :: Vx(2) = 0.d0
        Real(8) :: Vy(2) = 0.d0
        Real(8) :: Vz(2) = 0.d0
        Real(8) :: Ax(2) = 0.d0
        Real(8) :: Ay(2) = 0.d0
        Real(8) :: Az(2) = 0.d0
    end type Tracking_Particle_information
    
    type Tracking_Particle_Temp
        integer(4) :: Shift_Max
        integer(4) :: Shift
        real(8),allocatable :: Temp_Information(:,:)
        !contains
        !procedure :: TP_Temp_Initial => Tracking_Particle_Temp_Initial
        !procedure :: DTP_Temp => Dump_Tracking_Particle_Temp
    end type Tracking_Particle_Temp
        
    contains
         subroutine Tracking_Particle_Initial(PB,TPI)   
             implicit none
             type(ParticleBundle),intent(inout) :: PB(:)
             type(Tracking_Particle_Information),intent(in) :: TPI
             integer(4) :: i,j
             do i=1,size(TPI%Particle_Type)
                 do j=1,size(TPI%Particle_Number)
                     PB(TPI%Particle_Type(i) + 1)%PO(TPI%Particle_Number(j))%X  =  TPI%X(i) / PB(TPI%Particle_Type(i) + 1)%XFactor
                     PB(TPI%Particle_Type(i) + 1)%PO(TPI%Particle_Number(j))%Vx = TPI%Vx(i) / PB(TPI%Particle_Type(i) + 1)%VFactor
                     PB(TPI%Particle_Type(i) + 1)%PO(TPI%Particle_Number(j))%Vy = TPI%Vy(i) / PB(TPI%Particle_Type(i) + 1)%VFactor
                     PB(TPI%Particle_Type(i) + 1)%PO(TPI%Particle_Number(j))%Vz = TPI%Vz(i) / PB(TPI%Particle_Type(i) + 1)%VFactor
                     PB(TPI%Particle_Type(i) + 1)%PO(TPI%Particle_Number(j))%Ax = TPI%Ax(i)
                     PB(TPI%Particle_Type(i) + 1)%PO(TPI%Particle_Number(j))%Ay = TPI%Ay(i)
                     PB(TPI%Particle_Type(i) + 1)%PO(TPI%Particle_Number(j))%Az = TPI%Az(i)
                 end do
             end do
             Return
         end subroutine Tracking_Particle_Initial
         
         subroutine Tracking_Particle(PB,TPI,CF,ti,tj,TP_Temp)
             implicit none
             integer(4),intent(in) :: ti,tj
             Type(ControlFlow),intent(in) :: CF
             type(ParticleBundle),intent(inout) :: PB(:)
             type(Tracking_Particle_Information),intent(in) :: TPI
             type(Tracking_Particle_Temp),intent(inout) :: TP_Temp 
             real(8) :: ValueTemp(10) = 0
             integer(4) :: i,j,k=1
               do i=1,size(TPI%Particle_Type)
                 do j=1,size(TPI%Particle_Number)
                     ValueTemp(k) = PB(TPI%Particle_Type(i) + 1)%PO(TPI%Particle_Number(j))%X * PB(TPI%Particle_Type(i) + 1)%XFactor
                     k = k + 1       
                 end do
               end do
               k = 1
               TP_Temp%Temp_Information(TP_Temp%shift,1) = Dble(ti + (tj - 1) * CF%Period) * CF%Dt
               do j=2,(size(TPI%Particle_Type) * size(TPI%Particle_Number)+1)
                   TP_Temp%Temp_Information(TP_Temp%Shift,j) = ValueTemp(j-1)
               end do
               if (TP_Temp%shift==TP_Temp%Shift_Max) then
                   call Dump_Tracking_Particle_Temp(TP_Temp)
               else
                   TP_Temp%shift = TP_Temp%Shift + 1
               end if
            return
         end subroutine Tracking_Particle
         

         
         subroutine  Tracking_Particle_Temp_Initial(CF,TPI,TP_Temp)
             implicit none
             Type(ControlFlow),intent(in) :: CF
             Type(Tracking_Particle_Information),intent(in) :: TPI
             type(Tracking_Particle_Temp),intent(inout) :: TP_Temp 
             integer(4) :: Array_Column
             TP_Temp%Shift_Max = 1 * CF%Period
             TP_Temp%Shift = 1
             Array_Column = (size(TPI%Particle_Type)*size(TPI%Particle_Number)+1)
             allocate(TP_Temp%Temp_Information(TP_Temp%Shift_Max,Array_Column))
             TP_Temp%Temp_Information = 1
             return
         end subroutine  Tracking_Particle_Temp_Initial
         
         subroutine Dump_Tracking_Particle_Temp(TP_Temp)
             implicit none
             type(Tracking_Particle_Temp),intent(inout) :: TP_Temp
             integer(4)  i,j,Array_Column
             character(len=50) :: Filename = "Particle_Trajectory_Information.dat"   
             Array_Column =  size(TP_Temp%Temp_Information,2)
             OPEN(Unit=100,file=filename,position='APPEND')
             do i=1,TP_Temp%Shift_Max
                write(100,FMt="(*(es21.14,1x))") (TP_Temp%Temp_Information(i,j),j=1,Array_Column)
             end do
              close(100)
             TP_Temp%Temp_Information = 1
             TP_Temp%Shift = 1
             return
         END subroutine Dump_Tracking_Particle_Temp

         Subroutine Diag_Tracking_Partcile_Initilalization(CF,PB,TPI,TP_Temp)
            Implicit none
            Class(ControlFlow), intent(in) :: CF
            type(ParticleBundle),intent(inout) :: PB(:)
            Type(Tracking_Particle_Information),intent(in) :: TPI
            type(Tracking_Particle_Temp),intent(inout) :: TP_Temp
            call Tracking_Particle_Initial(PB,TPI)
            call Tracking_Particle_Temp_Initial(CF,TPI,TP_Temp)
            return  
        End Subroutine Diag_Tracking_Partcile_Initilalization
 
  End Module ModuleTrackingParticle