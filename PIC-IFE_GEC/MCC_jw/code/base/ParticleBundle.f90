!to Bale different particelone to one sorts
    Module ModuleParticleBundle
      Use ModuleParticleOne
      Use ModuleSpecyOne
      
     Implicit none
     Type :: ParticleBundle
           Integer(4) :: NPar=0,NParNormal=1000    !NPar resprent that the number of particle,NParNorma is how much memory you creat for paritcle
           Real(8) :: XFactor,VFactor              !different parmeter to calculate Coordinate transformation
           Logical :: LXScaled=.False.,LVScaled=.False. !to judge if its be normalized
           Real(8) :: Charge,Mass,Weight
           Real(8) :: Dx,Dt
           Logical :: UnequalWeightFlag = .FALSE.   !The variable denote that whether the weight is equal on every particle in Axis-Symmetric coordinates. 
           Real(8) :: WQOne !> ab.ZWZ for unequal weight
            !> ab.ZWZ for checking particle generation and loss
           Integer(8) :: nGen(0:3)=0    !> 0: total gen | 1: SEE on left boundary| 2: SEE on right boundary | 3: MCC 
           Integer :: nLoss(0:5)=0   !> 0: total gen | 1-4: intersection with each boundary | 5: MCC annihilation 
           Type(SpecyOne),Pointer :: SO    !To avoid be copy of a sort of particale information
           Type(ParticleOne),Allocatable :: PO(:)    !different particle but their are one sorts
           Contains
           Procedure :: AllInit=>AllInitializationParticleBundle
           Procedure :: AddOne=>AddParticleOneParticleBundle
           Procedure :: DelOne=>DelParticleOneParticleBundle
           !Procedure :: NumRes2=>ParticleBundleNumberRescale2
           !Convert in different unit systems fo different sorts particle
           Procedure :: PosRes=>PositionRescaleParticleBundle
           Procedure :: VelRes=>VelocityRescaleParticleBundle
           
           Procedure :: MoveES=>MoveElectrostaticParticleBundle		 
			   procedure :: MoveES1st=>MoveElectrostaticParticleBundle1st
			   procedure :: MoveES2nd=>MoveElectrostaticParticleBundle2nd			
			   
			  Procedure :: ExMoveES=>ExMoveElectrostaticParticleBundle
           Procedure :: MoveEM=>MoveElectromagneticParticleBundle
           Procedure :: WeightP2C=>WeightP2CParticleBundle
			  Procedure :: WeightP2C2=>WeightP2CParticleBundle2
           
           Procedure :: Dump=>DumpParticleBundle  !save partical's information to File
           Procedure :: Load=>LoadParticleBundle  !read partical's information from File
           !to ecpand or contact the particle array
           Procedure :: Norm=>ParticleBundleNormalization
           !Procedure :: RescaleFieldOne=>RescaleFieldOneParticleBundle
        End Type ParticleBundle
           
        Type :: ParticleBundleIndex
           Integer(4) :: NPar=0,NParNormal=100000
           Type(ParticleOneIndex),Allocatable :: POI(:)
           Contains
           Procedure :: Init=>InitializationParticleBundleIndex
           Procedure :: AddOne=>AddParticleOneParticleBundleIndex
           Procedure :: DelOne=>DelParticleOneParticleBundleIndex
           !Procedure :: NumRes2=>ParticleBundleNumberRescale2
           !Procedure :: Dump=>DumpParticleBundleIndex
         End Type ParticleBundleIndex

    contains
!DIR$ ATTRIBUTES FORCEINLINE ::  MoveElectrostaticParticleBundle,  MoveElectromagneticParticleBundle, WeightP2CParticleBundle      
             Subroutine AllInitializationParticleBundle(PB,SO,CF)
		        Class(ParticleBundle), intent(inout) :: PB
                Type(SpecyOne), intent(in),Target :: SO 
                Class(ControlFlow), intent(in) :: CF
                Integer(4) :: i,NPArMax
                Real(8) :: XFactor,VFactor
                Logical ::  Status
                PB%dx = CF%dx
                PB%dt = CF%dt
                PB%XFactor = PB%dx
                PB%VFactor = PB%dx / PB%dt
                
                PB%SO => SO     !to point at SO you have inputted
                PB%Charge = PB%SO%Charge
                PB%Mass = PB%SO%Mass
                PB%Weight = CF%InitDensity / Dble(CF%ParticlePerGrid)

                PB%NParNormal = CF%ParticlePerGrid * (CF%Nx-1) * PB%SO%InitDensity
                PB%NPar = PB%NParNormal           
                NPArMax = Ceiling(2.0 * PB%NParNormal)
                If(Allocated(PB%PO)) Deallocate(PB%PO)
                Allocate(PB%PO(NPArMax))
                If (CF%ReStartParticles) Then
                    Do i=1,PB%Npar
                         PB%LXScaled = .True.
                         Call PB%PO(i)%PosInit(Dble(CF%NxL),Dble(CF%NxU))
                         PB%LVScaled = .True.
                         Call PB%PO(i)%VelMaxInit(PB%SO%Mass,PB%SO%InitTemperature)
                         VFactor = 1.d0 / PB%VFactor
                         Call PB%PO(i)%VelRes(VFactor)
                         Call PB%PO(i)%AccInpInit()
                    End Do
                 Else
                    Call LoadParticleBundle(PB,Status)
                         If (Status) Then
                            Do i=1,PB%NPar
                                PB%LXScaled = .True.
                                Call PB%PO(i)%PosInit(Dble(CF%NxL),Dble(CF%NxU))
                                PB%LVScaled = .True.
                                Call PB%PO(i)%VelMaxInit(PB%SO%Mass,PB%SO%InitTemperature)
                                VFactor = 1.d0 / PB%VFactor
                                Call PB%PO(i)%VelRes(VFactor)
                                Call PB%PO(i)%AccInpInit()
                             End Do
                        Else
                            XFactor = 1.d0 / PB%XFactor
                            Call PB%PosRes()
                            VFactor = 1.d0 / PB%VFactor
                            Call PB%VelRes()
                        End If
                 End If         
             End Subroutine AllInitializationParticleBundle
             
            Subroutine  AddParticleOneParticleBundle(PB,PO)
             Implicit none
                 Class(ParticleBundle),intent(inout) ::  PB
                 Class(ParticleOne),intent(in) ::  PO
                  PB%NPar=PB%NPar+1
                  PB%PO(PB%NPar)=PO
             return   
             End  Subroutine  AddParticleOneParticleBundle
     
             Subroutine  DelParticleOneParticleBundle(PB,NDel)
             Implicit none
                 Class(ParticleBundle),intent(inout) ::  PB
                 Integer(4),intent(in) ::  NDel
                 PB%PO(NDel)=PB%PO(PB%NPar)  !MOVE the last particleone to the No.NDel one to replace the No.Ndel
                 PB%NPar=PB%NPar-1           !
             return   
             End  Subroutine  DelParticleOneParticleBundle
             
             subroutine PositionRescaleParticleBundle(PB)
		       Class(ParticleBundle), intent(inout) :: PB
               Integer(4) :: i
               Real(8) :: XFactor
               If (PB%LXScaled) Then
                   XFactor=PB%XFactor
                   PB%LXScaled=.False.
               Else
                   XFactor=1.d0/PB%XFactor
                   PB%LXScaled=.True.
                End If
                Do i=1,PB%NPar
                      Call PB%PO(i)%PosRes(XFactor)
                End Do
             end subroutine PositionRescaleParticleBundle
             
              subroutine VelocityRescaleParticleBundle(PB)
		       Class(ParticleBundle), intent(inout) :: PB
               Integer(4) :: i
               Real(8) :: VFactor
               If (PB%LVScaled) Then
                   VFactor=PB%VFactor
                   PB%LVScaled=.False.
               Else
                   VFactor=1.d0/PB%VFactor
                   PB%LVScaled=.True.
                End If
                Do i=1,PB%NPar
                      Call PB%PO(i)%VelRes(VFactor)
                ENd DO
              end subroutine VelocityRescaleParticleBundle
              
            subroutine MoveElectrostaticParticleBundle(PB,FG)
		       Class(ParticleBundle), intent(inout) :: PB
               Class(Field), intent(inout) :: FG
               Integer(4) :: i
               Real(8) :: Ex,Ey,Ez,EFactor
               EFactor=PB%Charge/PB%Mass*PB%dt/(PB%dx/PB%dt)
                Do i=1,PB%NPar
                      Call PB%PO(i)%WeightC2PES(FG,Ex)
                      Ex=Ex*EFactor
                      Call PB%PO(i)%MoveES(Ex)
                End Do
				End subroutine MoveElectrostaticParticleBundle  
				
			subroutine MoveElectrostaticParticleBundle1st(PB)
		       Class(ParticleBundle), intent(inout) :: PB
               Integer(4) :: i
                Do i=1,PB%NPar
                      Call PB%PO(i)%MoveES1st()
                End Do
			End subroutine MoveElectrostaticParticleBundle1st
			
			
			  subroutine MoveElectrostaticParticleBundle2nd(PB,FG)
		       Class(ParticleBundle), intent(inout) :: PB
               Class(Field), intent(inout) :: FG
               Integer(4) :: i
               Real(8) :: Ex,Ey,Ez,EFactor
					EFactor=PB%Charge/PB%Mass*PB%dt/(PB%dx/PB%dt)
                Do i=1,PB%NPar
                      Call PB%PO(i)%WeightC2PES(FG,Ex)
                      Ex=Ex*EFactor
                      Call PB%PO(i)%MoveES2nd(Ex)
					 End Do
					 return
			  End subroutine MoveElectrostaticParticleBundle2nd  

				subroutine ExMoveElectrostaticParticleBundle(PB,FG)
				   IMPLICIT NONE
		          Class(ParticleBundle), intent(inout) :: PB
               Class(Field), intent(inout) :: FG
					!Type(ParticleBoundaryOne),intent(in) :: PBDO
               Integer(4) :: i
               Real(8) :: Ex,Ey,Ez,EFactor,EnergyFactor
               EFactor=PB%Charge/PB%Mass*PB%dt/(PB%dx/PB%dt)
                Do i=1,PB%NPar
                      Call PB%PO(i)%WeightC2PES(FG,Ex)
                      Ex=Ex*EFactor
                      Call PB%PO(i)%ExMoveES(Ex)
					 End Do
				End subroutine ExMoveElectrostaticParticleBundle  
			

            Subroutine MoveElectromagneticParticleBundle(PB,FG)
		       Class(ParticleBundle), intent(inout) :: PB
               Type(Field), intent(in) :: FG
               Integer(4) :: i
               Real(8) :: Ex,Ey,Ez,Bx,By,Bz,EFactor,BFactor
               EFactor=PB%Charge/PB%Mass*PB%dt/(PB%dx/PB%dt)
               BFactor=PB%Charge/PB%Mass*PB%dt
                Do i=1,PB%NPar
                      Call PB%PO(i)%WeightC2PEM(FG,Ex,Ey,Ez,Bx,By,Bz)
                      Ex=Ex*EFactor
                      Ey=Ey*EFactor
                      Ez=Ez*EFactor
                      Bx=Bx*BFactor
                      By=By*BFactor
                      Bz=Bz*BFactor
                      Call PB%PO(i)%MoveEM(Ex,Ey,Ez,Bx,By,Bz)
                End Do
            End subroutine MoveElectromagneticParticleBundle
            
            subroutine WeightP2CParticleBundle(PB,FO)
		       Class(ParticleBundle), intent(inout) :: PB
               Type(FieldOne), intent(inout) :: FO
               Integer(4) :: i
               Real(8) :: RhoFactor,ChiFactor
               FO%RhoOne=0.d0
               FO%ChiOne=0.d0
                Do i=1,PB%NPar
                      Call PB%PO(i)%WeightP2C(FO)
                End Do
                RhoFactor=PB%Charge*PB%Weight
                FO%RhoOne=FO%RhoOne*RhoFactor
                ChiFactor=0.5d0*PB%Charge/PB%Mass*PB%dt*PB%dt/Epsilon
                FO%ChiOne=FO%RhoOne*ChiFactor
                Return
				End subroutine WeightP2CParticleBundle
				

         subroutine WeightP2CParticleBundle2(PB,FO)
		       Class(ParticleBundle), intent(inout) :: PB
               Type(FieldOne), intent(inout) :: FO
               Integer(4) :: i
               Real(8) :: RhoFactor,ChiFactor,energyFactor
               FO%RhoOne=0.d0
               FO%ChiOne=0.d0
					FO%EnergyOne=0.d0
					!FO%JxOne=0.d0
					FO%TOne=0.d0
                Do i=1,PB%NPar
                      Call PB%PO(i)%WeightP2C2(FO)
                End Do
                RhoFactor=PB%Charge*PB%Weight
                FO%RhoOne=FO%RhoOne*RhoFactor
                ChiFactor=0.5d0*PB%Charge/PB%Mass*PB%dt*PB%dt/Epsilon
                FO%ChiOne=FO%RhoOne*ChiFactor
					 FO%ChiOne(1)=2.0*FO%ChiOne(1)
					 FO%ChiOne(FO%Nx)=2.0*FO%ChiOne(FO%Nx)
					 
					 energyFactor = 0.5*PB%Mass*PB%weight*PB%VFactor*PB%VFactor
					 FO%EnergyOne= FO%EnergyOne*energyFactor
					 FO%EnergyOne(1)= 2.0*FO%EnergyOne(1)*PB%Dx
					 FO%EnergyOne(FO%Nx)= 2.0*FO%EnergyOne(FO%Nx)*PB%Dx
					 do i=1,FO%Nx
                   if (ABS(FO%RhoOne(i))>0.d0) then
                        FO%TOne(i) = FO%EnergyOne(i)/ABS(FO%RhoOne(i))
                      end if
					end do
					!FO%JxOne=FO%JxOne*PB%Vfactor*PB%Weight*PB%charge
					!FO%JxOne(1)=2.0*FO%JxOne(1)
					!FO%JxOne(FO%Nx)=2.0*FO%JxOne(FO%Nx)
					Return
				End subroutine WeightP2CParticleBundle2
						 
						 
            subroutine DumpParticleBundle(PB,Mode)
                implicit none
                Class(ParticleBundle),intent(inout) :: PB
                Type(ParticleBundle), Allocatable :: TempPB
                Integer(4),intent(in)  :: Mode
                Character(len=99) :: Filename
                Integer(4) :: i,NameIndex
                Real(8) :: XFactor,YFactor
                Integer(4),save :: NameTimer=1
                If (Mode==0) Then
                    NameIndex = DefaultNameIndex
                else
                    NameIndex = DefaultNameIndexInit + ModeMultiplier * Mode + NameTimer
                    NameTimer = NameTimer + 1
                End If 
                
                If (.Not. Allocated(TempPB))  Allocate(TempPB)
                
                TempPB = PB
                Call TempPB%PosRes
                Call TempPB%VelRes
                Write(filename,*) NameIndex,trim(TempPB%SO%Name),".dat"
                Write(*,*) 'Saving ', trim(TempPB%SO%Name), TempPB%NPar,' Please Wait...' 
                  !open (10,file=filename)
                !Write(10,*) TempPB%NPar,TempPB%NParNormal
                !Write(10,*) TempPB%Charge,TempPB%Mass,TempPB%Weight
                !!Dim=Sizeof(TempPB%PO(1))/8_4
                !do i=1,TempPB%NPar
                !     Write(10,FMt="(*(es21.14,1x))")  TempPB%PO(i)
                !End do
                !close(10)
                 
                open (10,FORM='UNFORMATTED',file=filename)
                    Write(10) TempPB%NPar,TempPB%NParNormal
                    Write(10) TempPB%Charge,TempPB%Mass,TempPB%Weight
                    !Dim=Sizeof(TempPB%PO(1))/8_4
                    do i=1,TempPB%NPar
                         Write(10)  TempPB%PO(i)
                    End do
                close(10)
                Write(*,*) 'Saving ', trim(TempPB%SO%Name),' Complete!' 
                
                If ( Allocated(TempPB)) Deallocate(TempPB)
                
                return
            end subroutine  DumpParticleBundle
   
            subroutine LoadParticleBundle(PB,Status)
                implicit none
                Class(ParticleBundle),intent(inout) :: PB
                Logical,intent(inout) ::  Status
                Logical :: alive
                Character(len=99) :: Filename
                Integer(4) :: i,NameIndex,NPArMax
                
                NameIndex=DefaultNameIndex
                Write(filename,*) NameIndex,trim(PB%SO%Name),".dat"
                Inquire(file=filename,exist=alive)
                If(alive) then
                       !Open (10,file=filename)
                       !Read(10,*) PB%NPar,PB%NParNormal
                       !Read(10,*)  PB%Charge, PB%Mass, PB%Weight
                       !Write(*,*) 'Loading ', trim(PB%SO%Name), PB%NPar,' Please Wait...' 
                       ! NPArMax=Ceiling(2.0*PB%NParNormal)
                       ! If(Allocated(PB%PO)) Deallocate(PB%PO)
                       ! Allocate(PB%PO(NPArMax))
                       !do i=1,PB%NPar
                       !      Read(10,*)  PB%PO(i)
                       !end do
                       !close(10)
                       
                       Open (10,FORM='UNFORMATTED',file=filename)
                            Read(10) PB%NPar,PB%NParNormal
                            Read(10)  PB%Charge, PB%Mass, PB%Weight
                            Write(*,*) 'Loading ', trim(PB%SO%Name), PB%NPar,' Please Wait...' 
                            NPArMax=Ceiling(2.0*PB%NParNormal)
                            If(Allocated(PB%PO)) Deallocate(PB%PO)
                            Allocate(PB%PO(NPArMax))
                            do i=1,PB%NPar
                                    Read(10)  PB%PO(i)
                            end do
                       close(10)
                       Status=.False.
                       Write(*,*) 'Loading ', trim(PB%SO%Name),' Complete!'   
                 else
                       Write(*,*) 'Can not find the file for ', trim(PB%SO%Name),' the particle will be randomly initilalized.'    
                       Status = .True. 
                 End if        
                return
            end subroutine  LoadParticleBundle
            
            Subroutine ParticleBundleNormalization(PB,NParNew) 
               Implicit none
               Class(ParticleBundle),intent(inout) :: PB
               Integer(4),intent(in) :: NParNew
               Type(ParticleOne) :: ParticleTemp
               Integer(4) :: i,NTemp,Ndiff,Index
               NTemp = PB%NPar
               PB%Weight = PB%Weight * dble(NTemp)/dble(NParNew)
               If(NParNew>NTemp) then
                   Ndiff = NParNew-NTemp
                   do  i=1,Ndiff
                       CALL RANDOM_NUMBER(R)
                       Index = Ceiling(R*NTemp)
                       ParticleTemp = PB%PO(i)
                       Call PB%AddOne(ParticleTemp)
                   end do
                else
                    Ndiff=NTemp-NParNew
                    do i=1,NDiff
                        CALL RANDOM_NUMBER(R)
                        Index = Ceiling(R*PB%NPar)
                        Call PB%DelOne(Index)
                    end do
                end if  
                return 
            end  Subroutine  ParticleBundleNormalization
            
           ! !> ab.ZWZ for unequal weight particle deleting
            Subroutine ParticleBundleNormalizationUnequalWeight(PB,NParNew,nWeight) 
                Implicit none
                Class(ParticleBundle),intent(inout) :: PB
                Integer(4),intent(in) :: NParNew
                Integer(4),intent(in) :: nWeight
                Type(ParticleOne) :: ParticleTemp
                Integer(4) :: i,j,k,NTemp,Ndiff,Index
                Integer(4) :: WeightCount(nWeight)
                Integer(4) :: iPart,iPartTemp
                Real(8) :: Res
               
                
                !Call ParticleWeightQuickSort(PB,1,PB%Npar)
                Call ParticleWeightShellSort(PB)
                
                WeightCount(:) = 0
                j = 1
                WeightCount(1) = 1
                Do i = 2,PB%Npar
                    If (Abs(PB%PO(i)%WQ - PB%PO(i-1)%WQ)<1E-5) Then
                        WeightCount(j) = WeightCount(j) + 1
                    Else
                        If (j < nWeight) Then
                            j = j + 1
                            WeightCount(j) = 1
                        Endif
                    End If
                End Do
                
                !iPart = 0
                !iPartTemp = 0
                !NTemp = PB%NPar
                !Do i = 1,nWeight
                !    If (WeightCount(i) > 50 ) Then  !> only if amount of particles with certain weight over 50
                !        Do j = 1,WeightCount(i)
                !            iPart = iPart + 1
                !            PB%PO(iPart)%W = PB%PO(iPart)%W * dble(NTemp)/dble(NParNew)
                !        End Do
                !        Res = WeightCount(i) * Abs(Real(NParNew-NTemp)/Real(NTemp))
                !        Ndiff = Res
                !        Res = Res - Ndiff
                !        CALL RANDOM_NUMBER(R)
                !        If (R < Res) Ndiff = Ndiff + 1
                !        
                !        If (NparNew > NTemp) Then
                !            Do  k=1,Ndiff
                !                CALL RANDOM_NUMBER(R)
                !                Index = iPartTemp + Ceiling(R*WeightCount(i))
                !                ParticleTemp = PB%PO(Index)
                !                Call PB%AddOne(ParticleTemp)
                !            End Do
                !        Else
                !            do k=1,NDiff
                !                CALL RANDOM_NUMBER(R)
                !                Index = iPartTemp + Ceiling(R*WeightCount(i))
                !                PB%PO(index)%W = -2000
                !            end do
                !        End If
                !    End If
                !    iPartTemp = iPartTemp + WeightCount(i)
                !End Do  
                
                iPart = PB%Npar + 1
                iPartTemp = PB%Npar 
                NTemp = PB%NPar
                Do i = nWeight,1,-1
                    iPartTemp = iPartTemp - WeightCount(i)
                    If (WeightCount(i) > 50 ) Then  !> only if amount of particles with certain weight over 50
                        Do j = 1,WeightCount(i)
                            iPart = iPart - 1 
                            PB%PO(iPart)%WQ = PB%PO(iPart)%WQ * dble(NTemp)/dble(NParNew)
                        End Do
                        Res = WeightCount(i) * Abs(Real(NParNew-NTemp)/Real(NTemp))
                        Ndiff = Res
                        Res = Res - Ndiff
                        CALL RANDOM_NUMBER(R)
                        If (R < Res) Ndiff = Ndiff + 1
                        
                        If (NparNew > NTemp) Then
                            Do  k=1,Ndiff
                                CALL RANDOM_NUMBER(R)
                                Index = iPartTemp + Ceiling(R*WeightCount(i))
                                ParticleTemp = PB%PO(Index)
                                Call PB%AddOne(ParticleTemp)
                            End Do
                        Else
                            do k=1,NDiff
                                CALL RANDOM_NUMBER(R)
                                Index = iPartTemp + Ceiling(R*WeightCount(i))
                                Call PB%DelOne(index)
                                Call SwapParticleOne(PB%PO(index),PB%PO(iPartTemp+WeightCount(i)))
                                WeightCount(i) = WeightCount(i) - 1
                            end do
                        End If
                    Else
                        iPart = iPartTemp + 1
                    End If
                    
                End Do  
                !    
                !Do i=PB%NPar,1,-1
                !    If (PB%PO(i)%W<0) then
                !        Call PB%DelOne(i)  
                !    End if
                !End do
                return 
           end  Subroutine  ParticleBundleNormalizationUnequalWeight
            
         !    Subroutine InitializationParticleBundleIndex(PBI)
		       !Class(ParticleBundleIndex), intent(inout) :: PBI
         !       If(Allocated(PBI%POI)) Deallocate(PBI%POI)
         !       Allocate(PBI%POI(PBI%NParNormal))
         !       PBI%NPar=0
         !       return
         !    End Subroutine InitializationParticleBundleIndex
             
             Subroutine InitializationParticleBundleIndex(PBI,PB,CollisionRatio)
		        Class(ParticleBundleIndex), intent(inout) :: PBI
		        Type(ParticleBundle), intent(in) :: PB
		        Real(8),intent(in) :: CollisionRatio
                PBI%NParNormal = Int(2*PB%Npar*CollisionRatio)+1
                If(Allocated(PBI%POI)) Deallocate(PBI%POI)
                Allocate(PBI%POI(PBI%NParNormal))
                PBI%NPar=0
                return
            End Subroutine InitializationParticleBundleIndex

            Subroutine  AddParticleOneParticleBundleIndex(PBI,POI)
             Implicit none
                 Class(ParticleBundleIndex),intent(inout) ::  PBI
                 Type(ParticleOneIndex),intent(in) ::  POI
                  PBI%NPar = PBI%NPar+1
                  Call PBI%POI(PBI%NPar)%Copy(POI)
             return   
             End  Subroutine  AddParticleOneParticleBundleIndex
     
             Subroutine  DelParticleOneParticleBundleIndex(PBI,NDel)
             Implicit none
                 Class(ParticleBundleIndex),intent(inout) ::  PBI
                 Integer(4),intent(in) ::  NDel
                 Call PBI%POI(NDel)%Copy(PBI%POI(PBI%NPar))
                 PBI%NPar=PBI%NPar-1
             return   
             End  Subroutine  DelParticleOneParticleBundleIndex
            
            
           ! 快速排序法的子程序
            RECURSIVE SUBROUTINE ParticleWeightQuickSort(PB,S,E)
                IMPLICIT NONE
    
                Class(ParticleBundle),intent(inout) :: PB
                Integer(4),intent(in) :: S ! 传入的参数, 这一组的类型起始位置
                Integer(4),intent(in) :: E ! 传入的参数, 这一组的类型结束位置
                Integer(4) :: L,R ! 用来找A(L)>K及A(R)<K时用的
                Real(8) :: K ! 记录键值A(S)
    
                ! 首先要先给定L,R的初值. L要从头开始,E则要从尾开始
                L=S
                R=E
                ! RIGHT值 > LEFT值 时才有必要进行排序 
                IF ( R<=L ) RETURN
    
                K=PB%PO(S)%WQ ! 设定键值
                DO WHILE(.TRUE.)
                    ! 找出PB%PO(L)%W>K的所在
                    DO WHILE( .TRUE. )
                        L=L+1
                        IF ( (PB%PO(L)%WQ > K) .OR. (L>=E) ) EXIT
                    END DO
                    ! 找出PB%PO(R)%W<K的所在
                    DO WHILE( .TRUE. )
                        !R=R-1
                        IF ( (PB%PO(R)%WQ < K) .OR. (R<=S) ) EXIT
                        R=R-1
                    END DO
        
                    ! 如果RIGHT 跑到 LEFT的左边时, 循环就该结束了
                    IF ( R <= L ) EXIT
                    ! exchange PO(L) and PO(R)
                    Call SwapParticleOne(PB%PO(L),PB%PO(R))
                END DO
    
                !print*,'before swap',' S=',S,' R=',R
                ! exchange PO(S) and PO(R)
                Call SwapParticleOne(PB%PO(S),PB%PO(R))
                ! 把R之前的数据重新分组,再做排序
                CALL ParticleWeightQuickSort(PB,S,R-1)
                ! 把R之后的数据重新分组,再做排序
                CALL ParticleWeightQuickSort(PB,R+1,E)
    
                RETURN
            END SUBROUTINE ParticleWeightQuickSort
             
           Subroutine ParticleWeightShellSort(PB)
                implicit none
                Class(ParticleBundle),intent(inout) :: PB
                Integer(4) :: h,i,j
                    
                !> 1. confirm the increase amount according to length
                h = 1
                Do While (h < Real(PB%Npar)/2)
                    h = 2 * h + 1
                End Do
                
                !> 2. shell sort
                Do While ( h >= 1) 
                    !> sort
                    !> 2.1. find the element waiting for being inserted
                    Do i = h+1, PB%Npar
                        !> 2.2. insert the element into array
                        j = i
                        Do While( j > h)
                            !> for inserting PB%PO(j), compare PB%PO(j) and PB%PO(j-h)
                            If ( PB%PO(j-h)%WQ > PB%PO(j)%WQ ) Then
                                !> exchange elements
                                Call SwapParticleOne(PB%PO(j-h),PB%PO(j))
                            Else
                                !> end the loop for the element having been inserted to appropriate position
                                Exit
                            End if
                            j = j - h
                        End Do
                    End Do
                    !> reduce h
                    h = h / 2
                End Do
                
               
           End Subroutine ParticleWeightShellSort

              
            !> ab.ZWZ 
            Subroutine Coalescing(PB,Nx,Ny)
                Use Domain_2D,Only: Vert_o, hxi
                Implicit none
                Real(8),parameter::Ratio=0.1d0
                Integer(4),intent(in)::Nx,Ny
                Integer(4)::i,j,k,N1,N2
                Type(ParticleBundle),intent(inout)::PB
                Type(ParticleBundle) :: ParticleTempBundle !这个用来存不需要合并的粒子和合并后的粒子。
                Type(ParticleOne),Allocatable ::  ParticleTemp(:,:),TempParticle1(:,:),TempParticle2(:,:),TempParticle3(:,:)
                Type(ParticleOne) ::  TempParticleA,TempParticleB,TempParticleC,TempParticleD,TempParticleE !ABC用来存合并前的粒子，DE用来存合并后的粒子
                Integer(4),Allocatable :: Counter(:,:)
                Real(8) :: WQNormal,WQ,MxFactor,ExFactor,MyFactor,EyFactor,MzFactor,EzFactor
                Real(8) :: Radius,PositionFactorZ,AccFactorX,AccFactorY,AccFactorZ
                Real(8) :: Rnum,AFai,CosFai,SinFai
                ParticleTempBundle=PB !因为下面的循环里面会有AddParticle操作，而循环要求PB%Npar不能改变
                ParticleTempBundle%NPar=0  
                Allocate(ParticleTemp(Nx,Ny))
                Allocate(TempParticle1(Nx,Ny))
                Allocate(TempParticle2(Nx,Ny))
                Allocate(TempParticle3(Nx,Ny))
                Allocate(Counter(Nx,Ny))
                
                Counter(1:Nx,1:Ny)=0

                Do i=1,PB%Npar
                    !N1=Ceiling(PB%PO(i)%Z)
                    !N2=Ceiling(PB%PO(i)%Radius) 
                    !N1=Ceiling(PB%PO(i)%X)
                    !N2=Ceiling(PB%PO(i)%Y) 
                    N1 = Int((PB%PO(i)%X - Vert_o(1))*hxi(1))
                    N2 = Int((PB%PO(i)%Y - Vert_o(2))*hxi(2))
    
                    !if(PB%PO(i)%Y<=1.d0) then
                    !    WQNormal=1.d0/4.d0*WQOne !*Density/Dble(PPG)*2.d0*PI*dx*dx*dx !这个地方是如何推出来的，还要再思考一下
                    !else
                    !    WQNormal=WQone*dble(N2) !*Density/Dble(PPG)*2.d0*PI*dx*dx*dx
                    !end if
                    WQNormal = PB%WQOne * (2*Dble(N2)-1)

                    if(PB%PO(i)%WQ<Ratio*WQNormal)then
                        Counter(N1,N2)=Counter(N1,N2)+1
                        if (counter(N1,N2)==1)  then
                            TempParticle1(N1,N2)=PB%PO(i)
                        else if (counter(N1,N2)==2)  then
                            TempParticle2(N1,N2)=PB%PO(i)
                        else if (counter(N1,N2)==3)  then
                            TempParticle3(N1,N2)=PB%PO(i)
                            TempParticleA=TempParticle1(N1,N2)
                            TempParticleB=TempParticle2(N1,N2)
                            TempParticleC=TempParticle3(N1,N2)
                            WQ=TempParticleA%WQ+TempParticleB%WQ+TempParticleC%WQ
                                             
                            TempParticleD%WQ=0.5d0*WQ
                            TempParticleE%WQ=0.5d0*WQ

                            !PositionFactor和AccFactor是用来计算坐标和加速度的因子，坐标和加速度的的算法是一样的
                            !PositionFactorX=TempParticleA%WQ*TempParticleA%X+TempParticleB%WQ*TempParticleB%X+TempParticleC%WQ*TempParticleC%X
                            !PositionFactorY=TempParticleA%WQ*TempParticleA%Y+TempParticleB%WQ*TempParticleB%Y+TempParticleC%WQ*TempParticleC%Y
                            !Radius=TempParticleA%WQ*TempParticleA%Radius+TempParticleB%WQ*TempParticleB%Radius+TempParticleC%WQ*TempParticleC%Radius
                            !PositionFactorZ=TempParticleA%WQ*TempParticleA%Z+TempParticleB%WQ*TempParticleB%Z+TempParticleC%WQ*TempParticleC%Z
                            Radius=TempParticleA%WQ*TempParticleA%Y+TempParticleB%WQ*TempParticleB%Y+TempParticleC%WQ*TempParticleC%Y
                            PositionFactorZ=TempParticleA%WQ*TempParticleA%X+TempParticleB%WQ*TempParticleB%X+TempParticleC%WQ*TempParticleC%X
                            AccFactorX=TempParticleA%WQ*TempParticleA%Ax+TempParticleB%WQ*TempParticleB%Ax+TempParticleC%WQ*TempParticleC%Ax
                            AccFactorY=TempParticleA%WQ*TempParticleA%Ay+TempParticleB%WQ*TempParticleB%Ay+TempParticleC%WQ*TempParticleC%Ay
                            AccFactorZ=TempParticleA%WQ*TempParticleA%Az+TempParticleB%WQ*TempParticleB%Az+TempParticleC%WQ*TempParticleC%Az
                            !if (WQ<0.d0.or.WQ==0.d0) then 
                            !    write(*,*) 'WQerrrrrr',i
                            !end if
                            !TempParticleD%X=PositionFactorX/WQ !+0.1d0*PositionFactorX  !合并后粒子的坐标老师让考虑一下类似二维情况下网格内粒子的电荷分布到四个格点的情况，要求网格内的粒子的坐标满足电荷分配到各个格点守恒的关系。
                            !TempParticleD%Y=PositionFactorY/WQ ! + 0.1d0*PositionFactorY
                            !TempParticleD%Radius=Radius/WQ
                            !TempParticleD%Z=PositionFactorZ/WQ ! + 0.1d0*PositionFactorZ
                            TempParticleD%Y=Radius/WQ
                            TempParticleD%X=PositionFactorZ/WQ ! + 0.1d0*PositionFactorZ
                            Call DRandom(Rnum)
                            AFai=2.d0*Pi*Rnum
                            CosFai=DCos(AFai)
                            SinFai=DSin(AFai)  
                            !TempParticleD%X=TempParticleD%Radius*CosFai
                            !TempParticleD%Y=TempParticleD%Radius*SinFai
                            TempParticleD%Z = AFai
                                             
                            !TempParticleE%X=PositionFactorX/WQ ! -0.1d0*PositionFactorX 
                            !TempParticleE%Y=PositionFactorY/WQ ! - 0.1d0*PositionFactorY
                            !TempParticleE%Radius=Radius/WQ
                            !TempParticleE%Z=PositionFactorZ/WQ ! - 0.1d0*PositionFactorZ
                            TempParticleE%Y=Radius/WQ
                            TempParticleE%X=PositionFactorZ/WQ ! - 0.1d0*PositionFactorZ
                            Call DRandom(Rnum)
                            AFai=2.d0*Pi*Rnum
                            CosFai=DCos(AFai)
                            SinFai=DSin(AFai)  
                            !TempParticleE%X=TempParticleE%Radius*CosFai
                            !TempParticleE%Y=TempParticleE%Radius*SinFai
                            TempParticleE%Z = AFai

                            TempParticleD%Ax=AccFactorX/WQ !+0.1d0*AccFactorX
                            TempParticleD%Ay=AccFactorY/WQ !+0.1d0*AccFactorY
                            TempParticleD%Az=AccFactorZ/WQ !+0.1d0*AccFactorZ
                            TempParticleE%Ax=AccFactorX/WQ !-0.1d0*AccFactorX
                            TempParticleE%Ay=AccFactorY/WQ !-0.1d0*AccFactorY
                            TempParticleE%Az=AccFactorZ/WQ !-0.1d0*AccFactorZ
                             
                            ExFactor=TempParticleA%WQ*TempParticleA%Vx*TempParticleA%Vx+TempParticleB%WQ*TempParticleB%Vx*TempParticleB%Vx+TempParticleC%WQ*TempParticleC%Vx*TempParticleC%Vx
                            MxFactor=TempParticleA%WQ*TempParticleA%Vx+TempParticleB%WQ*TempParticleB%Vx+TempParticleC%WQ*TempParticleC%Vx
                            EyFactor=TempParticleA%WQ*TempParticleA%Vy*TempParticleA%Vy+TempParticleB%WQ*TempParticleB%Vy*TempParticleB%Vy+TempParticleC%WQ*TempParticleC%Vy*TempParticleC%Vy
                            MyFactor=TempParticleA%WQ*TempParticleA%Vy+TempParticleB%WQ*TempParticleB%Vy+TempParticleC%WQ*TempParticleC%Vy
                            EzFactor=TempParticleA%WQ*TempParticleA%Vz*TempParticleA%Vz+TempParticleB%WQ*TempParticleB%Vz*TempParticleB%Vz+TempParticleC%WQ*TempParticleC%Vz*TempParticleC%Vz
                            MzFactor=TempParticleA%WQ*TempParticleA%Vz+TempParticleB%WQ*TempParticleB%Vz+TempParticleC%WQ*TempParticleC%Vz

                            TempParticleD%Vx=(MxFactor-Dsqrt(abs(WQ*ExFactor-MxFactor*MxFactor)))/WQ
                            TempParticleD%Vy=(MyFactor-Dsqrt(abs(WQ*EyFactor-MyFactor*MyFactor)))/WQ
                            TempParticleD%Vz=(MzFactor-Dsqrt(abs(WQ*EzFactor-MzFactor*MzFactor)))/WQ
                            TempParticleE%Vx=(MxFactor+Dsqrt(abs(WQ*ExFactor-MxFactor*MxFactor)))/WQ
                            TempParticleE%Vy=(MyFactor+Dsqrt(abs(WQ*EyFactor-MyFactor*MyFactor)))/WQ
                            TempParticleE%Vz=(MzFactor+Dsqrt(abs(WQ*EzFactor-MzFactor*MzFactor)))/WQ
                            !write(*,*) WQ*EzFactor,MzFactor*MzFactor,WQ*EzFactor-MzFactor*MzFactor
                                             
                            !if (IsNan(TempParticleD%Radius)==.true..or.IsNan(TempParticleE%Vx)==.true.) then
                            !       write(*,*) 'Before errrrrr',i
                            !   endif
                            !if (WQ*ExFactor-MxFactor*MxFactor<0.d0.or.WQ*EyFactor-MyFactor*MyFactor<0.d0) then
                            !    write(*,*) 'Beforessss errrrrr',i,WQ*ExFactor,MxFactor*MxFactor
                            !endif
                                             
                                             
                            !TempParticleD%Radius=dsqrt(TempParticleD%X*TempParticleD%X+TempParticleD%Y*TempParticleD%Y+TempParticleD%Z*TempParticleD%Z)
                            !TempParticleE%Radius=dsqrt(TempParticleE%X*TempParticleE%X+TempParticleE%Y*TempParticleE%Y+TempParticleE%Z*TempParticleE%Z)
                            !Call AddParticle(TempParticleD,ParticleTempBundle)
                            !Call AddParticle(TempParticleE,ParticleTempBundle)
                            Call ParticleTempBundle%AddOne(TempParticleD)
                            Call ParticleTempBundle%AddOne(TempParticleE)
                            Counter(N1,N2)=0
                            !Write (*,*) 'Ex',ExFactor,0.5d0*WQ*(TempParticleD%Vx*TempParticleD%Vx+TempParticleE%Vx*TempParticleE%Vx)
                            ! Write (*,*) 'Mx',MxFactor,0.5d0*WQ*(TempParticleD%Vx+TempParticleE%Vx)
                            !Write (*,*) 'Qx',WQ,(TempParticleD%WQ+TempParticleE%Wq)
                            !pause

                        end if
                    else   
                        !Call AddParticle(PB%PO(i),ParticleTempBundle)   
                        Call ParticleTempBundle%AddOne(PB%PO(i))
                    end if
  
                end do

                Do j=1,Ny
                    Do k=1,Nx    
                        If (counter(k,j)==1) then
                            !If (TempParticle1(k,j)%Z>0.d0) then
                                !Call AddParticle(TempParticle1(k,j),ParticleTempBundle)
                                Call ParticleTempBundle%AddOne(TempParticle1(k,j))
                            !end if
                        Else if (counter(k,j)==2) then
                            !If (TempParticle1(k,j)%Z>0.d0) then 
                                !call AddParticle(TempParticle1(k,j),ParticleTempBundle)
                                Call ParticleTempBundle%AddOne(TempParticle1(k,j))
                            !end if
                            !If (TempParticle2(k,j)%Z>0.d0) then 
                                !call AddParticle(TempParticle2(k,j),ParticleTempBundle)
                                Call ParticleTempBundle%AddOne(TempParticle2(k,j))
                            !end if
                               
                        End if 
                    End do
                End do 
                PB=ParticleTempBundle
                Do i=1,PB%Npar
                   !if (IsNan(TempParticleD%Radius)==.true..or.IsNan(TempParticleE%Vx)==.true.) then
                   If (IsNan(TempParticleD%Y)==.true..or.IsNan(TempParticleE%Vx)==.true.) then
                       write(*,*) 'After errrrrr',i
                   Endif
                Enddo
            
                Deallocate(ParticleTemp,TempParticle1,TempParticle2,TempParticle3,Counter)
            end subroutine Coalescing
    
            Subroutine WeightRenormalization(PB)
                Use Domain_2D,Only: Vert_o, hxi
                !Use ParticleModule 
                implicit none
                Type(ParticleBundle),intent(inout) :: PB
                Type(ParticleOne) :: TempParticle
                Real(8) :: WQ
                real(8) :: RTemp
                Integer(4) :: i,j,N,Ncr  
                Do i=1,PB%Npar
                    !RTemp=dsqrt(PB%PO(i)%X*PB%PO(i)%X+PB%PO(i)%Y*PB%PO(i)%Y) 
                    !N=Int(RTemp)
                    N = Int((PB%PO(i)%Y - Vert_o(2))*hxi(2))
                    !If(N>=1) then
                    !        WQ=WqOne*Dble(N)
                    !Else
                    !        Wq=1.d0/6.d0*WqOne
                    !End if    
                    Wq = (2*Dble(N)-1)*PB%WQOne
                    Ncr=Int(PB%PO(i)%WQ/WQ)
                    If(Ncr>1) then
                        do j=1,Ncr-1
                            PB%PO(i)%WQ=PB%PO(i)%WQ-WQ
                            TempParticle=PB%PO(i)
                            TempParticle%WQ=WQ
                            !Call AddParticle(TempParticle,PB)
                            Call PB%AddOne(TempParticle)
                        End do
                    End if
                End do
                return
            End  subroutine WeightRenormalization
            
END Module ModuleParticleBundle
