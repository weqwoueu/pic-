Module DiagnosticsMomentum
    Use ModuleGrid
    Use ModuleField
    Use ModuleParticleBundle
    Use ModuleParticleBoundary
    Implicit none
    Type ParticleMomentumOne!(Nx)!!!!
                        !Integer(4),Len :: Nx=100
                        !Integer(4) :: IOIndex=1
                        !Integer(4) :: XStart=0,XEnd=NxMax-1
                        Integer(4) :: Nx=NxMax,Timer=0
                        Real(8) :: Dx=Inputdx,Dt=Inputdt
                        !Real(8) :: QdM
                        Real(8) ::  RhoOne(1:NxMax),ChiOne(1:NxMax)
                        Real(8) ::  JxOne(1:NxMax),JyOne(1:NxMax),JzOne(1:NxMax)
                        Real(8) ::  TOne(1:NxMax),EnergyOne(1:NxMax)
                   EndType ParticleMomentumOne
    contains

     Subroutine  DiagParticleFieldPeriod(GD,NSpecy,PB,FG,Mode)
         Implicit none
         Class(*),intent(inout)  :: GD
         Integer(4),intent(in) ::  NSpecy
         Type(ParticleBundle),intent(inout) :: PB(0:NSpecy)
         Type(Field),intent(in) :: FG
         Integer(4),intent(in) ::  Mode
         !Integer(4),parameter :: NSpecyMax=2_4 
         Type(ParticleMomentumOne),save :: TempPMO(0:NspacyMax)
         Integer(4) :: i,j,Shift
         Real(8) :: HeatingRate(FG%Nx),EfieldPower(FG%Nx),FieldheatingRate(FG%Nx),CollisionPowerlossRate(FG%Nx),DeltaEpGain(FG%Nx),DeltaEpLoss(FG%Nx),DesplaceCurrent(FG%Nx)
         Select Type (GD)
             Type is (Grid1D(*,*))
                 Select Case (Mode)
                    case(0)
                        Shift=1
                        Do i=0, NSpecy
                            Call WeightingParticleMomentum(PB(i),TempPMO(i),FG)
                               Call GD%Update(GD%Nx,TempPMO(i)%RhoOne,Shift) !带电粒子数量密度
										 Call GD%Update(GD%Nx,TempPMO(i)%JxOne,Shift)   !带电粒子电荷流
										 Call GD%Update(GD%Nx,TempPMO(i)%EnergyOne/JtoeV,Shift)
                               !Call GD%Update(GD%Nx,TempPMO(i)%TOne,Shift)
								   	  HeatingRate=TempPMO(i)%JxOne*FG%Ex
								   	 Call GD%Update(GD%Nx,HeatingRate,Shift)
								 End do
									 DesplaceCurrent = Epsilon*(FG%Ex-FG%Ex2)/FG%Dt
									 Call GD%Update(GD%Nx,DesplaceCurrent,Shift)
									 EfieldPower = DesplaceCurrent*(FG%Ex+FG%Ex2)/2
									 Call GD%Update(GD%Nx,EfieldPower,Shift)
                            Call GD%Update(GD%Nx,FG%Ex,Shift)
                            Call GD%Update(GD%Nx,FG%Phi,Shift)
                            Call GD%Update(GD%Nx,FG%Rho,Shift)
                        GD%Timer=GD%Timer+1
                     Case(1) 
                             Call GD%Rescale
                             Call GD%Dump(1)
                             GD%Value=0.d0
                             GD%Timer=0
                     Case(2)
                         Call GD%Dump(0)
                     case default
                         
                     End Select
                Type is (Grid2D(*,*,*))
                    Select Case (Mode)
                        case(0)
                            Shift=1
                             Do i=0, NSpecy
										 !Call WeightingParticleMomentum(PB(i),TempPMO)
                               Call GD%Update(GD%Nx,TempPMO(i)%RhoOne,Shift)
										 Call GD%Update(GD%Nx,TempPMO(i)%JxOne,Shift)
										 Call GD%Update(GD%Nx,TempPMO(i)%EnergyOne/JtoeV,Shift)
                               !Call GD%Update(GD%Nx,TempPMO(i)%TOne,Shift)
								   	  HeatingRate=TempPMO(i)%JxOne*FG%Ex
								   	 Call GD%Update(GD%Nx,HeatingRate,Shift)
									 End do
									 DesplaceCurrent = Epsilon*(FG%Ex-FG%Ex2)/FG%Dt
									 Call GD%Update(GD%Nx,DesplaceCurrent,Shift)
									 EfieldPower = DesplaceCurrent*(FG%Ex+FG%Ex2)/2
									 Call GD%Update(GD%Nx,EfieldPower,Shift)
                            Call GD%Update(GD%Nx,FG%Ex,Shift)
                            Call GD%Update(GD%Nx,FG%Phi,Shift)
                            Call GD%Update(GD%Nx,FG%Rho,Shift)
                            GD%Timer=GD%Timer+1
                        Case(1) 
                                 Call GD%Rescale
                                 Call GD%Dump(1)
                                 GD%Value=0.d0
                                 GD%Timer=0
                         Case(2)
                             Call GD%Dump(0) 
                         case default
                     End Select     
                End select
        Return    
     End Subroutine  DiagParticleFieldPeriod

  subroutine WeightingParticleMomentum(PB,PMO,FG)
                implicit none
                Type(ParticleBundle),intent(inout) :: PB
                Type(ParticleMomentumOne),intent(inout) :: PMO
					 Type(Field),intent(in) :: FG
					 Type(ParticleOne) :: TempPO
                Real(8) :: RhoFactor,ChiFactor,JFactor,EnergyFactor
					 Real(8) :: Vdirft(FG%Nx)
					 Real(8) ::Ex,Ey,Ez,Bx,By,Bz,EFactor,BFactor
                Real(8) :: S1,S2,Energy
                Integer(4) :: i,N
                PMO%RhoOne=0.d0
                PMO%JxOne=0.d0
					 PMO%EnergyOne  = 0.d0
                PMO%TOne=0.d0
					 Vdirft = 0.d0
					 EFactor=PB%Charge/PB%Mass*PB%dt/(PB%dx/PB%dt)
					 BFactor=PB%Charge/PB%Mass*PB%dt
                do i=1,PB%Npar
                   N=Ceiling(PB%PO(i)%X)
                   S1=Dble(N)-PB%PO(i)%X
                   S2=1.d0-S1
                   PMO%RhoOne(N)=PMO%RhoOne(N)+S1
                   PMO%RhoOne(N+1)=PMO%RhoOne(N+1)+S2
					!!!!!!!Volocity diagnosis needed update
						 Ex=FG%Ex(N)*S1+FG%Ex(N+1)*S2 	!CALL PB%PO(i)WeightC2PES(FG,Ex)
						 Ex=Ex*EFactor
						 TempPO = PB%PO(i)
						 if(FG%Electromagnetic) then
                      Ey=FG%Ey(N)*S1+FG%Ey(N+1)*S2
                      Ez=FG%Ez(N)*S1+FG%Ez(N+1)*S2
                      Bx=FG%Bx(N)*S1+FG%Bx(N+1)*S2	!call PB%PO(i)WeightC2PEM(FG,Ex,Ey,Ez,Bx,By,Bz) 
                      By=FG%By(N)*S1+FG%By(N+1)*S2
                      Bz=FG%Bz(N)*S1+FG%Bz(N+1)*S2
                      Ey=Ey*EFactor
                      Ez=Ez*EFactor
                      Bx=Bx*BFactor
                      By=By*BFactor
                      Bz=Bz*BFactor
							 call TempPO%MoveEM(Ex,Ey,Ez,Bx,By,Bz)
						 else
							 call TempPO%MoveES(Ex)
						 end if
						 !TempVx = 0.5*(PB%PO(i)%Vx + TempVx)   !control the corrected volocity of particle
                   PMO%JxOne(N)=PMO%JxOne(N)+S1*TempPO%Vx 
                   PMO%JxOne(N+1)=PMO%JxOne(N+1)+S2*TempPO%Vx					
                   Energy= PB%PO(i)%Vx*PB%PO(i)%Vx+PB%PO(i)%Vy*PB%PO(i)%Vy+PB%PO(i)%Vz*PB%PO(i)%Vz
                   PMO%EnergyOne(N)=PMO%EnergyOne(N)+S1*Energy
                   PMO%EnergyOne(N+1)=PMO%EnergyOne(N+1)+S2*Energy
                end do
                RhoFactor=PB%Weight
                PMO%RhoOne=PMO%RhoOne*RhoFactor
					 PMO%RhoOne(1)=2.d0*PMO%RhoOne(1)
                PMO%RhoOne(PMO%Nx)=2.d0*PMO%RhoOne(PMO%Nx)
					
                ChiFactor=0.5d0*PB%Charge/PB%Mass*PB%dt*PB%dt/Epsilon*PB%Charge
                PMO%ChiOne=PMO%RhoOne*ChiFactor
					 PMO%ChiOne(1)=2.d0*PMO%ChiOne(1)
                PMO%ChiOne(PMO%Nx)=2.d0*PMO%ChiOne(PMO%Nx)
					 
                JFactor=PB%Charge*PB%Weight*PB%VFactor
                PMO%JxOne=PMO%JxOne*JFactor
					 PMO%JxOne(1)=2.d0*PMO%JxOne(1)
                PMO%JxOne(PMO%Nx)=2.d0*PMO%JxOne(PMO%Nx)
					
					 
                EnergyFactor=0.5*PB%mass*PB%Vfactor*PB%Vfactor*PB%Weight
					 PMO%EnergyOne=PMO%EnergyOne*EnergyFactor
					 PMO%EnergyOne(1) = 2.0*PMO%EnergyOne(1)
					 PMO%EnergyOne(PMO%Nx) = 2.0*PMO%EnergyOne(PMO%Nx)
					 
                PMO%TOne=PMO%TOne*EnergyFactor/JtoeV
					 PMO%TOne(1)=2.d0*PMO%TOne(1)
                PMO%TOne(PMO%Nx)=2.d0*PMO%TOne(PMO%Nx)
                do i=1,PMO%Nx
                    if (PMO%RhoOne(i)>0.d0) then
                        PMO%TOne(i) = PMO%TOne(i)/PMO%RhoOne(i)
								PMO%EnergyOne(i) = PMO%EnergyOne(i)/PMO%RhoOne(i)
                    end if
                end do
                return
  end subroutine WeightingParticleMomentum
    End Module DiagnosticsMomentum
!!   
!!    
!