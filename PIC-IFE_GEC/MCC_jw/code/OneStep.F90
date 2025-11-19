Module ModuleOneStep
      Use ModuleParticleBundle
      use Particle_Injection
      Use ModuleSpecyOne
      Use ModuleMCCPublic
      Use ModuleOneStepField
      Use ModuleMCCInitialization
		Use ModuleSEE
		 Use ModuleRecombination
        Use ModuleReactionOnePegasus
      
      !Use MoveModule
!      Use Diagnostics
      Implicit none
      ! This section defines the particles.
                  Type(ControlFlow),save :: ControlFlowGlobal
      
                  Type(ParticleBundle),save,Allocatable ::  ParticleGlobal(:)
                  Type(ParticleBoundaryOne),save,Allocatable  :: ParticleBDOneGlobal(:)
                  
                  Type(ElectronIonRecombination),save :: EIRGlobal
                  Type(IonIonRecombination),save :: IIRGlobal
                  !Type(ParticleBoundary),save  :: ParticleBDGlobal
                
    contains
    Subroutine AllInitilalization()
             Implicit none
             Integer(4) :: i
             Call InitializationControlFlow(ControlFlowGlobal)
             Call GasInitPegasus(ControlFlowGlobal)
             Call InitializationField(ControlFlowGlobal)
             Allocate(ParticleGlobal(0:ControlFlowGlobal%Ns))
             Allocate(ParticleBDOneGlobal(0:ControlFlowGlobal%Ns))
             DO i=0,ControlFlowGlobal%Ns
                    Call ParticleGlobal(i)%AllInit(SpecyGlobal(i),ControlFlowGlobal)
                    Call ParticleBDOneGlobal(i)%AllInit(ParticleGlobal(i),ControlFlowGlobal)
				 End do
             Call MCCBundleInit(ControlFlowGlobal,SpecyGlobal,GasGlobal)
             Call EIRGlobal%Init(ControlFlowGlobal,GasGlobal)
             Call EIRGlobal%InitConstant(ParticleGlobal(:))
             Call IIRGlobal%Init(ControlFlowGlobal,GasGlobal)
             Call IIRGlobal%InitConstant(ParticleGlobal(:))
            Return  
    End Subroutine AllInitilalization

    Subroutine OneStep()
            Implicit none
               Integer(4) :: i,j
               do i=0,ControlFlowGlobal%Ns
                     Call ParticleGlobal(i)%MoveES(FieldGlobal)
							!Call ParticleGlobal(i)%MoveEM(FieldGlobal)
                     Call ParticleAborption(ParticleGlobal(i),ParticleBDOneGlobal(i))    !Update particle boundary
                     !if (i==0) then
							!	call ElectronIndcedSEE(ParticleGlobal(i),ParticleBDOneGlobal(i),0)
							!  call ElectronIndcedSEE(ParticleGlobal(i),ParticleBDOneGlobal(i),1)
							!else 
							!	call IonIndcedSEE(ParticleGlobal(i),ParticleBDOneGlobal(i),0)
							!  call IonIndcedSEE(ParticleGlobal(i),ParticleBDOneGlobal(i),1)
							!end if
                     !if(i==Inject_Particle_type) then
                     !    call Particle_Beam_InjectionS(ControlFlowGlobal,ParticleGlobal(i))
                     !end if
                     Call ParticleGlobal(i)%WeightP2C(FieldOneGlobal(i))   !caculate No. of particle to cell !density and current update 
					end do
				 	FieldGlobal%Ex2=FieldGlobal%Ex
                    Call FieldOneStep(ControlFlowGlobal%Ns,FieldOneGlobal,FieldGlobal,FieldBoundaryGlobal,FieldSolverGlobal,ParticleGlobal,ParticleBDOneGlobal)
				    Call MCC(ControlFlowGlobal%Ns,ControlFlowGlobal%Ng,ParticleGlobal,SpecyGlobal,GasGlobal,MCCBundleGlobal) 
					
                    If (ControlFlowGlobal%withRecombination) Then
                    Call EIRGlobal%Onestep(ParticleGlobal,FieldOneGlobal)
                    Call IIRGlobal%Onestep(ParticleGlobal,FieldOneGlobal)
                    End If
                    
      !              if (mod(ControlFlowGlobal%Timer,50)==0) then
						! do i=0,ControlFlowGlobal%Ns
						!     Call ParticleGlobal(i)%WeightP2C2(FieldOneGlobal(i))
						! END DO
						! Call IonIonRecombinationRate(ControlFlowGlobal,ParticleGlobal,FieldOneGlobal,2,1,50.d0*ControlFlowGlobal%Dt) 
						! Call IonIonRecombinationRate(ControlFlowGlobal,ParticleGlobal,FieldOneGlobal,2,3,50.d0*ControlFlowGlobal%Dt) 
      !            Call ElectronIonRecombinationRate(ControlFlowGlobal,ParticleGlobal,FieldOneGlobal,2,50.d0*ControlFlowGlobal%Dt) 
						! 
						!!Call IonIonRecombinationRateOld(ControlFlowGlobal%Nx,ParticleGlobal(2),ParticleGlobal(1),ParticleGlobal(1)%Weight,50*ControlFlowGlobal%Dt)
      !!            Call IonIonRecombinationRateOld(ControlFlowGlobal%Nx,ParticleGlobal(2),ParticleGlobal(3),ParticleGlobal(3)%Weight,50*ControlFlowGlobal%Dt)
      !!            Call ElectronIonRecombinationRateOld(ControlFlowGlobal%Nx,FieldOneGlobal(0),ParticleGlobal(2),ParticleGlobal(0),ParticleGlobal(0)%Weight,50*ControlFlowGlobal%Dt) 
      !          end if
              return
	 End subroutine OneStep


        Subroutine OneStepRestart()
               !Use FileIO 
               Implicit none
               Integer(4) :: i
               !Write(*,*) Period,ParticleGlobal%NPar,"Period"!,FieldBoundaryGlobal%Qmin,FieldBoundaryGlobal%Qmax
               Call DumpFieldSolver(FieldSolverGlobal,0)
               Call DumpField(FieldGlobal,0)
               Call DumpFieldOne(ControlFlowGlobal%Ns,FieldOneGlobal,0)
               do i=0,ControlFlowGlobal%Ns
                     Call ParticleGlobal(i)%Dump(0)
               End do
               Call FieldBoundayFinalization(FieldBoundaryGlobal)
              return
        End  subroutine OneStepRestart
 
    End Module ModuleOneStep
    
    
