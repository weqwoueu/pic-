Module ModuleDiagOneStep
    Use ModuleControlFlow
    Use DiagnosticsCollisionRate
    Use DiagnosticsEEPF
    Use DiagnosticsMomentum
    !Use ModuleOneStepField
    !Use ModuleOneStep
    Use ModuleTrackingParticle
    use DiagnosticsTestParticle
    use DiagParticleBouncary
    Use ModuleMCCInterface
    Implicit none
  !    Type(Grid1D(Nx=NxMax,Ns=30)),save :: G1DDiagParticleField
  !    Type(Grid2D(Nx=NxMax,Ny=100,Ns=30)),save :: G2DDiagParticleField
  !    
		!Type(Grid1D(Nx=1,Ns=16)),save :: G1DDiagParticleBD
  !    Type(Grid2D(Nx=1,Ny=100,Ns=16)),save :: G2DDiagParticleBD 
  !
  !    Type(Grid1D(Nx=1000,Ns=6)),save :: G1DDiagParticleEEDF
  !    Type(Grid2D(Nx=1000,Ny=100,Ns=6)),save :: G2DDiagParticleEEDF
		!Type(Grid1D(Nx=1000,Ns=6)),save :: G1DDiagParticleIEDF
  !    Type(Grid2D(Nx=1000,Ny=100,Ns=6)),save :: G2DDiagParticleIEDF
      
    Type(Grid1D(Nx=1000,Ns=2)),save,Allocatable :: G1DDiagParticleEPF(:)
    !Type(Grid1D(Nx=1000,Ns=6)),save :: G1DDiagParticleIEPF
		
		!Type(Grid1D(Nx=NxMax,Ns=16)),save :: G1DDiagParticleCREG  !electron - gas    
  !    Type(Grid2D(Nx=NxMax,Ny=100,Ns=16)),save :: G2DDiagParticleCREG
  !    Type(Grid1D(Nx=NxMax,Ns=4)),save :: G1DDiagParticleCRIG  !Ion - gas
  !    Type(Grid2D(Nx=NxMax,Ny=100,Ns=4)),save :: G2DDiagParticleCRIG
  !    !Tracing Particle Information in Diagnosis OneStep
  !    type(Tracking_Particle_Information),save :: Tracking_Particle_Information1
  !    type(Tracking_Particle_Temp) :: TP_Temp_Global
      
      !Type(ParticleElectrode),save :: ParticleElectrodeGlobal(0:NsMax)
      
      contains
      Subroutine DiagInitilalization(CF)
            Implicit none
            Class(ControlFlow), intent(in) :: CF
            Integer(4) :: i
            !Integer(4),intent(in) ::  Ns
            !Type(ParticleBundle),intent(in) :: PB(0:Ns)
                 !Call G1DDiagParticleField%Init(CF)
                 !Call G2DDiagParticleField%Init(CF)
					  !call G1DDiagParticleBD%Init(CF)
					  !call G2DDiagParticleBD%Init(CF)
       !
       !          Call G1DDiagParticleEEDF%Init(CF)
       !          Call G2DDiagParticleEEDF%Init(CF)
                 Allocate(G1DDiagParticleEPF(0:CF%Ns))
                 Do i = 0,CF%Ns
                    Select case(i)
                    Case(0)
                        Call G1DDiagParticleEPF(i)%Init(CF,'ElectronEPF',0.1_8)       
                    Case(1)
                        Call G1DDiagParticleEPF(i)%Init(CF,'IonEPF',1._8)  
                    End Select
                 EndDo
                 
					  !Call G1DDiagParticleIEDF%Init(CF)
       !          Call G2DDiagParticleIEDF%Init(CF)
					  !
					  !Call G1DDiagParticleCREG%Init(CF)
       !          Call G2DDiagParticleCREG%Init(CF)
					  !Call G1DDiagParticleCRIG%Init(CF)
       !          Call G2DDiagParticleCRIG%Init(CF)
                 !Call Diag_Tracking_Partcile_Initilalization(CF,ParticleGlobal,Tracking_Particle_Information1,TP_Temp_Global)          
          return  
      End Subroutine DiagInitilalization
      
    Subroutine  DiagOneStep(CF)
        Implicit none 
        Class(ControlFlow), intent(in) :: CF
        Integer(4) :: i,j

             !Call DiagParticleFieldPeriod(G1DDiagParticleField,ControlFlowGlobal%Ns,ParticleGlobal,FieldGlobal,0)
             !Call DiagParticleFieldPeriod(G2DDiagParticleField,ControlFlowGlobal%Ns,ParticleGlobal,FieldGlobal,0)
				 !call DiagParticleBoundaryOne(G1DDiagParticleBD,ParticleBDOneGlobal,ControlFlowGlobal,0)
				 !call DiagParticleBoundaryOne(G2DDiagParticleBD,ParticleBDOneGlobal,ControlFlowGlobal,0)
     !        Call DiagParticleEDFOne(G1DDiagParticleEEDF,ParticleGlobal(0),ParticleBDOneGlobal(0),0)
     !        Call DiagParticleEDFOne(G2DDiagParticleEEDF,ParticleGlobal(0),ParticleBDOneGlobal(0),0)
        Do i = 0,CF%Ns
            Call DiagParticleEPFOne(G1DDiagParticleEPF(i),ParticleGlobal(i),0)
        End Do
				 !Call DiagParticleEDFOne(G1DDiagParticleIEDF,ParticleGlobal(2),ParticleBDOneGlobal(2),0)
     !        Call DiagParticleEDFOne(G2DDiagParticleIEDF,ParticleGlobal(2),ParticleBDOneGlobal(2),0)
				 !Call DiagParticleCollisionRateOne(G1DDiagParticleCREG,ControlFlowGlobal,ParticleGlobal(0),MCCBundleGlobal(0,1),0)
     !        Call DiagParticleCollisionRateOne(G2DDiagParticleCREG,ControlFlowGlobal,ParticleGlobal(0),MCCBundleGlobal(0,1),0)
				 !Call DiagParticleCollisionRateOne(G1DDiagParticleCRIG,ControlFlowGlobal,ParticleGlobal(1),MCCBundleGlobal(1,1),0)
     !        Call DiagParticleCollisionRateOne(G2DDiagParticleCRIG,ControlFlowGlobal,ParticleGlobal(1),MCCBundleGlobal(1,1),0)
             !call Tracking_Particle(ParticleGlobal,Tracking_Particle_Information1,ControlFlowGlobal,i,j,TP_Temp_Global)
             !call Diag_Particle_Test_In_Field(ParticleGlobal(0),ParticleBDOneGlobal(0),ControlFlowGlobal,i,j,FieldGlobal)
        return  
      End Subroutine DiagOneStep

       !
       Subroutine DiagOneStepFinal(CF)
            Implicit none 
            Class(ControlFlow), intent(in) :: CF
            Integer(4) :: i
             !Call DiagParticleFieldPeriod(G1DDiagParticleField,ControlFlowGlobal%Ns,ParticleGlobal,FieldGlobal,1)
             !Call DiagParticleFieldPeriod(G2DDiagParticleField,ControlFlowGlobal%Ns,ParticleGlobal,FieldGlobal,1)
				 !call DiagParticleBoundaryOne(G1DDiagParticleBD,ParticleBDOneGlobal,ControlFlowGlobal,1)
				 !call DiagParticleBoundaryOne(G2DDiagParticleBD,ParticleBDOneGlobal,ControlFlowGlobal,1)
     !        Call DiagParticleEDFOne(G1DDiagParticleEEDF,ParticleGlobal(0),ParticleBDOneGlobal(0),1)
     !        Call DiagParticleEDFOne(G2DDiagParticleEEDF,ParticleGlobal(0),ParticleBDOneGlobal(0),1)
				 !Call DiagParticleEDFOne(G1DDiagParticleIEDF,ParticleGlobal(2),ParticleBDOneGlobal(2),1)
     !        Call DiagParticleEDFOne(G2DDiagParticleIEDF,ParticleGlobal(2),ParticleBDOneGlobal(2),1)
     !        Call DiagParticleCollisionRateOne(G1DDiagParticleCREG,ControlFlowGlobal,ParticleGlobal(0),MCCBundleGlobal(0,1),1)
     !        Call DiagParticleCollisionRateOne(G2DDiagParticleCREG,ControlFlowGlobal,ParticleGlobal(0),MCCBundleGlobal(0,1),1)
				 !Call DiagParticleCollisionRateOne(G1DDiagParticleCRIG,ControlFlowGlobal,ParticleGlobal(1),MCCBundleGlobal(1,1),1)
				 !Call DiagParticleCollisionRateOne(G2DDiagParticleCRIG,ControlFlowGlobal,ParticleGlobal(1),MCCBundleGlobal(1,1),1)
            Do i = 0,CF%Ns
                Call DiagParticleEPFOne(G1DDiagParticleEPF(i),ParticleGlobal(i),1)
            End Do
          return  
       End Subroutine DiagOneStepFinal
       
      Subroutine DiagOneStepFinal2()
       Implicit none  
          !Call DiagParticleTestParticle(ParticleGlobal(0),ParticleBDOneGlobal(0),FieldGlobal,1)              
          return  
        End Subroutine DiagOneStepFinal2

  End Module ModuleDiagOneStep