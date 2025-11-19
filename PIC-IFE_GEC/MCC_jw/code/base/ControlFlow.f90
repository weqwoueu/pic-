Module ModuleControlFlow
     Use Constants
     Implicit none
            Integer(4),Parameter :: NxMax=512
            Real(8),parameter :: ZLength=0.02d0
            Real(8),parameter :: Inputdx=ZLength/dble(NxMax-1)
            !Real(8),parameter :: Inputdt=1.25d-11 !4.d-12
            !Real(8),parameter :: Inputdt=5d-11 !4.d-12
            !Real(8),parameter :: Inputdt=1.0d-10 !4.d-12
            Real(8),parameter :: Inputdt=2.0d-10 !4.d-12
              
            Integer(4),Parameter  ::  DefaultNameIndex=10000   !Index
            Integer(4),Parameter  ::  DefaultNameIndexInit=20000
            Integer(4),Parameter  ::  ModeMultiplier=1000
				Integer(4),Parameter  ::  NspacyMax=9
 
     !  Ns--NSpecy,Ng---NGas
     Type ControlFlow
           Real(8)  :: Dx=Inputdx,Dt=Inputdt
           Integer(4) :: ParticlePerGrid=1000
           Real(8)  :: InitDensity=1.d15
           !Specify how many types of particles are , value 0 indicate that ion type is 0 , no ions
			  Real(8)  :: ElectrodeArea = 1.0d0
			  Real(8)  :: ElectrodeArea2 = 2.0d0
           Integer(4) :: Ns=2,Ng=0,Nt=50 !Ns => N species type of particle
           Integer(4) :: Nx=NxMax,NxL=0,NxU=NxMax-1
           Integer(4) :: Timer=0,Period=0
           Integer(4) :: NRun=1000,NDiagShort=10,NDiagLong=100  !Some diagosis needs lower Period but some needs more
           Logical :: ReStartParticles=.False.
           Logical :: withRecombination=.true.
           !contains
              !procedure :: Init=>InitializationControlFlow
           End Type ControlFlow
       contains     
        Subroutine InitializationControlFlow(CF)
               Type(ControlFlow) :: CF
                Logical :: Alive
                Character(len=99) :: Filename
                NAMELIST /ControlFlow/ CF
                !Filename="./input/controlflow.txt"
                Filename="./MCC_jw/input/controlflow.txt"   !$ mb.ZWZ 2021/12/19
                  Inquire(file=Filename,exist=alive)
                   If(alive) then
                       OPEN(20,FILE=Filename)           !$InitDensity, NRun, NDiagShort, NDiagLong, ParticlePerGrid
                       Read(20,NML=ControlFlow)
                       Close (20)
                   Else
                      Write(*,*)  "ControlFlow Load", Trim(Filename),"ERROR! The program will abort!"
                      !Stop
                   ENd If
        End subroutine InitializationControlFlow
End Module ModuleControlFlow