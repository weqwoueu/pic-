Module ModuleSpecyOne
     !Use TypeModule
     Use Constants
     Use ModuleControlFlow
     Implicit none
     Integer(4),parameter,private :: NsMax=9
     Integer(4),parameter :: NgMax=3 
     
     Type  SpecyOnePhysical      ! Define a physical particle.
                       Character(99) :: Name="Ar+"
                       Integer(4) :: NAtom=1,Charge=0   !MCModel=1,
                       Real(8) :: Mass=0.d0,Radius=1.d0,polarizability=0.d0!,Test(2,2)=0.d0!,test2(3)=1.d0
     End  Type  SpecyOnePhysical
     
    Type  SpecyOnePegasus     ! Define a physical particle.
        Integer(4) :: Index=0                 !粒子编号
        Real(8) :: Mass=0.d0     !粒子质量(AMU单位)
        Integer(4) :: Charge=0   !粒子单位电荷
        Real(8) :: Radius=1.d0  !粒子半径（pm单位）
        Real(8) :: Polarizability=0.d0 !极化率(化学反应性碰撞，单位为Bohr半径的倍数)
        Character(CharLenthMax) :: Name="Ar+" !粒子符号
        contains
        Procedure :: Load=>LoadSpecyOnePegasus
        Procedure :: Convert2Specy=>ConvertSpecyOnePegasus2SpecyOnePhysical
        Procedure :: Convert2Gas=>ConvertSpecyOnePegasus2GasPhysical
        Generic :: Convert=>Convert2Specy,Convert2Gas
    EndType SpecyOnePegasus 
     
     Type(SpecyOnePhysical),parameter :: ElectronPhysical=SpecyOnePhysical('Electron',1,-1,ElectronMass/AtomicMass,0.d0) 

    Type GasPhysical    ! Define a physical gas.('Electron',1,-1,9.1095d-31)
        Character(99) :: Name="Ar"
        ! CSModel is the crosssection data source index: 1-lxcat; 2-pegasus 
        Integer(4) ::  MCModel=1,NAtom,Ns,CSModel
        Real(8) ::  Mass,Radius=0.d0,BetaMax=0.d0,polarizability=0.d0
        Type (SpecyOnePhysical) :: SP(NsMax)
        Type (SpecyOnePegasus) :: SOP(0:NsMax)
    End Type  GasPhysical

    ! Model is the background gas mode,
    !   1X the density and temperature will not change over time:
    !      10  --the density and temperature are both constant over space, default; 
    !      11  --the density and temperature are both NOT constant over space ;
    !   2X Neautral  depleting is considered, but the energy balance or heating is not
    !   3X Plasma+continous Fluid Gas
   !    4X Plasma + DSMC Gas

    Type ::  SpecyOne     ! Define a physical particle.
                Character(99) :: Name
                Integer(4) :: SpecyIndex=-1,GasIndex=-1
                Real(8) :: Charge,Mass,Radius,NAtom
                Real(8) :: InitDensity,Density,InitTemperature,Temperature
    End  Type  SpecyOne
    
    Type ::  GasOne    ! Define a physical particle.
                Character(99) :: Name
                Integer(4) ::  MCModel=-1,Ns,IndexStart,CSModel=0
                Real(8) :: Mass,Radius,NAtom,BetaMax,polarizability
                Real(8) :: InitDensity,Density,InitTemperature,Temperature
                Type (SpecyOnePegasus) :: SOP(0:NsMax)
    End  Type  GasOne
    



     Type(SpecyOne),Allocatable,save :: SpecyGlobal(:)
     Type(GasOne),Allocatable,save :: GasGlobal(:)
     
     contains
             Subroutine  GasInit(CF)!Gas)
                  Implicit none
                  !Type(GasOne),Intent(inout) :: Gas
                  !Type(GasOne) :: GO
                  Class(ControlFlow), intent(inout) :: CF
                  Type(GasPhysical) :: GP
                  Type(GasPhysical),Allocatable :: TempGP(:)
                  Character(len=99),dimension(NgMax) :: GasName
                  Integer(4) :: Ng=0,Ns=0
                  Real(8) :: PGas,TGas,TElectron,GasRatio(NgMax)

                  Integer(4) ::i,j,k,Index=0
                  logical :: Alive

  
                  Character(len=99) :: Filename
                  NAMELIST /Ngfile/ Ng,GasName,PGas,TGas,TElectron,GasRatio
                  NAMELIST /GasPhysicsfile/ GP
                  !NAMELIST /Gasfile/ GasGlobal
                  
                  Filename="./input/gas.txt"
                  Inquire(file=Filename,exist=alive)
                   If(alive) then
                       OPEN(20,FILE=Filename)
                       Read(20,NML=Ngfile)
                       Close (20)
                   ENd If
                   
                   CF%Ng=Ng
                   Allocate(GasGlobal(Ng))
                   Allocate(TempGP(Ng))
                   
                   do i=1, Ng
                          Filename="./gas/"//Trim(GasName(i))//"/"//Trim(GasName(i))//".txt"
                          Inquire(file=Filename,exist=alive)
                          If(alive) then
                               OPEN(20,FILE=Filename)
                               Read(20,NML=GasPhysicsfile)
                               Close (20)
                          Else

                          End If
                          TempGP(i)=GP
                          Ns=Ns+GP%Ns
                   End Do
                   
                   
                   CF%Ns=Ns
                   Allocate(SpecyGlobal(0:Ns))
                   
                   Call SpecyOnePhysics2SpecyOne(SpecyGlobal(0),ElectronPhysical)
                   SpecyGlobal(0)%SpecyIndex=0
                   SpecyGlobal(0)%GasIndex=0
                   
                   SpecyGlobal(0)%InitDensity=1.d0
                   SpecyGlobal(0)%Density=SpecyGlobal(0)%InitDensity
                   SpecyGlobal(0)%InitTemperature = TElectron
                   SpecyGlobal(0)%Temperature=SpecyGlobal(0)%InitTemperature

                   Do i=1,Ng
                         Call GasPhysics2Gas(GasGlobal(i),TempGP(i))

                         GasGlobal(i)%InitTemperature=TGas
                          GasGlobal(i)%Temperature=GasGlobal(i)%InitTemperature
                          GasGlobal(i)%InitDensity=(PGas*mTorrtoPa)/(kB*GasGlobal(i)%InitTemperature)*GasRatio(i)
                          GasGlobal(i)%Density=GasGlobal(i)%InitDensity
                          GasGlobal(i)%IndexStart=Index
                          
                         Do j=1,GasGlobal(i)%Ns
                              Index=Index+1
                              SpecyGlobal(Index)%SpecyIndex=Index
                              SpecyGlobal(Index)%GasIndex=i
                              SpecyGlobal(Index)%InitDensity=GasRatio(i)
                              SpecyGlobal(Index)%Density = SpecyGlobal(Index)%InitDensity
                              SpecyGlobal(Index)%InitTemperature = TGas
                              SpecyGlobal(Index)%Temperature=SpecyGlobal(Index)%InitTemperature
                              Call SpecyOnePhysics2SpecyOne(SpecyGlobal(Index),TempGP(i)%SP(j))
                         End Do
                   ENd do
                   DeAllocate(TempGP)
                   return
             End  subroutine  GasInit
             
             Subroutine  SpecyOnePhysics2SpecyOne(SO,SP)
                  Implicit none
                  Type(SpecyOne),Intent(inout)  :: SO
                  Type(SpecyOnePhysical),Intent(in)  :: SP
                  SO%Name=SP%Name
                  !SO%MCModel=SP%MCModel
                  SO%Charge=Dble(SP%Charge)*ElectronCharge
                  SO%Mass=SP%Mass*AtomicMass    ! reletive mass(base: 1/12 of carbon atom) -> real mass
                  SO%NAtom=Dble(SP%NAtom)
                  SO%Radius=SP%Radius
                  return
             End  subroutine  SpecyOnePhysics2SpecyOne

            Subroutine  GasPhysics2GasPegasus(GO,GP)
                  Implicit none
                  !Type(GasOne),Intent(inout) :: Gas
                  Type(GasOne),Intent(inout)  :: GO
                  Type(GasPhysical),Intent(in)  :: GP
                  Integer(4) ::  i
                          GO%Name=GP%Name
                          GO%CSModel=GP%CSModel
                          GO%MCModel=GP%MCModel
                          GO%Ns=GP%Ns
                          GO%Mass=GP%Mass*AtomicMass
                          GO%NAtom=Dble(GP%NAtom)
                          GO%Radius=GP%Radius
                          GO%BetaMax=GP%BetaMax
                          GO%polarizability=14.d0*a0*a0*a0
                          GO%SOP=GP%SOP
                    return
            End  subroutine  GasPhysics2GasPegasus
            
            Subroutine  GasPhysics2Gas(GO,GP)
                  Implicit none
                  !Type(GasOne),Intent(inout) :: Gas
                  Type(GasOne),Intent(inout)  :: GO
                  Type(GasPhysical),Intent(in)  :: GP
                  Integer(4) ::  i
                          GO%Name=GP%Name
                          GO%CSModel=GP%CSModel
                          GO%MCModel=GP%MCModel
                          GO%Ns=GP%Ns
                          GO%Mass=GP%Mass*AtomicMass
                          GO%NAtom=Dble(GP%NAtom)
                          GO%Radius=GP%Radius
                          GO%BetaMax=GP%BetaMax
                          GO%polarizability=14.d0*a0*a0*a0
                    return
            End  subroutine  GasPhysics2Gas
            
            subroutine ConvertSpecyOnePegasus2SpecyOnePhysical(SOP,SP)
                Class(SpecyOnePegasus),intent(inout) :: SOP
                Type(SpecyOnePhysical),Intent(inout)  :: SP
                !SP%Index=SOP%Index
                SP%Mass=SOP%Mass
                SP%Charge=SOP%Charge
                SP%Radius=SOP%Radius*1.d-10
                SP%Polarizability=SOP%Polarizability
                SP%Name=SOP%Name
                
                return
        end subroutine ConvertSpecyOnePegasus2SpecyOnePhysical
        
        subroutine ConvertSpecyOnePegasus2GasPhysical(SOP,GP)
                Class(SpecyOnePegasus),intent(inout) :: SOP
                Type(GasPhysical) :: GP
                !SP%Index=SOP%IndexReal(8) ::  Mass,Radius=0.d0,BetaMax=0.d0,polarizability=0.d0
                GP%Name=SOP%Name
                GP%Mass=SOP%Mass
                GP%Radius=SOP%Radius*1.d-10
                return
        end subroutine ConvertSpecyOnePegasus2GasPhysical
        
        subroutine LoadSpecyOnePegasus(SOP,IONumber)
                Class(SpecyOnePegasus),intent(inout) :: SOP
                Integer(4) :: IONumber
                Read(IONumber,*) SOP%Index,SOP%Mass,SOP%Charge,SOP%Radius,SOP%Polarizability,SOP%Name
                return
        end subroutine LoadSpecyOnePegasus
  End Module  ModuleSpecyOne