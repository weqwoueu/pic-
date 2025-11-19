      Module ModuleParticleOne
        Use ModuleField 
        Implicit none
        Type :: ParticleOne
		       Real(8) :: X,Y,Z,Vx,Vy,Vz,Ax,Ay,Az,Wq
               Integer :: Location ! the element where the particle is located
        contains
          procedure :: ParLocate => ParticleOneLocalization
        
		  procedure :: PosInit=>PositionRandomInitializationParticleOne  !"ParticleOne%PosInit(   )" to use
          
          procedure :: VelInpInit=>VelocityInputInitializationParticleOne
          procedure :: VelMaxInit=>VelocityMaxwellianInitializationParticleOne
          procedure :: VelRanInit=>VelocityRandomInitializationParticleOne
 
          procedure :: AccInpInit=>AccelerationInputInitializationParticleOne
          
          !Convert in different unit systems
          procedure :: PosRes=>PositionRescaleParticleOne
          procedure :: VelRes=>VelocityRescaleParticleOne
          procedure :: AccRes=>AccelerationRescaleParticleOne
          !to caculate one particle energy
           procedure :: Energy=>EnergyParticleOne
          
          procedure :: Copy=>CopyParticleOne
          procedure :: Swap=>SwapParticleOne
          ! value cell to position
          procedure :: WeightP2C=>WeightP2CParticleOne
			 procedure :: WeightP2C2=>WeightP2CParticleOne2
          !Interpolate the field amount to the particle position based on the grid point position
          procedure :: WeightC2PES=>WeightC2PElectrostaticParticleOne
          !two kind ways to push paritcle
          procedure :: MoveES=>MoveElectrostaticParticleOne
			  procedure :: MoveES1st=>MoveElectrostaticParticleOne1st
			  procedure :: MoveES2nd=>MoveElectrostaticParticleOne2nd
			 procedure :: ExMoveES=>ExMoveElectrostaticParticleOne


          procedure :: WeightC2PEM=>WeightC2PElectromagneticParticleOne
          procedure :: MoveEM=>MoveElectromagneticParticleOne
      End Type ParticleOne
     
      Type,extends(ParticleOne) :: ParticleOneIndex
          Integer(4) ::  Index    !Mark each particle
          contains
          procedure :: IndexInit=>IndexInitializationParticleOneIndex
          procedure :: Copy=>CopyParticleOneIndex
          procedure :: Swap=>SwapParticleOneIndex
       End Type ParticleOneIndex
     
    contains
!DIR$ ATTRIBUTES FORCEINLINE ::  PositionRandomInitializationParticleOne
!DIR$ ATTRIBUTES FORCEINLINE ::  VelocityInputInitializationParticleOne,VelocityMaxwellianInitializationParticleOne,VelocityRandomInitializationParticleOne
!DIR$ ATTRIBUTES FORCEINLINE ::  AccelerationInputInitializationParticleOne
!DIR$ ATTRIBUTES FORCEINLINE ::  PositionRescaleParticleOne,  VelocityRescaleParticleOne,AccelerationRescaleParticleOne
!DIR$ ATTRIBUTES FORCEINLINE ::  CopyParticleOne,EnergyParticleOne
!DIR$ ATTRIBUTES FORCEINLINE ::   WeightP2CParticleOne
!DIR$ ATTRIBUTES FORCEINLINE ::  WeightC2PElectrostaticParticleOne,MoveElectrostaticParticleOne
!DIR$ ATTRIBUTES FORCEINLINE ::   WeightC2PElectromagneticParticleOne,MoveElectromagneticParticleOne 
        
         subroutine ParticleOneLocalization(PO, Index)
            Class(ParticleOne), intent(inout) :: PO    !abbreviate particleone to PO,use class because of extends
            Integer, intent(in) :: Index
                PO%Location = Index
         End subroutine
    
         subroutine PositionRandomInitializationParticleOne(PO,XL,XU,YL,YU,ZL,ZU,Radius_L,Radius_U,Theta_L,Theta_U)
            Use Domain_2D, Only: delta_global ,region_type
            Class(ParticleOne), intent(inout) :: PO    !abbreviate particleone to PO,use class because of extends
            Real(8),intent(in),optional :: XL,XU,YL,YU,ZL,ZU,Radius_L,Radius_U,Theta_L,Theta_U
            real (8):: Radius,Theta,a,b,c
            
            
            a=XU
            b=3.2768
            if (region_type==0) then 
            
             CALL RANDOM_NUMBER(R)
            PO%X=XL+(XU-XL)*R   !XL=> X,Lower
            If (delta_global == 0) Then
                CALL RANDOM_NUMBER(R)
                c=0.0
                c=YU*(1.0d0 - EXP**(-(((PO%X - a)**2) / b)))
                PO%Y=YL+(YU-YL)*R
                PO%Z=0.
                
                !if(c<PO%Y)then 
                !    PO%X=2000
                !end if
                
                
            Elseif (delta_global == 1) Then
                CALL RANDOM_NUMBER(R)
                PO%Y=SQRT(YL**2+(YU+YL)*(YU-YL)*R)
                CALL RANDOM_NUMBER(R)
                PO%Z=2*PI*R
            Endif
            
         elseif (region_type==1)then 
             
         CALL RANDOM_NUMBER(R4)
         CALL RANDOM_NUMBER(R5)
         
         
         Radius=SQRT (Radius_L*Radius_L+(Radius_U*Radius_U-Radius_L*Radius_L)*R4)
         Theta=Theta_L+(Theta_U-Theta_L)*R5
         
         PO%X=Radius*COS (Theta)
         PO%Y=Radius*SIN (Theta)
         PO%Z=0.
             
            endif 
             
             
            
         End subroutine
         
        


    
         subroutine VelocityInputInitializationParticleOne(PO,Vx,Vy,Vz)
		    Class(ParticleOne), intent(inout) :: PO
            Real(8),intent(in),optional :: Vx,Vy,Vz
               PO%Vx=Vx
               PO%Vy=Vy
               PO%Vz=Vz
         end subroutine VelocityInputInitializationParticleOne
         
         Subroutine VelocityMaxwellianInitializationParticleOne(PO,Mass,Temperature)
		    Class(ParticleOne), intent(inout) :: PO
            Real(8),intent(in) :: Mass,Temperature
            Real(8) :: V,Beta,FuncA,FuncB
            Real(8) :: Theta,CosTheta,SinTheta,Fai 
            Beta=1.d0/(kB*Temperature)
            FuncA=1.d0
            FuncB=0.d0
            do while(FuncA>FuncB)
                CALL RANDOM_NUMBER(R)
                    FuncA=R*R
                CALL RANDOM_NUMBER(R)
                    FuncB=-exp*R*Dlog(R)
            end do
            V=DSqrt(-3.d0*Dlog(R)/Beta/Mass)
            Call VelocityRandomInitializationParticleOne(PO,V)  !to get the random direction
            return 
         End subroutine VelocityMaxwellianInitializationParticleOne
         
         Subroutine VelocityRandomInitializationParticleOne(PO,V)
		    Class(ParticleOne), intent(inout) :: PO
               Real(8),intent(in) ::  V
               Real(8) :: Fai,CosFai,SinFai
               Real(8) :: Theta,CosTheta,FcosTheta,SinTheta
               CALL RANDOM_NUMBER(R)
               CosTheta=1.d0-2.d0*R
                SinTheta=Dsqrt(1.d0-cosTheta*cosTheta)
                CALL RANDOM_NUMBER(R)
                Fai=2.d0*PI*R
                 PO%Vx=V*CosTheta
                PO%Vy=V*SinTheta*DCos(Fai)
                PO%Vz=V*SinTheta*Dsin(Fai)
               return
         end subroutine VelocityRandomInitializationParticleOne
         
         Subroutine AccelerationInputInitializationParticleOne(PO,Ax,Ay,Az)
            Class(ParticleOne), intent(inout) :: PO
            Real(8),intent(in),optional :: Ax,Ay,Az
            If (Present(Ax)) Then
                     PO%Ax=Ax
                Else
                     PO%Ax=0.d0
                End If
                If (Present(Ay)) Then
                     PO%Ay=Ay
                Else
                     PO%Ay=0.d0
                End If
                If (Present(Az)) Then
                     PO%Az=Az
                Else
                     PO%Az=0.d0
                End If
         End subroutine  AccelerationInputInitializationParticleOne
         
         Subroutine PositionRescaleParticleOne(PO,XFactor)
             Class(ParticleOne), intent(inout) :: PO
             Real(8),intent(in) :: XFactor
             PO%X=PO%X*XFactor 
         End subroutine PositionRescaleParticleOne

         Subroutine VelocityRescaleParticleOne(PO,VFactor)
             Class(ParticleOne), intent(inout) :: PO
             Real(8),intent(in) :: VFactor
             PO%Vx=PO%Vx*VFactor
             PO%Vy=PO%Vy*VFactor
             PO%Vz=PO%Vz*VFactor
         End subroutine VelocityRescaleParticleOne
         
         Subroutine AccelerationRescaleParticleOne(PO,AFactor)
             Class(ParticleOne), intent(inout) :: PO
             Real(8),intent(in) :: AFactor
             PO%Vx=PO%Ax*AFactor
             PO%Vy=PO%Ay*AFactor
             PO%Vz=PO%Az*AFactor
         End subroutine AccelerationRescaleParticleOne

         subroutine CopyParticleOne(POD,POC)
            Class(ParticleOne), intent(inout) :: POD
            Class(ParticleOne), intent(in) :: POC
            !INteger(4),optional :: Ndel
            Select Type (POD)
                 Type is (ParticleOne)
                      Select Type (POC)
                          Type is (ParticleOne)
                          POD=POC
                     ENd Select
            ENd Select
         End subroutine  CopyParticleOne
         
        subroutine SwapParticleOne(POD,POC)
            Class(ParticleOne), intent(inout) :: POD
            Class(ParticleOne), intent(inout) :: POC
            Type (ParticleOne) :: POT
            Select Type (POD)
                 Type is (ParticleOne)
                      Select Type (POC)
                        Type is (ParticleOne)
                          POT=POC
                          POC=POD
                          POD=POT
                     ENd Select
            ENd Select
         End subroutine  SwapParticleOne

         Function EnergyParticleOne(PO,Mass,VFactor) 
              Implicit none
              Class(ParticleOne),intent(in) :: PO
              Real(8),intent(in) :: Mass
              Real(8),intent(in),optional :: VFactor
              Real(8) :: EnergyParticleOne
              EnergyParticleOne=0.5d0*Mass*(PO%Vx*PO%Vx+PO%Vy*PO%Vy +PO%Vz*PO%Vz)
              If (Present(VFactor)) EnergyParticleOne=EnergyParticleOne*VFactor*VFactor
             return 
         End  Function EnergyParticleOne
         
         Subroutine WeightP2CParticleOne(PO,FO)
            Class(ParticleOne), intent(in) :: PO
            Type(FieldOne), intent(inout) :: FO
            Integer(4) :: Nx
            Real(8) :: S1,S2
                    Nx=Ceiling(PO%X)
                    S1=Dble(Nx)-PO%X
                    S2=1.d0-S1
                    FO%RhoOne(Nx)=FO%RhoOne(Nx)+S1
                    FO%RhoOne(Nx+1)= FO%RhoOne(Nx+1)+S2
			End subroutine  WeightP2CParticleOne
			
			Subroutine WeightP2CParticleOne2(PO,FO)
            Class(ParticleOne), intent(in) :: PO
            Type(FieldOne), intent(inout) :: FO
            Integer(4) :: Nx
            Real(8) :: S1,S2,Energy
                    Nx=Ceiling(PO%X)
                    S1=Dble(Nx)-PO%X
                    S2=1.d0-S1
                    FO%RhoOne(Nx)=FO%RhoOne(Nx)+S1
                    FO%RhoOne(Nx+1)= FO%RhoOne(Nx+1)+S2
						  Energy = PO%Vx*PO%Vx+PO%Vy*PO%Vy+PO%Vz*PO%Vz
						  FO%EnergyOne(Nx)=FO%EnergyOne(Nx)+S1*Energy
                    FO%EnergyOne(Nx+1)=FO%EnergyOne(Nx+1)+S2*Energy
						  !FO%JxOne(Nx)=FO%JxOne(Nx)+S1*PO%Vx
						  !FO%JxOne(Nx+1)=FO%JxOne(Nx+1)+S2*PO%Vx
			End subroutine  WeightP2CParticleOne2
			
         
         subroutine WeightC2PElectrostaticParticleOne(PO,FG,Ex)
            Class(ParticleOne), intent(inout) :: PO
            Type(Field), intent(in) :: FG
            Real(8),intent(out) :: Ex
            Integer(4) :: Nx
            Real(8) :: S1,S2
                    Nx=Ceiling(PO%X)
                    S1=Dble(Nx)-PO%X
                    S2=1.d0-S1
                    Ex=FG%Ex(Nx)*S1+FG%Ex(Nx+1)*S2
         End subroutine  WeightC2PElectrostaticParticleOne
         
         subroutine WeightC2PElectromagneticParticleOne(PO,FG,Ex,Ey,Ez,Bx,By,Bz)
            Class(ParticleOne), intent(inout) :: PO
            Type(Field), intent(in) :: FG
            Real(8),intent(out) :: Ex,Ey,Ez,Bx,By,Bz
            Integer(4) :: Nx
            Real(8) :: S1,S2
                    Nx=Ceiling(PO%X)
                    S1=Dble(Nx)-PO%X
                    S2=1.d0-S1
                    Ex=FG%Ex(Nx)*S1+FG%Ex(Nx+1)*S2
                    Ey=FG%Ey(Nx)*S1+FG%Ey(Nx+1)*S2
                    Ez=FG%Ez(Nx)*S1+FG%Ez(Nx+1)*S2
                    Bx=FG%Bx(Nx)*S1+FG%Bx(Nx+1)*S2
                    By=FG%By(Nx)*S1+FG%By(Nx+1)*S2
                    Bz=FG%Bz(Nx)*S1+FG%Bz(Nx+1)*S2
         End subroutine  WeightC2PElectromagneticParticleOne

         Subroutine MoveElectrostaticParticleOne(PO,Ex)
            Class(ParticleOne), intent(inout) :: PO
            Real(8),intent(in),optional :: Ex
                PO%Vx=PO%Vx+0.5d0*Ex
                PO%X=PO%X+0.5d0*Ex
                PO%Ax=0.5d0*(PO%Ax+Ex)
                PO%Vx=PO%Vx+0.5d0*PO%Ax
                PO%X=PO%X+PO%Vx
			End subroutine  MoveElectrostaticParticleOne
			
			Subroutine MoveElectrostaticParticleOne1st(PO)
            Class(ParticleOne), intent(inout) :: PO
                !PO%Vx=PO%Vx+0.5d0*Ex
                !PO%X=PO%X+0.5d0*Ex
                !PO%Ax=0.5d0*(PO%Ax+Ex)
                PO%Vx=PO%Vx+0.5d0*PO%Ax
                PO%X=PO%X+PO%Vx
			End subroutine  MoveElectrostaticParticleOne1st
			
						
			Subroutine MoveElectrostaticParticleOne2nd(PO,Ex)
            Class(ParticleOne), intent(inout) :: PO
            Real(8),intent(in),optional :: Ex
				Real(8) :: S1,S2
				integer(4) :: Nx
                PO%Vx=PO%Vx+0.5d0*Ex
                PO%X=PO%X+0.5d0*Ex
                PO%Ax=0.5d0*(PO%Ax+Ex)
			End subroutine  MoveElectrostaticParticleOne2nd		
			
			Subroutine ExMoveElectrostaticParticleOne(PO,Ex)
            Class(ParticleOne), intent(inout) :: PO
            Real(8),intent(in),optional :: Ex
				!Type(ParticleBoundaryOne),intent(in) :: PBDO
				Real(8) :: S1,S2
				integer(4) :: Nx
                PO%Vx=PO%Vx+Ex
                PO%X=PO%X+PO%Vx
               RETURN
			End subroutine  ExMoveElectrostaticParticleOne 
         
         subroutine MoveElectromagneticParticleOne(PO,Ex,Ey,Ez,Bx,By,Bz)
            Class(ParticleOne), intent(inout) :: PO
            Real(8),intent(in) :: Ex,Ey,Ez,Bx,By,Bz
            Real(8) :: SVelocity(1:3)=0.d0
            Real(8) :: Omega(1:3)=0.d0,OmegaT,TransB(3,3)=0.d0
                       Omega(1)=Bx*0.5d0
                       Omega(2)=By*0.5d0
                       Omega(3)=Bz*0.5d0
                       OmegaT=1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
                       
                       TransB(1,1)=1.d0+Omega(1)*Omega(1)
                       TransB(2,1)=Omega(1)*Omega(2)+Omega(3)
                       TransB(3,1)=Omega(1)*Omega(3)-Omega(2)
                       TransB(1,2)=Omega(1)*Omega(2)-Omega(3)
                       TransB(2,2)=1.d0+Omega(2)*Omega(2)
                       TransB(3,2)=Omega(2)*Omega(3)+Omega(1)
                       TransB(1,3)=Omega(1)*Omega(3)+Omega(2)
                       TransB(2,3)=Omega(2)*Omega(3)-Omega(1)
                       TransB(3,3)=1.d0+Omega(3)*Omega(3)
                       TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT

                       SVelocity(1)=(TransB(1,1)*Ex+TransB(2,1)*Ey+TransB(3,1)*Ez)*0.5d0
                       SVelocity(2)=(TransB(1,2)*Ex+TransB(2,2)*Ey+TransB(3,2)*Ez)*0.5d0
                       SVelocity(3)=(TransB(1,3)*Ex+TransB(2,3)*Ey+TransB(3,3)*Ez)*0.5d0

                       PO%Vx=PO%Vx+SVelocity(1)
                       PO%Vy=PO%Vy+SVelocity(2)
                       PO%Vz=PO%Vz+SVelocity(3)
                       
                       PO%X=PO%X+SVelocity(1)
                       PO%Y=PO%Y+SVelocity(2)
                       PO%Z=PO%Z+SVelocity(3)

                       PO%Ax=0.5d0*(PO%Ax+Ex) 
                       PO%Ay=0.5d0*(PO%Ay+Ey)
                       PO%Az=0.5d0*(PO%Az+Ez)
                      
                       SVelocity(1)=PO%Vx+0.5d0*PO%Ax+PO%Vy*Omega(3)-PO%Vz*Omega(2)
                       SVelocity(2)=PO%Vy+0.5d0*PO%Ay+PO%Vz*Omega(1)-PO%Vx*Omega(3)
                       SVelocity(3)=PO%Vz+0.5d0*PO%Az+PO%Vx*Omega(2)-PO%Vy*Omega(1)
                       
                       PO%Vx=TransB(1,1)*SVelocity(1)+TransB(2,1)*SVelocity(2)+TransB(3,1)*SVelocity(3)
                       PO%Vy=TransB(1,2)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(3,2)*SVelocity(3)
                       PO%Vz=TransB(1,3)*SVelocity(1)+TransB(2,3)*SVelocity(2)+TransB(3,3)*SVelocity(3)
                       
                       PO%X=PO%X+PO%Vx
                       PO%Y=PO%Y+PO%Vy
                       PO%Z=PO%Z+PO%Vz
                       
         End subroutine  MoveElectromagneticParticleOne
  
        subroutine IndexInitializationParticleOneIndex(POI,Index)
            Class(ParticleOneIndex), intent(inout) :: POI
            Integer(4),intent(in) :: Index
            POI%Index=Index
         End subroutine IndexInitializationParticleOneIndex
         
         subroutine CopyParticleOneIndex(POD,POC)
            Class(ParticleOneIndex), intent(inout) :: POD
            Class(ParticleOne), intent(in) :: POC
            Select Type (POD)
                 Type is (ParticleOneIndex)
                      Select Type (POC)
                          Type is (ParticleOneIndex)
                          POD=POC
                     ENd Select
            ENd Select
         End subroutine  CopyParticleOneIndex
         
        Subroutine SwapParticleOneIndex(POD,POC)
            Class(ParticleOneIndex), intent(inout) :: POD
            Class(ParticleOne), intent(inout) :: POC
            Type (ParticleOneIndex) :: POT
            Select Type (POD)
                 Type is (ParticleOneIndex)
                      Select Type (POC)
                        Type is (ParticleOneIndex)
                          POT=POC
                          POC=POD
                          POD=POT
                     ENd Select
            ENd Select
         End subroutine  SwapParticleOneIndex
END Module ModuleParticleOne