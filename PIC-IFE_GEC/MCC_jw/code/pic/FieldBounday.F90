Module ModuleFieldBoundary
  Use ModuleControlFlow
  Use ModuleField
  Use ModuleParticleBoundary
  !Use ModuleParticleBoundary
  Implicit none
     Type FieldBoundary
                        !Integer(4) :: XStart=0,XEnd=1
                         Integer(4) :: FieldBoundaryModel=11
                         Integer(4) :: Timer=0,Period
								 Integer(4) :: BiasMode = 0  !0->No bias Voltage  1->Power electrode bias 2->asymatic area bias
                         Real(8) :: Dt=Inputdt,AreaRatio
                         Real(8) :: Frequency(2)=(/13.56d6,13.56d6/),Voltage(2)=(/300.d0,300.d0/)
                         Real(8) :: V1=0.d0,V2=0.d0
                         Real(8) :: Vdc=0.d0!-20.d0  Self bios voltage

                         Contains
                         procedure :: Init=>InitilalizationFieldBounday
                         procedure :: Updater=>UpdaterFieldBounday
                          
     EndType FieldBoundary

    contains
     Subroutine InitilalizationFieldBounday(FB,CF)
               Implicit none
               Class(FieldBoundary),intent(inout) :: FB
               Type(ControlFlow),intent(inout) :: CF
               Real(8) :: Dt
               Character(len=99) :: Filename
               Logical :: alive
               integer(4) :: i

               Write(filename,*) "10000","FieldBoundary",".dat"
               Inquire(file=filename,exist=alive)
               If(alive) then
                   Open (10,file=filename)
                   Read(10,*) FB%Timer,FB%V1,FB%V2,FB%Vdc
                   Close(10)
               Else
                   Write(*,*) 'Can not find the file for ', filename,' the FieldBoundary will be set to Zero.' 
                   FB%Timer=0
                   FB%Vdc=0.d0
                   FB%V1=0.d0
                   FB%V2=0.d0
               End If
               CF%Period=Int(1.d0/FB%Frequency(1)/CF%dt)
               FB%Dt=CF%dt
               FB%Period=CF%Period
               FB%AreaRatio=CF%ElectrodeArea/CF%ElectrodeArea2
               return
     End Subroutine InitilalizationFieldBounday
     
     Subroutine UpdaterFieldBounday(FB,PB,FG,PBDO)
            implicit none
            Class(FieldBoundary),intent(inout) :: FB
            class(Field),intent(in) :: FG
            CLASS(ParticleBundle),intent(in) ::  PB(:)
            Class(ParticleBoundaryOne),intent(inOUT) :: PBDO(:)
				
					if (FB%Timer==FB%Period) then 
						FB%Timer=0
					end if
					
				if (FB%BiasMode/=0) then
				   call UpdateBiasCorrection(FB,PBDO,FB%BiasMode)
				end if
            Select case (FB%FieldBoundaryModel)
				
                      case (11)
                               FB%V1=FB%Voltage(1)*DCOs(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))+FB%Vdc
                               FB%V2=0.d0
                               !Write (*,*) FB%V1,FB%V2
                      case (21)
                               FB%V1=FB%Voltage(1)*DCos(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer)) + FB%Voltage(2)*DCos(2*PI*FB%Frequency(2)*FB%dt*Dble(FB%Timer))+FB%Vdc
                               FB%V2=0.d0
                      case(22)
                               FB%V1=FB%Voltage(2)*DCos(2*PI*FB%Frequency(2)*FB%dt*Dble(FB%Timer))
                               FB%V2=FB%Voltage(1)*DCos(2*PI*FB%Frequency(1)*FB%dt*Dble(FB%Timer))  
                      Case(31)

              end  Select
              FB%Timer=FB%Timer+1
            return
     end subroutine  UpdaterFieldBounday

     
     Subroutine FieldBoundayFinalization(FB)
               Implicit none
               Type(FieldBoundary),intent(inout) :: FB
               Character(len=99) :: Filename
               Write(filename,*) "10000","FieldBoundary",".dat"
               Open (10,file=filename)
               Write(10,*) FB%Timer,FB%V1,FB%V2,FB%Vdc
               Close(10)
               Write(*,*) "Saving ",filename," Please wait..."
               return
     End Subroutine FieldBoundayFinalization
     !
     Subroutine UpdateBiasCorrection(FB,PBDO,mode)
               Implicit none
               Type(FieldBoundary),intent(inout) :: FB
               Type(ParticleBoundaryOne),intent(inOUT) :: PBDO(:)
					 Integer(4), intent(in) :: mode
               Integer(4) :: i,NS
               Real(8) :: NetCharge,NetchargeGround
					Real(8) :: ChargeFluxRatio
					Character(len=99) :: Filename
					Ns=SIZE(PBDO(:))
					 do i=1,Ns
						 if(PBDO(i)%PBLower%SO%Charge>0.d0) then
                        PBDO(i)%CountMin = PBDO(i)%CountMin + (PBDO(i)%SEECountMinOne+PBDO(i)%FEECountMinOne+PBDO(i)%TEECountMinOne)
						 else
							  PBDO(i)%CountMin = PBDO(i)%CountMin - (PBDO(i)%SEECountMinOne+PBDO(i)%FEECountMinOne+PBDO(i)%TEECountMinOne)
						 end if
						 if(PBDO(i)%PBUpper%SO%Charge>0.d0) then
                        PBDO(i)%CountMax = PBDO(i)%CountMax + (PBDO(i)%SEECountMaxOne+PBDO(i)%FEECountMaxOne+PBDO(i)%TEECountMaxOne)
						 else
							  PBDO(i)%CountMax = PBDO(i)%CountMax - (PBDO(i)%SEECountMaxOne+PBDO(i)%FEECountMaxOne+PBDO(i)%TEECountMaxOne)
						 end if
					 End Do
                If(Mod(FB%Timer,FB%Period)==0) Then
                   NetCharge=0.d0
						 NetChargeGround=0.d0
                   do i=1,Ns
                        NetCharge = NetCharge + PBDO(i)%CountMin*PBDO(i)%PBLower%SO%Charge*PBDO(i)%PBLower%Weight
								NetChargeGround = NetCharge + PBDO(i)%CountMax*PBDO(i)%PBUpper%SO%Charge*PBDO(i)%PBUpper%Weight
                        PBDO(i)%CountMin=0
                        PBDO(i)%CountMax=0
						 End Do
						 
						  Select case(MODE)
						  case(0)
						  case(1)
							   If (NetCharge<0.d0) then
                            FB%Vdc=FB%Vdc-1.d0
                        Else if (NetCharge>0.d0) Then
                            FB%Vdc=FB%Vdc+1.d0
					      	 ENd IF
						  case(2)
							  If ((-FB%AreaRatio*NetCharge)<NetChargeGround) then
                          FB%Vdc=FB%Vdc+1.d0
                       Else if ((-FB%AreaRatio*NetCharge)>NetChargeGround) Then
                          FB%Vdc=FB%Vdc-1.d0
							  ENd IF							  
						  case(3)
							  
						  end select

                  Write(filename,*) "1SelfBiasVoltage",".dat"
                    Open (10,file=filename,position='APPEND')
						      Write(10,FMt="(*(es21.14,1x))")  dble(FB%Timer)*FB%Dt,FB%V1,FB%V2,FB%Vdc
                       !Write(10,*) FB%Timer,FB%V1,FB%V2,FB%Vdc
                    Close(10)
                  End IF
                  return
     End Subroutine UpdateBiasCorrection
     
end  Module ModuleFieldBoundary

              
              