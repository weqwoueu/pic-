Module DiagnosticsCollisionRate
    Use ModuleGrid
    Use ModuleParticleBundle
    Use ModuleTypeMCC
    Implicit none
    Integer(4),Parameter  :: NReactionMax=100
                    Type ParticleCollisionRate
                        Integer(4) :: Nx=NxMax
                        Integer(4) :: NReaction
                        Real(8) :: CollisionRate(NxMax,NReactionMax)
								Real(8) :: CollisionPowerRate(NxMax,NReactionMax)
                    EndType ParticleCollisionRate
    contains
    
    subroutine DiagParticleCollisionRateOne(GD,CF,PB,MCCB,Mode)
            Implicit none
            Class(*),intent(inout)  :: GD
            Type(ParticleBundle),intent(in) :: PB
				Type(ControlFlow),intent(in) :: CF
            Type(MCCBundle),intent(in) ::  MCCB
            Integer(4),intent(in) ::  Mode
            Type(ParticleCollisionRate),save :: PCR
            Integer(4) :: Shift,i
            Select Type (GD)
                  Type is (Grid1D(*,*))
                     Select Case (Mode)
                        case(0)
                            Shift=1
                            Call WeightingParticleCollisionRate2(PB,CF,PCR,MCCB)
                            do i=1,MCCB%NReaction
                                 Call GD%Update(GD%Nx,PCR%CollisionRate(:,i),Shift)
											Call GD%Update(GD%Nx,PCR%CollisionPowerRate(:,i),Shift)
                            End Do
                            GD%Timer=GD%Timer+1
									 
                         Case(1) 
                                 Call GD%Rescale()
                                 Call GD%Dump(1)
                                 Call GD%Reset()
                         Case(2)
                             Call GD%Dump(0)
                         case default
                         End Select
                    Type is (Grid2D(*,*,*))
                        Select Case (Mode)
                            case(0)
                                Shift=1
                                !Call WeightingParticleCollisionRate2(PB,CF,PCR,MCCB)
                                do i=1,MCCB%NReaction
                                      Call GD%Update(GD%Nx,PCR%CollisionRate(:,i),Shift)
												  Call GD%Update(GD%Nx,PCR%CollisionPowerRate(:,i),Shift)
                                End Do
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
            return !sort()
    end subroutine DiagParticleCollisionRateOne
    
       subroutine WeightingParticleCollisionRate(PB,PCR,MCCB)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(ParticleCollisionRate),intent(inout) :: PCR
                Type(MCCBundle),intent(in) ::  MCCB
                Real(8) :: EnergyInterval,Energy,S1,S2
                Real(8) :: TempProbility(MCCB%NReaction)
                Integer(4) :: i,j,Nc,Np,Index,Upper,Center,Lower,NReaction
                PCR%CollisionRate=0.d0
                NReaction=MCCB%NReaction
                PCR%NReaction=NReaction
                do j=1,PB%Npar
                   Energy=PB%PO(j)%Energy(PB%Mass,PB%VFactor)/JtoeV
                   Associate(Emin=>MCCB%EnergyMin,Eint=>MCCB%EnergyInterval,Nr=>MCCB%NReaction,Ns=>MCCB%NSigma,Probility=>MCCB%Probility,EnergyMax=>MCCB%EnergyMax)
                   TempProbility=0.d0
                   Nc=Int((Log10(Energy)-Emin)/Eint)
                   If (Nc<1) Then
                        TempProbility=Probility(1:Nr,1)
                    ELse If (Nc<Ns) Then
                       S2=(Log10(Energy)-Emin)/Eint-Dble(Nc)
                       S1=1.d0-S2
                       TempProbility(1:Nr)=Probility(1:Nr,Nc)*S1+Probility(1:Nr,Nc+1)*S2
                    Else
                      TempProbility(1:Nr)=Probility(1:Nr,Ns)*Energy/EnergyMax
                   End If
                
                   TempProbility=TempProbility*MCCB%SigmaMax*MCCB%GO%Density
                   Do i=NReaction,2
                         If ((TempProbility(i)-TempProbility(i-1))>=MinReal) Then
                              TempProbility(i)=TempProbility(i)-TempProbility(i-1)
                         End If
                   End do
                   End Associate  
                   Np=Ceiling(PB%PO(j)%X)
                   Do  i=1, NReaction
                          PCR%CollisionRate(Np,i)=PCR%CollisionRate(Np,i)+TempProbility(i)
                   End DO
                end do
                return
		 end subroutine WeightingParticleCollisionRate
		 
		 	 subroutine WeightingParticleCollisionRate2(PB,CF,PCR,MCCB)
            implicit none
			   Type(ParticleBundle),intent(in) :: PB
				Class(ControlFlow), intent(in) :: CF
            Type(ParticleCollisionRate),intent(inout) :: PCR
            Type(MCCBundle),intent(in) ::  MCCB
            Real(8) :: EnergyInterval,Energy,S1,S2,Residue
				Real(8) :: TempWeight,EnergyTemp(3)=0.d0
            Real(8) :: TempProbility(MCCB%NReaction),TempCollisionRate(3)
            Integer(4) :: i,j,k,N,Np,Index,Upper,Center,Lower,NReaction,NCollisionOne
				TempWeight = PB%Weight!*CF%Dx*CF%ElectrodeArea
            PCR%CollisionRate = 0.d0
				PCR%CollisionPowerRate = 0.d0
            NReaction=MCCB%NReaction
            PCR%NReaction=NReaction
				NCollisionOne=Int(PB%NPar*MCCB%CollisionRatio)
            Residue=PB%NPar*MCCB%CollisionRatio-Dble(NCollisionOne)
            Call Random_Number(R)
            If (R<Residue) Then
                   NCollisionOne=NCollisionOne+1
				End If
				EnergyTemp=0.d0
            do j=1,NCollisionOne
					Call Random_Number(R)
               k=(j+1)+(PB%NPar-(j+1))*R
               Energy=PB%PO(k)%Energy(PB%Mass,PB%VFactor)/JtoeV
               Associate(Emin=>MCCB%EnergyMin,Eint=>MCCB%EnergyInterval,Nr=>MCCB%NReaction,Ns=>MCCB%NSigma,Probility=>MCCB%Probility,EnergyMax=>MCCB%EnergyMax)
                TempProbility=0.d0		
                N=Ceiling((Log10(Energy)-Emin)/Eint)
                If (N<1) Then
                    TempProbility=Probility(1:Nr,1)
					 ELse If (N<Ns) Then
                     S1=Dble(N)-(Log10(Energy)-Emin)/Eint
                     S2=1.d0-S1
                     !If(S1<0.d0.or.S1>1.d0.or.S2<0.d0.or.S2>1.d0) Write(*,*) "ERRRRRR S1S2SelectProbility"
                     TempProbility(1:Nr)=Probility(1:Nr,N)*S1+Probility(1:Nr,N+1)*S2
                Else
                     TempProbility(1:Nr)=Probility(1:Nr,Ns)*Energy/EnergyMax
                End If
                DO i=1,Nr
                       If (Energy<MCCB%Reaction(i)%Threshold) Then
                            TempProbility(i)=0.d0
                       End If
                 End Do
                 Call Random_Number(R)
					  Np=Ceiling(PB%PO(k)%X)
					  S1=Dble(Np)-PB%PO(k)%X
                 S2=1.d0-S1
                 Do i=1,Nr
                       If (R<TempProbility(i)) Then		
									 PCR%CollisionRate(Np,i)=PCR%CollisionRate(Np,i)+S1
									 PCR%CollisionRate(Np+1,i)=PCR%CollisionRate(Np+1,i)+S2
									 exit
							  End If
                  End Do	 
                  End Associate  
				end do
				do i=1,MCCB%NReaction
				   PCR%CollisionPowerRate(:,i) =  MCCB%Reaction(i)%Threshold*JtoeV*PCR%CollisionRate(:,i)
				end do
				PCR%CollisionRate = PCR%CollisionRate*TempWeight/CF%Dt
				PCR%CollisionPowerRate = PCR%CollisionPowerRate*TempWeight/CF%Dt
				!if(PB%Mass==ElectronMass) then
				!	PCR%ParticleAdded(1) = sum(PCR%CollisionRate(:,1))*TempWeight
				!	PCR%ParticleAdded(2) = sum(PCR%CollisionRate(:,2))*TempWeight
				!	PCR%ParticleAdded(3) = sum(PCR%CollisionRate(:,3))*TempWeight
				!	PCR%ParticleAdded = PCR%ParticleAdded*TempWeight
				!	PCR%EnergyLoss(1) = 0.d0 * PCR%ParticleAdded(1)*JtoeV
				!	PCR%EnergyLoss(2) = 11.5d0 * PCR%ParticleAdded(2)*JtoeV
				!	PCR%EnergyLoss(3) = 15.8d0 * PCR%ParticleAdded(3)*JtoeV
				!else
				!	PCR%ParticleAdded = 0.d0
				!	PCR%EnergyLoss(1) = 0.5d0*EnergyTemp(1)*TempWeight*JtoeV
				!	PCR%EnergyLoss(2) = 1.0d0*EnergyTemp(2)*TempWeight*JtoeV
				!	PCR%EnergyLoss(3)	= 0.d0
				!end if
            return
        end subroutine WeightingParticleCollisionRate2

 End Module DiagnosticsCollisionRate 