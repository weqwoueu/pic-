Module DiagnosticsEEPF
    use ModuleParticleBoundary
    Use ModuleGrid
    Use ModuleParticleBundle
    Implicit none

   Integer(4),Parameter,Private :: NeMax=1000
   Type ParticleEDF
         Integer(4) :: Ne=NeMax
         Real(8) :: EnergyInterval,Mass
         Real(8) :: EDF(NeMax),EDFNormalized(NeMax)
   EndType ParticleEDF
   
   Type ParticleEPF
         Integer(4) :: Ne=NeMax, NMax, NMaxAver
         Real(8) :: EnergyInterval,Mass
         Real(8) :: EPF(NeMax),EPFNormalized(NeMax)
         !Real(8) :: AverEPF(NeMax)=0,AverEPFNormalized(NeMax)=0
   EndType ParticleEPF
    contains
        subroutine DiagParticleEDFOne(GD,PB,PBDO,Mode)
            Implicit none
            Class(*),intent(inout)  :: GD
            Type(ParticleBundle),intent(in) :: PB
			   Type(ParticleBoundaryOne),intent(in) :: PBDO
            Integer(4),intent(in) ::  Mode
            Type(ParticleEDF),SAVE :: PEDF(3,0:NspacyMax)
            Integer(4) :: Shift,i,ParticleType=0

            If (PB%SO%SpecyIndex==0) Then
					 ParticleType=0
					 PEDF%EnergyInterval = 1.d-1
            Else
                PEDF%EnergyInterval = 1.d0
					 ParticleType=PB%SO%SpecyIndex
            End If
            PEDF%Mass=PB%Mass
            Select Type (GD)
                  Type is (Grid1D(*,*))
                     Select Case (Mode)
                        case(0)
                            PEDF%Ne=GD%Nx
                            Shift=1
                            Call WeightingParticleEDF(PB,PEDF(1,ParticleType))
                            Call GD%Update(GD%Nx,PEDF(1,ParticleType)%EDF,Shift)
                            Call GD%Update(GD%Nx,PEDF(1,ParticleType)%EDFNormalized,Shift)
									 !AEDF
									 Call WeightingParticleEDF2(PBDO%PBLower,PEDF(2,ParticleType),1,PBDO%CountMinOne)
                            Call GD%Update(GD%Nx,PEDF(2,ParticleType)%EDF,Shift)
                            Call GD%Update(GD%Nx,PEDF(2,ParticleType)%EDFNormalized,Shift)
									 !Shift=Shift-2
									 Call WeightingParticleEDF2(PBDO%PBUpper,PEDF(3,ParticleType),1,PBDO%CountMaxOne)
                            Call GD%Update(GD%Nx,PEDF(3,ParticleType)%EDF,Shift)
                            Call GD%Update(GD%Nx,PEDF(3,ParticleType)%EDFNormalized,Shift)
									 !see EEDF
									 !Call WeightingParticleEDF2(PBDO%PBLower,PEDF(4,ParticleType),PBDO%CountMinOne+1,PBDO%CountMinOne+PBDO%SEECountMinOne)
          !                  Call GD%Update(GD%Nx,PEDF(4,ParticleType)%EDF,Shift)
          !                  Call GD%Update(GD%Nx,PEDF(4,ParticleType)%EDFNormalized,Shift)
									 !Shift=Shift-2
									 !Call WeightingParticleEDF2(PBDO%PBUpper,PEDF(5,ParticleType),PBDO%CountMaxOne+1,PBDO%CountMaxOne+PBDO%SEECountMaxOne)
          !                  Call GD%Update(GD%Nx,PEDF(5,ParticleType)%EDF,Shift)
          !                  Call GD%Update(GD%Nx,PEDF(5,ParticleType)%EDFNormalized,Shift)
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
             !                  Call WeightingParticleEDF(PB,PEDF(1,ParticleType))
                               Call GD%Update(GD%Nx,PEDF(1,ParticleType)%EDF,Shift)
                               Call GD%Update(GD%Nx,PEDF(1,ParticleType)%EDFNormalized,Shift)
									    !AEDF
									    Call WeightingParticleEDF2(PBDO%PBLower,PEDF(2,ParticleType),1,PBDO%CountMinOne)
                               Call GD%Update(GD%Nx,PEDF(2,ParticleType)%EDF,Shift)
                               Call GD%Update(GD%Nx,PEDF(2,ParticleType)%EDFNormalized,Shift)
									    !Shift=Shift-2
									    Call WeightingParticleEDF2(PBDO%PBUpper,PEDF(3,ParticleType),1,PBDO%CountMaxOne)
                               Call GD%Update(GD%Nx,PEDF(3,ParticleType)%EDF,Shift)
                               Call GD%Update(GD%Nx,PEDF(3,ParticleType)%EDFNormalized,Shift)
									    !see EEDF
									    !Call WeightingParticleEDF2(PBDO%PBLower,PEDF(4,ParticleType),PBDO%CountMinOne+1,PBDO%CountMinOne+PBDO%SEECountMinOne)
          !                  Call GD%Update(GD%Nx,PEDF(4,ParticleType)%EDF,Shift)
          !                  Call GD%Update(GD%Nx,PEDF(4,ParticleType)%EDFNormalized,Shift)
									 !Shift=Shift-2
									 !Call WeightingParticleEDF2(PBDO%PBUpper,PEDF(5,ParticleType),PBDO%CountMaxOne+1,PBDO%CountMaxOne+PBDO%SEECountMaxOne)
          !                     Call GD%Update(GD%Nx,PEDF(5,ParticleType)%EDF,Shift)
          !                     Call GD%Update(GD%Nx,PEDF(5,ParticleType)%EDFNormalized,Shift)
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
            return
        end subroutine DiagParticleEDFOne
        
        Subroutine DiagParticleEPFOne(GD,PB,Mode)
            Implicit none
            Class(*),intent(inout)  :: GD
            Type(ParticleBundle),intent(in) :: PB
            Integer(4),intent(in) ::  Mode
            Type(ParticleEPF),SAVE :: PEPF(1,0:NspacyMax)
            Integer(4) :: Shift,i,ParticleType=0
            
            If (PB%SO%SpecyIndex==0) Then
				ParticleType=0
				PEPF%EnergyInterval = 1.d-1
            Else
                ParticleType=PB%SO%SpecyIndex
                PEPF%EnergyInterval = 1.d0
            End If
            Select Type (GD)
            Type is (Grid1D(*,*))
                Select Case (Mode)
                Case(0)
                    PEPF%Ne=GD%Nx
                    Shift=1
                    Call WeightingParticleEPF(PB,PEPF(1,ParticleType))
                    Call GD%Update(GD%Nx,PEPF(1,ParticleType)%EPF,Shift)
                    Call GD%Update(GD%Nx,PEPF(1,ParticleType)%EPFNormalized,Shift)
                    GD%Timer=GD%Timer+1
                Case(1) 
                    Call GD%Rescale
                    Call GD%Dump(1)
                    GD%Value=0.d0
                    GD%Timer=0
                Case(2)
                    !Call GD%Dump(0)
                case default
                End Select
            Type is (Grid2D(*,*,*))
                Select Case (Mode)
                Case(0)
                    GD%Timer=GD%Timer+1
                Case(1) 
                    Call GD%Rescale
                    Call GD%Dump(1)
                    GD%Value=0.d0
                    GD%Timer=0
                Case(2)
                    Call GD%Dump(0)
                Case default
                End Select     
            End select
            return
        End subroutine DiagParticleEPFOne
        
             
        subroutine WeightingParticleEDF(PB,PEDF)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(ParticleEDF),intent(inout) :: PEDF
                Real(8) :: Energy,Frac
                Integer(4) :: i,N
                
                PEDF%EDF=0.d0
                PEDF%EDFNormalized=0.d0
                
                Frac=1.d0/DBLE(PB%Npar)
                do i=1,PB%Npar
                   Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)/JtoeV
                   N=Ceiling( Energy/PEDF%EnergyInterval)
                   If (N>=1.and.N<PEDF%Ne) Then
                           PEDF%EDF(N) = PEDF%EDF(N)+PB%Weight*PB%Dx
                           PEDF%EDFNormalized(N)=PEDF%EDFNormalized(N)+Frac
                   End IF
                end do
                return
		  end subroutine WeightingParticleEDF     
		  
		  subroutine WeightingParticleEDF2(PB,PEDF,NumLow,NumMax)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(ParticleEDF),intent(inout) :: PEDF
					 integer,intent(in) :: NumLow
					 integer,intent(in) :: NumMax
                Real(8) :: Energy,Frac
                Integer(4) :: i,N
                
                PEDF%EDF=0.d0
                PEDF%EDFNormalized=0.d0
                Frac=1.d0/DBLE(PB%Npar)
                do i=NumLow,NumMax
                   Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)/JtoeV
                   N=Ceiling( Energy/PEDF%EnergyInterval)
                   If (N>=1.and.N<PEDF%Ne) Then
                           PEDF%EDF(N) = PEDF%EDF(N)+PB%Weight*PB%Dx
                           PEDF%EDFNormalized(N)=PEDF%EDFNormalized(N)+Frac
                   End IF
                end do
                return
          end subroutine WeightingParticleEDF2     
          
          
          subroutine WeightingParticleEPF(PB,PEPF)
                implicit none
                Type(ParticleBundle),intent(in) :: PB
                Type(ParticleEPF),intent(inout) :: PEPF
                Real(8) :: Energy,Frac,NParReal
                Integer(4) :: i,N
                
                PEPF%EPF=0.d0
                PEPF%EPFNormalized=0.d0
                
                If (PB%UnequalWeightFlag) Then
                    NparReal = 0.
                    do i=1,PB%Npar
                        NParReal = NparReal + PB%PO(i)%Wq
                    end do
                    do i=1,PB%Npar
                        Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)/JtoeV
                        N=Ceiling( Energy/PEPF%EnergyInterval)
                        If (N>=1.and.N<PEPF%Ne) Then
                            PEPF%EPF(N) = PEPF%EPF(N)+PB%PO(i)%Wq
                            Frac = PB%PO(i)%Wq/NParReal/Sqrt(Energy)
                            PEPF%EPFNormalized(N)=PEPF%EPFNormalized(N)+Frac
                        End IF
                    end do
                Else
                    Do i=1,PB%Npar
                        Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)/JtoeV
                        N=Ceiling( Energy/PEPF%EnergyInterval)
                        If (N>=1.and.N<PEPF%Ne) Then
                            PEPF%EPF(N) = PEPF%EPF(N)+PB%Weight
                            Frac=1.d0/DBLE(PB%Npar)/Sqrt(Energy)
                            PEPF%EPFNormalized(N)=PEPF%EPFNormalized(N)+Frac
                        End IF
                    end do
                Endif
                !PEPF%AverEPF = PEPF%AverEPF + PEPF%EPF
                !PEPF%AverEPFNormalized = PEPF%AverEPFNormalized + PEPF%EPFNormalized
                return
          end subroutine WeightingParticleEPF  
		  
        !Subroutine AveragingEPF(PEPF,Timer)
        !    implicit none
        !    Type(ParticleEPF),intent(inout) :: PEPF
        !    Integer(4),intent(in) :: Timer
        !    PEPF%AverEPF = PEPF%AverEPF/Timer
        !    PEPF%AverEPFNormalized = PEPF%AverEPFNormalized/Timer
        !End Subroutine AveragingEPF
        !  
        !Subroutine ZeroOutEPF(PEPF)
        !    implicit none
        !    Type(ParticleEPF),intent(inout) :: PEPF
        !    Integer(4) :: i
        !    Do i = 1, PEPF%Ne
        !        PEPF%EPF(i) = 0.
        !        PEPF%AverEPF(i) = 0
        !        PEPF%AverEPFNormalized(i) = 0
        !    End Do
        !End Subroutine ZeroOutEPF
End Module DiagnosticsEEPF
	
	
	!Module DiagnosticsEEPF
!    Use ModuleGrid
!    Use ModuleParticleBundle
!    Implicit none
!
!                   Integer(4),Parameter,Private :: NeMax=1000
!                    Type ParticleEDF
!                        Integer(4) :: Ne=NeMax
!                        Real(8) :: EnergyInterval,Mass
!                        Real(8) :: EDF(NeMax),EDFNormalized(NeMax)
!                    EndType ParticleEDF
!    contains
!             subroutine DiagParticleEDFOne(GD,PB,Mode)
!            Implicit none
!            Class(*),intent(inout)  :: GD
!            Type(ParticleBundle),intent(in) :: PB
!            Integer(4),intent(in) ::  Mode
!            Type(ParticleEDF) :: PEDF
!            Integer(4) :: Shift
!
!            If (PB%SO%SpecyIndex==0) Then
!                PEDF%EnergyInterval=0.1d0
!            Else
!                PEDF%EnergyInterval=1.d0
!            End If
!            PEDF%Mass=PB%Mass
!
!            Select Type (GD)
!                  Type is (Grid1D(*,*))
!                     Select Case (Mode)
!                        Case(-1)
!                            !Call GridInitialization(GD,PB%Period,PEDF%EnergyInterval,PB%Dt)
!                        case(0)
!                            
!                            PEDF%Ne=GD%Nx
!                            Shift=1
!                            Call WeightingParticleEDF(PB,PEDF)
!                            Call GD%Update(GD%Nx,PEDF%EDF,Shift)
!                            Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
!                            GD%Timer=GD%Timer+1
!                         Case(1) 
!                                 Call GD%Rescale
!                                 Call GD%Dump(1)
!                                 GD%Value=0.d0
!                                 GD%Timer=0
!                         Case(2)
!                             Call GD%Dump(0)
!                         case default
!                         End Select
!                    Type is (Grid2D(*,*,*))
!                        Select Case (Mode)
!                            Case(-1)
!                                !Call GridInitialization(GD,PB%Period,PEDF%EnergyInterval,PB%Dt)
!                            case(0)
!                                Shift=1
!                                Call WeightingParticleEDF(PB,PEDF)
!                                Call GD%Update(GD%Nx,PEDF%EDF,Shift)
!                                Call GD%Update(GD%Nx,PEDF%EDFNormalized,Shift)
!                                GD%Timer=GD%Timer+1
!                            Case(1) 
!                                     Call GD%Rescale
!                                     Call GD%Dump(1)
!                                     GD%Value=0.d0
!                                     GD%Timer=0
!                             Case(2)
!                                 Call GD%Dump(0)
!                             case default
!                         End Select     
!                    End select
!            return
!             end subroutine DiagParticleEDFOne
!             
!        subroutine WeightingParticleEDF(PB,PEDF)
!                implicit none
!                Type(ParticleBundle),intent(in) :: PB
!                Type(ParticleEDF),intent(inout) :: PEDF
!                Real(8) :: Energy,Frac
!                Integer(4) :: i,N
!                
!                PEDF%EDF=0.d0
!                PEDF%EDFNormalized=0.d0
!                
!                Frac=1.d0/DBLE(PB%Npar)
!                do i=1,PB%Npar
!                   Energy=PB%PO(i)%Energy(PB%Mass,PB%VFactor)/JtoeV
!                   N=Ceiling( Energy/PEDF%EnergyInterval)
!                   If (N>=1.and.N<PEDF%Ne) Then
!                           !PEDF%EDF(N)=PEDF%EDF(N)+PB%Weight/PB%XMax
!                           PEDF%EDFNormalized(N)=PEDF%EDFNormalized(N)+Frac
!                   End IF
!                end do
!                return
!    end subroutine WeightingParticleEDF     
!End Module DiagnosticsEEPF