Module DiagParticleBouncary
	use ModuleParticleBoundary
	use ModuleGrid
	use ModuleSEE
	USE ModuleControlFlow
	implicit none
	
	 Type DiagParticleBoundaryOneData
          Real(8) :: ParticleNumber(2),ParticleFlux(2),ParticleCurrent(2),Energy(2),AverageEnergy(2),PowerCurrent(2),PowerFlux(2)
			 Real(8) :: SEEParticleNumber(2),SEEParticleFlux(2),SEEParticleCurrent(2),SeeEnergy(2),SEEAverageEnergy(2),SEEPowerCurrent(2),SEEPowerFlux(2)
			 Real(8) :: TEEParticleNumber(2),TEEParticleFlux(2),TEEParticleCurrent(2),TeeEnergy(2),TEEAverageEnergy(2),TEEPowerCurrent(2),TEEPowerFlux(2)
    EndType DiagParticleBoundaryOneData
	
	contains
	subroutine DiagParticleBoundaryOne(GD,PBDO,CF,Mode)
	    !class(Grid0D(*)),intent(inout) :: GD
	    Class(*),intent(inout)  :: GD
		 Class(ParticleBoundaryOne),intent(in) :: PBDO(:)
		 tYPE(ControlFlow), intent(in) :: CF
		 integer(4),intent(in) :: Mode
		 Type(DiagParticleBoundaryOneData),save :: TempDPBDO(NspacyMax)
		 integer(4) :: i,shift
		  Select Type (GD)
             Type is (Grid1D(*,*))
                 Select Case (Mode)
                    case(0)
                        Shift=1
                        do i=1,CF%Ns+1
				                    call DiagParticleBoundaryCaculate(PBDO(i),CF,TempDPBDO(i))
				                    Call GD%Update(1,TempDPBDO(i)%ParticleFlux(1),Shift)
										  Call GD%Update(1,TempDPBDO(i)%ParticleFlux(2),Shift)
				                    !Call GD%Update(1,TempDPBDO(i)%SEEParticleFlux(1),Shift)
				                    !Call GD%Update(1,TempDPBDO(i)%SEEParticleFlux(2),Shift)
				        !            if(i==1) then
				        !                !Call GD%Update(1,SUM(TempDPBDO(i)%TEEParticleFlux(:)),Shift)
										  !end if
				                      !Power
			                      !Call GD%Update(1,SUM(TempDPBDO(i)%PowerFlux(:)),Shift)
										   Call GD%Update(1,TempDPBDO(i)%PowerFlux(1),Shift)
										 Call GD%Update(1,TempDPBDO(i)%PowerFlux(2),Shift)
				                   !Call GD%Update(1,TempDPBDO(i)%SEEPowerFlux(1),Shift)
				                   !Call GD%Update(1,TempDPBDO(i)%SEEPowerFlux(2),Shift)
				                   !if(i==1) then
				                   !     !Call GD%Update(1,SUM(TempDPBDO(i)%TEEPowerFlux(:)),Shift)
				                   !end if
			               end do
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
                            do i=1,CF%Ns+1
				                    !call DiagParticleBoundaryCaculate(PBDO(i),CF,TempDPBDO(i))
				                    Call GD%Update(1,TempDPBDO(i)%ParticleFlux(1),Shift)
										  Call GD%Update(1,TempDPBDO(i)%ParticleFlux(2),Shift)
				                    !Call GD%Update(1,TempDPBDO%SEEParticleFlux(1),Shift)
				                    !Call GD%Update(1,TempDPBDO%SEEParticleFlux(2),Shift)
				        !            if(i==1) then
				        !                !Call GD%Update(1,SUM(TempDPBDO%TEEParticleFlux(:)),Shift)
										  !end if
				                      !Power
			                      Call GD%Update(1,TempDPBDO(i)%PowerFlux(1),Shift)
										 Call GD%Update(1,TempDPBDO(i)%PowerFlux(2),Shift)
				                   !Call GD%Update(1,TempDPBDO(i)%SEEPowerFlux(1),Shift)
				                   !Call GD%Update(1,TempDPBDO(i)%SEEPowerFlux(2),Shift)
				                   !if(i==1) then
				                   !     !Call GD%Update(1,SUM(TempDPBDO(i)%TEEPowerFlux(:)),Shift)
				                   !end if
			               end do
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
		 
		 !Select Case (Mode)
		 !case(0)
			! Shift=1
			! GD%Timer=GD%Timer+1
		 !case(1)
			!   Call GD%Rescale
   !         Call GD%Dump(1)
   !          GD%Value=0.d0
   !          GD%Timer=0
		 !case(2)
			! Call GD%Dump(0)
		 !end  Select
	  return
	end subroutine DiagParticleBoundaryOne
	
	subroutine DiagParticleBoundaryCaculate(PBDO,CF,DPBDO)
	   implicit none
		Type(DiagParticleBoundaryOneData),INTENT(INOUT) :: DPBDO
		Class(ParticleBoundaryOne),intent(in) :: PBDO
		Class(ControlFlow), intent(in) :: CF
		integer(4) :: i,j
		
		DPBDO%ParticleNumber   = 0.d0
		DPBDO%SeeParticleNumber= 0.d0
		DPBDO%TeeParticleNumber= 0.d0
		DPBDO%ParticleCurrent  = 0.d0
		DPBDO%SeeParticleCurrent = 0.d0
		DPBDO%TeeParticleCurrent = 0.d0
		DPBDO%ParticleFlux    = 0.d0
		DPBDO%SeeParticleFlux = 0.d0
		DPBDO%TeeParticleFlux = 0.d0
		DPBDO%Energy    = 0.d0
		DPBDO%SEEEnergy = 0.d0
		DPBDO%TEEEnergy = 0.d0
		DPBDO%AverageEnergy   = 0.d0
		DPBDO%SEEAverageEnergy= 0.d0
		DPBDO%TEEAverageEnergy= 0.d0
		DPBDO%PowerCurrent    = 0.d0
		DPBDO%SeePowerCurrent = 0.d0
		DPBDO%TeePowerCurrent = 0.d0
		DPBDO%PowerFlux       = 0.d0
		DPBDO%SeePowerFlux    = 0.d0
		DPBDO%TeePowerFlux    = 0.d0
		!Lower electrode
		do i=1,PBDO%CountMinOne
			DPBDO%Energy(1)= DPBDO%Energy(1)+PBDO%PBLower%PO(i)%Energy(PBDO%PBLower%Mass,PBDO%PBLower%VFactor)
		end do
		do i=1,PBDO%SEECountMinOne
			j=PBDO%CountMinOne+i
			DPBDO%SeeEnergy(1)= DPBDO%SeeEnergy(1)+PBDO%PBLower%PO(j)%Energy(electronmass,PBDO%PBLower%VFactor)
		end do
		do i=1,PBDO%TEECountMinOne
			j=PBDO%CountMinOne+PBDO%SEECountMinOne+i
			DPBDO%TeeEnergy(1)= DPBDO%TeeEnergy(1)+PBDO%PBLower%PO(j)%Energy(electronmass,PBDO%PBLower%VFactor)
		end do
		!Upper electrode
		do i=1,PBDO%CountMaxOne
			DPBDO%Energy(2) = DPBDO%Energy(2) + PBDO%PBUpper%PO(i)%Energy(PBDO%PBUpper%Mass,PBDO%PBUpper%VFactor)
		end do
		do i=1,PBDO%SEECountMaxOne
			j=PBDO%CountMaxOne+i
			DPBDO%SEEEnergy(2)= DPBDO%SEEEnergy(2)+PBDO%PBUpper%PO(j)%Energy(electronmass,PBDO%PBUpper%VFactor)
		end do
		do i=1,PBDO%TEECountMaxOne
			j=PBDO%CountMaxOne+PBDO%SEECountMaxOne+i
			DPBDO%TEEEnergy(2)= DPBDO%TEEEnergy(2)+PBDO%PBUpper%PO(j)%Energy(electronmass,PBDO%PBUpper%VFactor)
		end do
		
		DPBDO%ParticleNumber(1)	= DBLE(PBDO%CountMinOne) * PBDO%PBLower%Weight*CF%Dx
		DPBDO%ParticleNumber(2)	= DBLE(PBDO%CountMAXOne) * PBDO%PBUpper%Weight*CF%Dx
		DPBDO%ParticleCurrent	= DPBDO%ParticleNumber / CF%Dt 
		DPBDO%ParticleFlux	= DPBDO%ParticleCurrent * CF%ElectrodeArea
		DPBDO%Energy(1)	   = DPBDO%Energy(1) * PBDO%PBLower%Weight*CF%Dx
		DPBDO%Energy(2)	   = DPBDO%Energy(2) * PBDO%PBUpper%Weight*CF%Dx
		DPBDO%PowerCurrent   = DPBDO%Energy / CF%Dt 
		DPBDO%PowerFlux      = DPBDO%PowerCurrent * CF%ElectrodeArea
		DPBDO%AverageEnergy	= DPBDO%Energy/ (JtoeV * DPBDO%ParticleNumber)
		
		DPBDO%SEEParticleNumber(1) = DBLE(PBDO%SEECountMinOne) * PBDO%PBLower%Weight*CF%Dx
		DPBDO%SEEParticleNumber(2)	= DBLE(PBDO%SEECountMaxOne) * PBDO%PBUpper%Weight*CF%Dx
		DPBDO%SEEParticleCurrent = DPBDO%SEEParticleNumber / CF%Dt 
		DPBDO%SEEParticleFlux	 = DPBDO%SEEParticleCurrent * CF%ElectrodeArea
		DPBDO%SEEEnergy		   = DPBDO%SEEEnergy * PBDO%PBLower%Weight*CF%Dx
		DPBDO%SEEPowerCurrent   = DPBDO%SEEEnergy / CF%Dt 
		DPBDO%SEEPowerFlux      = DPBDO%SEEPowerCurrent * CF%ElectrodeArea
		DPBDO%SeeAverageEnergy  = DPBDO%seeEnergy/ (JtoeV*DPBDO%SEEParticleNumber)
		
		DPBDO%TEEParticleNumber(1) = DBLE(PBDO%TEECountMinOne) * PBDO%PBLower%Weight*CF%Dx
		DPBDO%TEEParticleNumber(2)	= DBLE(PBDO%TEECountMaxOne) * PBDO%PBUpper%Weight*CF%Dx
		DPBDO%TEEParticleCurrent = DPBDO%TEEParticleNumber / CF%Dt 
		DPBDO%TEEParticleFlux	 = DPBDO%TEEParticleCurrent * CF%ElectrodeArea
		DPBDO%TEEEnergy		   = DPBDO%TEEEnergy * PBDO%PBLower%Weight*CF%Dx
		DPBDO%TEEPowerCurrent   = DPBDO%TEEEnergy / CF%Dt 
		DPBDO%TEEPowerFlux      = DPBDO%TEEPowerCurrent * CF%ElectrodeArea
		DPBDO%TeeAverageEnergy  = DPBDO%TeeEnergy/ (JtoeV*DPBDO%SEEParticleNumber)
	  return
	end subroutine DiagParticleBoundaryCaculate
	
end module DiagParticleBouncary