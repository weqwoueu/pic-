module ModuleSEE
  !use ModuleParticleBoundary
  Use Constants
  use ModuleSpecyOne
  Use ModuleParticleBundle
   implicit none	
	!Horv¨¢th B, Daksha M, Korolov I, Derzsi A and Schulze J 2017 The role of electron induced secondary electron emission from SiO 2 surfaces in capacitively coupled radio frequency plasmas operated at low pressures Plasma Sources Science and Technology 26 124001
   !Type(Gas),parameter :: SEEGas=(Gas(1,0,9.1095d-31,3.0*11605.d0,0.d0,0.d0)) 

	 Type SecondaryElectron
		 Integer(4) :: Number = 0
		 Integer(4) :: SEETpye = 0_4       !0->NoSEE 1->real see, 2->elastic reflection 3->inelastic reflection 4->ion
		  Real(8) ::  EnergyRange(4) = (/20.0d0, 0.d0, 0.0d0, 5.0d0/)
		             !1 real see Energy=0-20eV 
		             !2 elastic reflection Energy=the energy bombarding the electrode
		             !3 inelastic reflection Energy=0~E(the energy bombarding the electrode)
		             !4 ion induced SEE  0-5eV 
		  Real(8) :: ks = 1  !%smoothness factor 0~2(ks=0 rough surface, ks=1(dull surface), ks=2(apolished surface)
        Real(8) :: deltamax0 = 2.5
        Real(8) :: E0 = 15.d0    !eV the threshold energy for the emission of electron induced SEs
        Real(8) :: Emax0 = 400  !the energy at the maximum yield
        Real(8) :: Gamae = 0.03 !elastic reflect control coefficient
        Real(8) :: Gamai = 0.07    !elastic reflect coefficient
        Real(8) :: Ee0 = 0.d0      !eV the threshold energy for the elastic reflect
        Real(8) :: Eemax = 5.0     !5-10eV the electron energy for maximun elastic reflect coefficient
                                   !Barral S, Makowski K and Peradzynsky Z 2003 Phys. Plasmas 10 4137
        Real(8) :: Yitaemax = 0.5  !when Eemax lies within 5-10eV
        Real(8) :: Deltae =5      !a parameter which controls the decay of Yita_e
		  Real(8) :: IonSEEcoefficient = 0.38d0     !a parameter which controls the decay of Yita_e
	end type SecondaryElectron
	CONTAINS
	
	!Subroutine  ElectronIndcedSEE(PB,PBDO,electrode)
 !              Implicit none
 !              class(ParticleBundle),intent(inout) :: PB
 !              class(ParticleBoundaryOne),intent(inout) :: PBDO
	!				Integer(4),intent(in) :: electrode
 !              type(ParticleOne) :: ParticleTemp
	!				Type(SecondaryElectron) :: SEE
 !              Integer(4) :: i,j,ntemp
 !              Real(8) :: EmissionEnergy
	!				 !the Lower electrode
	!				if(electrode == 0) then
	!				 do i = 1,PBDO%CountMinOne
	!					SEE%Number = 0.d0
	!					call SEEEType(PBDO%PBLower%PO(i),PBDO%PBLower%VFactor,SEE)
	!					PBDO%SEECountMinOne = PBDO%SEECountMinOne+SEE%Number
	!					j=0
	!					DO while (j<SEE%Number)
	!						Select case(SEE%SEETpye)
	!						case(0)
	!						case(1)!	 real secondary electrons(0-20eV)
	!								Call Random_NUMBER(R)
	!								 EmissionEnergy = R*20.0d0*JtoeV
	!								!call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO,electrode)
	!								 call SEEAddToBundle2(SEE%EnergyRange(1),ParticleTemp,PB,PBDO,electrode)
	!						case(2)  !elastic reflected
	!								 EmissionEnergy = SEE%EnergyRange(2)
	!								 call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO,electrode)
	!						case(3)!inelastic reflected electron (0~electron energy)
	!								 Call Random_NUMBER(R)
	!								 EmissionEnergy = R*SEE%EnergyRange(3)
	!								 call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO,electrode)
	!								 !call SEEAddToBundle2(SEE%EnergyRange(3),ParticleTemp,PB,PBDO,electrode)			 
	!						End Select 
	!						j=j+1
	!					  end do
	!				  end do
	!				end if
	!				
	!				if (electrode == 1) then !the Upper electrode
	!				do i = 1,PBDO%CountMaxOne
	!					SEE%Number = 0.d0
	!					call SEEEType(PBDO%PBUpper%PO(i),PBDO%PBUpper%VFactor,SEE)
	!					PBDO%SEECountMaxOne = PBDO%SEECountMaxOne+SEE%Number
	!					j=0
	!					DO while (j<SEE%Number)
	!						Select case(SEE%SEETpye)
	!						   case(0)!No see
	!						   case(1)!	 real secondary electrons(0-20eV)
	!						     Call Random_NUMBER(R)
	!								EmissionEnergy = R*20.0d0*JtoeV
	!								!call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO)
	!								call SEEAddToBundle2(SEE%EnergyRange(1),ParticleTemp,PB,PBDO,electrode)  !Maxwell distribution
	!						   case(2)  !elastic reflected
	!								 Call Random_NUMBER(R)
	!								 EmissionEnergy = SEE%EnergyRange(2)
	!								 call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO,electrode)
	!						   case(3)!inelastic reflected electron (0~electron energy)
	!								  Call Random_NUMBER(R)
	!								 EmissionEnergy = R*SEE%EnergyRange(3)
	!								 call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO,electrode)   !Random distribution
	!								 !call SEEAddToBundle2(SEE%EnergyRange(3),ParticleTemp,PB,PBDO,electrode)  !Maxwell distribution
	!						   End Select 
	!							j=j+1
	!						end do
	!				  end do
	!				end if
 !           return
	!end subroutine ElectronIndcedSEE
	!
	!Subroutine  IonIndcedSEE(PB,PBDO,electrode)
 !              Implicit none
 !              class(ParticleBundle),intent(inout) :: PB
 !              class(ParticleBoundaryOne),intent(inout) :: PBDO
	!				Integer(4),intent(in) :: electrode
 !              type(ParticleOne) :: ParticleTemp
	!				Type(SecondaryElectron) :: SEE
 !              Integer(4) :: i,j,ntemp
 !              Real(8) :: Residue,EmissionEnergy
	!				if (electrode == 0)  then !the Lower electrode
	!				 do i = 1,PBDO%CountMinOne
	!					SEE%Number = 0.d0
	!					SEE%SEETpye=4
	!					Call Random_NUMBER(R)
	!					if (R<SEE%IonSEEcoefficient) then
	!						SEE%Number = 1
	!						Call Random_NUMBER(R)
	!						EmissionEnergy = R*SEE%EnergyRange(4)*JtoeV
	!						!call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO,electrode)
	!						call SEEAddToBundle2(SEE%EnergyRange(4),ParticleTemp,PB,PBDO,electrode)
	!					 end if
	!					 PBDO%SEECountMinOne = PBDO%SEECountMinOne+SEE%Number
	!				  end do
	!				end if
	!				if (electrode == 1) then !the Upper electrode
	!				do i = 1,PBDO%CountMaxOne
	!					SEE%Number = 0.d0
	!					SEE%SEETpye=4
	!					Call Random_NUMBER(R)
	!					if (R<SEE%IonSEEcoefficient) then
	!						SEE%Number = 1
	!						Call Random_NUMBER(R)
	!						!EmissionEnergy = R*SEE%EnergyRange(4)*JtoeV
	!						!call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO,electrode)   !Random distribution
	!						call SEEAddToBundle2(SEE%EnergyRange(4),ParticleTemp,PB,PBDO,electrode) !Maxwell distribution
	!					end if
	!					PBDO%SEECountMaxOne = PBDO%SEECountMaxOne+SEE%Number
	!				end do
	!				end if
 !           return
	!end subroutine IonIndcedSEE
	!
	!subroutine SEEAddToBundle1(Energy,PO,PB,PBDO,Electrode)
	!   implicit none
	!	Real(8),intent(in) ::  Energy
	!	class(ParticleBundle),intent(inout) :: PB
 !     class(ParticleBoundaryOne),intent(inout) :: PBDO
 !     type(ParticleOne),intent(inout) :: PO
	!	Integer(4),intent(in) :: Electrode
	!	Real(8) :: velocity,Xfactor,Vfactor
	!	Xfactor=1/PB%XFactor
	!	Vfactor=1/PB%VFactor
	!	velocity = sqrt(2*Energy/ElectronMass)
	!	call PO%VelRanInit(velocity)
	!	call PO%VelRes(Vfactor)
	!	Call Random_NUMBER(R)
	!	Select case(Electrode)
	!	case(0)
	!	   PO%Vx = abs(PO%Vx)
	!		PO%X = PBDO%XMin+R*PO%Vx
	!	case(1)
	!		PO%Vx = -abs(PO%Vx)
	!		PO%X = PBDO%XMax+R*PO%Vx
	!	end Select
	!	PO%Ax = 0.0
 !     PO%Ay = 0.0
 !     PO%Az = 0.0
	!	call PB%AddOne(PO)
	!	if(Electrode==0) then
	!	  Call PBDO%PBLower%AddOne(PO)
	!	else
	!		Call PBDO%PBUpper%AddOne(PO)
	!	end if
	!	return
	!end subroutine SEEAddToBundle1
	!
	!subroutine SEEAddToBundle2(EnergyeV,PO,PB,PBDO,Electrode) !Maxewll distribution
	!   implicit none
	!	Real(8),intent(in) ::  EnergyeV
	!	class(ParticleBundle),intent(inout) :: PB
 !     class(ParticleBoundaryOne),intent(inout) :: PBDO
 !     type(ParticleOne),intent(inout) :: PO
	!	Integer(4),intent(in) :: Electrode
	!	Real(8) :: velocity,Xfactor,Vfactor
	!	Xfactor=1/PB%XFactor
	!	Vfactor=1/PB%VFactor
	!	call PO%VelMaxInit(ElectronMass,11605*EnergyeV)
	!	call PO%VelRes(Vfactor)
	!	Call Random_NUMBER(R)
	!	Select case(Electrode)
	!	case(0)
	!	   PO%Vx = abs(PO%Vx)
	!		PO%X = PBDO%XMin+R*PO%Vx
	!	case(1)
	!		PO%Vx = -abs(PO%Vx)
	!		PO%X = PBDO%XMax+R*PO%Vx
	!	end Select
	!	PO%Ax = 0.0
 !     PO%Ay = 0.0
 !     PO%Az = 0.0
	!	call PB%AddOne(PO)
	!	if(Electrode==0) then
	!	  Call PBDO%PBLower%AddOne(PO)
	!	else
	!		Call PBDO%PBUpper%AddOne(PO)
	!	end if
	!	return
	!end subroutine SEEAddToBundle2
	
	subroutine SEEEType(PO,Vfactor,SEE)
	   implicit none
		Type(ParticleOne),intent(in):: PO
		TYPE(SecondaryElectron),intent(inout) :: SEE
		real(8),intent(in) :: Vfactor
		Real(8) :: energy,angle,Emax, deltamax, Womiga,Womiga1,Womiga2,k
		Real(8) :: R,Yitai,deltav,Yitae,Xigema,Totaldelta, EnergyeV
		INTEGER(4) :: I,tempNum
		energy = PO%Energy(ElectronMass,VFactor)
		EnergyeV = energy/JtoeV
		angle = asin(PO%Vx/Dsqrt(PO%Vx*PO%Vx+PO%Vy*PO%Vy+PO%Vz*PO%Vz))
		!open (10,file="SEEdData.dat",position='APPEND')
		!DO i=0,5000
		!energy = dble(i)
		!EnergyeV = energy
		angle = 0.d0  !
		Emax = SEE%Emax0*(1+(SEE%ks*angle*angle)/pi)
		deltamax =  SEE%deltamax0*(1+SEE%ks* angle* angle/(2.0d0*pi))
        Womiga = (EnergyeV - SEE%E0)/(Emax-SEE%E0)
        Womiga1 = (EnergyeV - SEE%Ee0)/(SEE%Eemax-SEE%Ee0)  !elastic reflect electron
        Womiga2 = (EnergyeV - SEE%Eemax)/SEE%Deltae		    !elastic reflect electron
        if (Womiga<1) then
            k=0.56d0      
        else
            k=0.25d0
        end if
		deltav = Womiga*Dexp(1-Womiga)
		if(deltav<0.d0) then
        deltav = 0.d0
		end if
		deltav = deltamax*deltav**k
		Yitai = SEE%Gamai*deltav
        if ((EnergyeV<SEE%Eemax).AND.(EnergyeV>SEE%Ee0)) then
            Yitae = SEE%Gamae*deltav + SEE%Yitaemax*Womiga1*dexp(1-Womiga1)
        else
            Yitae = SEE%Gamae*deltav + SEE%Yitaemax*(1+Womiga2)*dexp(-Womiga2)
	    end if
	    Xigema = (1-Yitae-Yitai)*deltav !The real SEE induced coefficient by electron
        Totaldelta = Yitae+Yitai+Xigema
	 !  write(10,FMt="(*(es21.14,1x))") EnergyeV,deltav, Totaldelta, Yitae, Yitai,Xigema
		!end do
		!close(10)
		!stop
	    Call Random_NUMBER(R)
		IF (R<Yitai) then
			SEE%SEETpye = 3
			SEE%Number = 1
			SEE%EnergyRange(3) = energy
		else
			IF (R<(Yitai+Yitae)) then
				SEE%SEETpye = 2
				SEE%Number = 1
				SEE%EnergyRange(2) = energy
			else
				IF (R<Totaldelta) then
					SEE%SEETpye = 1
					tempNum = int(Xigema)
					Call Random_NUMBER(R)
					if (R<(Xigema-dble(tempNum))) then
							tempNum = tempNum+1
					end if
					SEE%Number = tempNum
				else
					SEE%SEETpye = 0
				end if
			end if
		end if  
		return
    end subroutine SEEEType

    !> -------------------------- ab.ZWZ 2022/3/16 ----------------------------
    Subroutine  ElectronIndcedSEEOne(PB,PO,InterPoint,electrode)
        Implicit none
        class(ParticleBundle),intent(inout) :: PB
        class(ParticleOne),intent(inout) :: PO
        Real(8),intent(in) :: InterPoint(2)
        Integer(4),intent(in) :: electrode
        
        type(ParticleOne) :: ParticleTemp
        Type(SecondaryElectron) :: SEE
        Integer(4) :: i,j,ntemp
        Real(8) :: EmissionEnergy
        
        Call ParticleTemp%Copy(PO)
        
		SEE%Number = 0.d0
		Call SEEEType(PO,PB%VFactor,SEE)
        j=0
        DO while (j<SEE%Number)
		    Select case(SEE%SEETpye)
		    case(0)
		    case(1) ! real secondary electrons(0-20eV)
			    Call Random_NUMBER(R)
                EmissionEnergy = R*20.0d0*JtoeV
			    !call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO,electrode)
                call SEEAddToBundleMaxwell(SEE%EnergyRange(1),ParticleTemp,PB,InterPoint,electrode)
		    case(2) !elastic reflected
			    EmissionEnergy = SEE%EnergyRange(2)
			    call SEEAddToBundleVRandom(EmissionEnergy,ParticleTemp,PB,InterPoint,electrode)
		    case(3) !inelastic reflected electron (0~electron energy)
			    Call Random_NUMBER(R)
			    EmissionEnergy = R*SEE%EnergyRange(3)
			    call SEEAddToBundleVRandom(EmissionEnergy,ParticleTemp,PB,InterPoint,electrode)
			    !call SEEAddToBundle2(SEE%EnergyRange(3),ParticleTemp,PB,PBDO,electrode)			 
            End Select
            j=j+1
        End Do	
        return
    end subroutine ElectronIndcedSEEOne
    
    Subroutine  IonIndcedSEEOne(PB,PO,InterPoint,electrode)
        Implicit none
        class(ParticleBundle),intent(inout) :: PB
        class(ParticleOne),intent(inout) :: PO
        Real(8),intent(in) :: InterPoint(2)
        Integer(4),intent(in) :: electrode
        type(ParticleOne) :: ParticleTemp
        Type(SecondaryElectron) :: SEE
        Integer(4) :: i,j,ntemp
        Real(8) :: Residue,EmissionEnergy
        
        Call ParticleTemp%Copy(PO)
	    SEE%Number = 0.d0
	    SEE%SEETpye=4
	    Call Random_NUMBER(R)
	    if (R<SEE%IonSEEcoefficient) then
		    SEE%Number = 1
		    Call Random_NUMBER(R)
		    EmissionEnergy = R*SEE%EnergyRange(4)*JtoeV
		    !call SEEAddToBundle1(EmissionEnergy,ParticleTemp,PB,PBDO,electrode)
		    call SEEAddToBundleMaxwell(SEE%EnergyRange(4),ParticleTemp,PB,InterPoint,electrode)
		end if
        return
	end subroutine IonIndcedSEEOne
	
    subroutine SEEAddToBundleMaxwell(EnergyeV,PO,PB,InterPoint,Electrode) !Maxwell distribution
	    implicit none
	    Real(8),intent(in) ::  EnergyeV
	    class(ParticleBundle),intent(inout) :: PB
        type(ParticleOne),intent(inout) :: PO
	    Real(8),intent(in) :: InterPoint(2)
	    Integer(4),intent(in) :: Electrode
	    Real(8) :: velocity,Xfactor,Vfactor
        Real(8) :: Vr,Vt !> ab.ZWZ
	    Xfactor=1/PB%XFactor
	    Vfactor=1/PB%VFactor
	    call PO%VelMaxInit(ElectronMass,11605*EnergyeV)
	    call PO%VelRes(Vfactor)
	    Call Random_NUMBER(R)
	    Select case(Electrode)
	    case(0)
		    PO%Vx = abs(PO%Vx)
		    PO%X = InterPoint(1)+R*PO%Vx
            PO%Y = InterPoint(2)
            PO%Z = 2*Pi*R
            PB%nGen(1) = PB%nGen(1) + 1 !> ab.ZWZ
	    case(1)
		    PO%Vx = -abs(PO%Vx)
		    PO%X = InterPoint(1)+R*PO%Vx
            PO%Y = InterPoint(2)
            PO%Z = 2*Pi*R
            PB%nGen(2) = PB%nGen(2) + 1 !> ab.ZWZ
        case(2)
            Vr =  PO%Vy*DCOS(PO%Z) + PO%Vz*DSIN(PO%Z)
            Vt = -PO%Vy*DSIN(PO%Z) + PO%Vz*DCOS(PO%Z)
            Vr = abs(Vr)
            PO%Vy = Vr*DCOS(PO%Z) - Vt*DSIN(PO%Z)
            PO%Vz = Vr*DSIN(PO%Z) + Vt*DCOS(PO%Z)
		    PO%X = InterPoint(1)
            PO%Y = InterPoint(2)+R*Vr
            PO%Z = 2*Pi*R
        case(3)
            Vr =  PO%Vy*DCOS(PO%Z) + PO%Vz*DSIN(PO%Z)
            Vt = -PO%Vy*DSIN(PO%Z) + PO%Vz*DCOS(PO%Z)
            Vr = -abs(Vr)
            PO%Vy = Vr*DCOS(PO%Z) - Vt*DSIN(PO%Z)
            PO%Vz = Vr*DSIN(PO%Z) + Vt*DCOS(PO%Z)
		    PO%X = InterPoint(1)
            PO%Y = InterPoint(2)+R*Vr
            PO%Z = 2*Pi*R
	    end Select
	    PO%Ax = 0.0
        PO%Ay = 0.0
        PO%Az = 0.0
	    call PB%AddOne(PO)
	    return
    end subroutine SEEAddToBundleMaxwell
    
    subroutine SEEAddToBundleVRandom(Energy,PO,PB,InterPoint,Electrode)
	    implicit none
	    Real(8),intent(in) ::  Energy
	    class(ParticleBundle),intent(inout) :: PB
        type(ParticleOne),intent(inout) :: PO
        Real(8),intent(in) :: InterPoint(2)
	    Integer(4),intent(in) :: Electrode
	    Real(8) :: velocity,Xfactor,Vfactor
        Real(8) :: Vr,Vt !> ab.ZWZ
	    Xfactor=1/PB%XFactor
	    Vfactor=1/PB%VFactor
	    velocity = sqrt(2*Energy/ElectronMass)
	    call PO%VelRanInit(velocity)
	    call PO%VelRes(Vfactor)
	    Call Random_NUMBER(R)
	    Select case(Electrode)
	    case(0)
		    PO%Vx = abs(PO%Vx)
		    PO%X = InterPoint(1)+R*PO%Vx
            PO%Y = InterPoint(2)
            PO%Z = 2*Pi*R
            PB%nGen(1) = PB%nGen(1) + 1 !> ab.ZWZ
	    case(1)
		    PO%Vx = -abs(PO%Vx)
		    PO%X = InterPoint(1)+R*PO%Vx
            PO%Y = InterPoint(2)
            PO%Z = 2*Pi*R
            PB%nGen(2) = PB%nGen(2) + 1 !> ab.ZWZ
        case(2)
            Vr =  PO%Vy*DCOS(PO%Z) + PO%Vz*DSIN(PO%Z)
            Vt = -PO%Vy*DSIN(PO%Z) + PO%Vz*DCOS(PO%Z)
            Vr = abs(Vr)
            PO%Vy = Vr*DCOS(PO%Z) - Vt*DSIN(PO%Z)
            PO%Vz = Vr*DSIN(PO%Z) + Vt*DCOS(PO%Z)
		    PO%X = InterPoint(1)
            PO%Y = InterPoint(2)+R*Vr
            PO%Z = 2*Pi*R
         case(3)
            Vr =  PO%Vy*DCOS(PO%Z) + PO%Vz*DSIN(PO%Z)
            Vt = -PO%Vy*DSIN(PO%Z) + PO%Vz*DCOS(PO%Z)
            Vr = -abs(Vr)
            PO%Vy = Vr*DCOS(PO%Z) - Vt*DSIN(PO%Z)
            PO%Vz = Vr*DSIN(PO%Z) + Vt*DCOS(PO%Z)
		    PO%X = InterPoint(1)
            PO%Y = InterPoint(2)+R*Vr
            PO%Z = 2*Pi*R
	    end Select
	    PO%Ax = 0.0
        PO%Ay = 0.0
        PO%Az = 0.0
	    call PB%AddOne(PO)
	    return
	end subroutine SEEAddToBundleVRandom
    
end module ModuleSEE