Module ModuleMCCInterface
    Use ModuleMCCInitialization
    Use ModuleReactionOnePegasus
    Use ModuleSpecyOne
    Use ModuleParticleBundle
    Use ModuleParticleBoundary
    Use DiagnosticsEEPF
    !|------------------------------------------------
    Use Particle_2D
    Use Constant_Variable_2D
    USE TimeControl, ONLY:dt, irestart
    USE Domain_2D
    Use Field_2D, Only:dens0
    Use Object_2D
    Use Object_Data_2D
    Implicit none
    
    Integer(4),private :: NCollision
    Type(ControlFlow),save :: ControlFlowGlobal
    Type(ParticleBundle),save,Allocatable ::  ParticleGlobal(:)
    !Type(ParticleEPF),save,Allocatable ::  ParticleEPFGlobal(:)
    
    !TYPE(ObjectType), Allocatable ::	objects(:)
    !Integer(4) :: N_objects
    
    contains
    
    Subroutine AllInitialization()
        Implicit none
        Integer(4) :: i
        
        Call InitializationControlFlow(ControlFlowGlobal)
        ControlFlowGlobal%Dt = dt*time_ref
        If (irestart==0) Then
            ControlFlowGlobal%ReStartParticles = .TRUE.
        End If
        Call GasInitPegasus(ControlFlowGlobal)      !$ SpecyGlobal and GasGlobal is initialized in the subroutine
        print *, "粒子种类数量 Ns = ", ControlFlowGlobal%Ns
        Allocate(ParticleGlobal(0:ControlFlowGlobal%Ns))
        !Allocate(ParticleEPFGlobal(0:ControlFlowGlobal%Ns))
        DO i=0,ControlFlowGlobal%Ns
            Call AllInitializationParticleBundleIFE(ParticleGlobal(i),SpecyGlobal(i),ControlFlowGlobal)
            !Call ParticleInformationPassingInit(ParticleGlobal(i),SpecyGlobal(i),i+1)
		End do
        Call MCCBundleInit(ControlFlowGlobal,SpecyGlobal,GasGlobal)
    
    End Subroutine AllInitialization
    
    Subroutine AllInitializationParticleBundleIFE(PB,SO,CF)
	    Class(ParticleBundle), intent(inout) :: PB
        Type(SpecyOne), intent(in),Target :: SO 
        Class(ControlFlow), intent(in) :: CF
        Integer(4) :: i,j,k,NPArMax
        Real(8) :: XFactor,VFactor
        Logical ::  Status
        Integer(4) :: isp
        Integer(4) :: MaxKappa=2   !1:maxwell, 0:kappa 2:poly
        Real(8) :: RegionVolume
        Real(8) :: XL,XU,YL,YU
        Integer(4) :: i_part,nParMesh=128
        isp = SO%SpecyIndex+1
        PB%dx = hx(1)*L_ref
        PB%dt = dt*time_ref
        PB%XFactor = L_ref
        PB%VFactor = v_ref
                
        PB%SO => SO     !to point at SO you have inputted
        PB%Charge = PB%SO%Charge
        PB%Mass = PB%SO%Mass          
        PB%Weight = affp_bjw(isp)
        
        !======= for wsy paper case ==========
        PB%NParNormal = 6000000
        NPArMax = PB%NParNormal
        If(Allocated(PB%PO)) Deallocate(PB%PO)
        Allocate(PB%PO(NPArMax))
        !======================================
        
        If (PB%UnequalWeightFlag .and. delta_global == 1) Then
            nParMesh = 100
            PB%NPar = nParMesh * (nx-1) * (ny-1)
            PB%NParNormal = 120 * (nx-1) * (ny-1)
            !PB%NPar = 2 * (nx-1) * (ny-1)
            !PB%NParNormal = 1 * (nx-1) * (ny-1)
            NPArMax = Ceiling(2.0 * PB%NParNormal)
            If(Allocated(PB%PO)) Deallocate(PB%PO)
            Allocate(PB%PO(NPArMax))
            RegionVolume = hx(1)*L_ref * PI*hx(2)**2 *L_ref**2
            PB%WQOne = dens0(isp)*n_ref * RegionVolume / nParMesh
            i_part = 0
            If (CF%ReStartParticles) Then
                Do j=1,ny-1
                    Do i=1,nx-1
                !Do i=1,nx-1
                !    Do j=1,ny-1
                        XL=(i-1)*hx(1)
                        XU=i    *hx(1)
                        YL=(j-1)*hx(2)
                        YU=j    *hx(2)
                        !RegionVolume = (XU-XL)*L_ref * PI*(YU**2-YL**2) *L_ref**2
                        Do k=1,nParMesh
                            i_part = i_part + 1
                            !PB%PO(i_part)%WQ = dens0(isp)*n_ref * RegionVolume / nParMesh
                            PB%PO(i_part)%WQ = (2*j-1) * PB%WQOne
                            PB%LXScaled = .True.
                            Call PB%PO(i_part)%PosInit(XL,XU,YL,YU)
                            PB%LVScaled = .True.
                            Call PB%PO(i_part)%VelMaxInit(PB%SO%Mass,PB%SO%InitTemperature)
                            VFactor = 1.d0 / PB%VFactor
                            Call PB%PO(i_part)%VelRes(VFactor)       ! normalization
                            Call PB%PO(i_part)%AccInpInit()
                        Enddo
                    Enddo
                Enddo
            Else
                Call LoadParticleBundle(PB,Status)
                If (Status) Then
                     Do j=1,ny-1
                        Do i=1,nx-1
                            XL=(i-1)*hx(1)
                            XU=i    *hx(1)
                            YL=(j-1)*hx(2)
                            YU=j    *hx(2)
                            RegionVolume = (XU-XL)*L_ref * PI*(YU**2-YL**2) *L_ref**2
                            Do k=1,nParMesh
                                i_part = i_part + 1
                                PB%PO(i_part)%WQ = dens0(isp)*n_ref * RegionVolume / nParMesh
                                PB%LXScaled = .True.
                                Call PB%PO(i_part)%PosInit(XL,XU,YL,YU)
                                PB%LVScaled = .True.
                                Call PB%PO(i_part)%VelMaxInit(PB%SO%Mass,PB%SO%InitTemperature)
                                VFactor = 1.d0 / PB%VFactor
                                Call PB%PO(i_part)%VelRes(VFactor)       ! normalization
                                Call PB%PO(i_part)%AccInpInit()
                            Enddo
                        Enddo
                    Enddo
                Else
                    XFactor = 1.d0 / PB%XFactor
                    Call PB%PosRes()
                    VFactor = 1.d0 / PB%VFactor
                    Call PB%VelRes()
                End If
            End If         
        Else
            If (delta_global == 0) Then
                RegionVolume = (dxmax-dxmin)*L_ref * (dymax-dymin) *L_ref
            Elseif (delta_global == 1) Then
                RegionVolume = (dxmax-dxmin)*L_ref * PI*(dymax-dymin)**2 *L_ref**2
            Endif
            !PB%NParNormal = CF%ParticlePerGrid * (CF%Nx-1) * PB%SO%InitDensity
            !PB%NPar = dens0(isp)*n_ref * RegionVolume / affp_bjw(isp)      
            
           PB%NPar = 5000 * (dxmaxmax-dxminmin) * (dymaxmax-dyminmin)     
            PB%Weight = affp_bjw(isp)
            !PB%Weight = dens0(isp)*n_ref * RegionVolume / PB%NPar 
            !print*,dens0(isp)*n_ref 
            PB%NParNormal = 5000 * (dxmaxmax-dxminmin) * (dymaxmax-dyminmin) ! wsy revise cause insufficient virtual memory   PB%NParNormal = 2048 * (nx-1) * (ny-1)
            
        
        
            NPArMax = Ceiling(3.0 * PB%NParNormal)
            If(Allocated(PB%PO)) Deallocate(PB%PO)
            Allocate(PB%PO(NPArMax))
            If (CF%ReStartParticles) Then
                Do i=1,PB%Npar
                    PB%LXScaled = .True.
                    !Call PB%PO(i)%PosInit(Dble(CF%NxL),Dble(CF%NxU))
                    Call PB%PO(i)%PosInit(dxminmin,dxmaxmax,dyminmin,dymaxmax) !改初始区域改这
                    
                    If (PB%PO(i)%x>dxmax) then
                       Call PB%DelOne(i)
                   End if
                    PB%LVScaled = .True.
                    if(MaxKappa==1)then 
                        Call PB%PO(i)%VelMaxInit(PB%SO%Mass,PB%SO%InitTemperature)
                        VFactor = 1.d0 / PB%VFactor
                        Call PB%PO(i)%VelRes(VFactor)       ! normalization
                    elseif (MaxKappa==0) then 
                        VFactor = 1.d0 / PB%VFactor
                        Call PB%PO(i)%VelKappaInit(PB%SO%Mass,PB%SO%InitTemperature,VFactor)
                    else 
                        VFactor = 1.d0 / PB%VFactor
                        Call PB%PO(i)%VelPolyInit(PB%SO%Mass,PB%SO%InitTemperature,VFactor)
                    endif 
                    Call PB%PO(i)%AccInpInit()
                    Call PB%PO(i)%ParLocate(0)
                End Do
            Else
                Call LoadParticleBundle(PB,Status)
                If (Status) Then
                    Do i=1,PB%NPar
                        PB%LXScaled = .True.
                        Call PB%PO(i)%PosInit(Dble(CF%NxL),Dble(CF%NxU))
                        PB%LVScaled = .True.
                        if(MaxKappa==1)then 
                            Call PB%PO(i)%VelMaxInit(PB%SO%Mass,PB%SO%InitTemperature)
                            VFactor = 1.d0 / PB%VFactor
                            Call PB%PO(i)%VelRes(VFactor)       ! normalization
                        else 
                            VFactor = 1.d0 / PB%VFactor
                            Call PB%PO(i)%VelKappaInit(PB%SO%Mass,PB%SO%InitTemperature,VFactor)
                        endif 
                        Call PB%PO(i)%AccInpInit()
                    End Do
                Else
                    XFactor = 1.d0 / PB%XFactor
                    Call PB%PosRes()
                    VFactor = 1.d0 / PB%VFactor
                    Call PB%VelRes()
                End If
            End If        
        Endif
    End Subroutine AllInitializationParticleBundleIFE
    
    
    Subroutine AdjustParticleBoundary(PB,PO,TimeMove,IvelFlag,IposFlag,N_objects,objects,delta,OriPosi,detaV)
        Use Domain_2D, Only: dxmin, dymin, dxmax, dymax
        Implicit none
        Class(ParticleBundle),intent(inout) :: PB
        Type(ParticleOne),intent(inout) :: PO
        Real(8),intent(inout) :: TimeMove
        Integer(4),intent(inout) :: IvelFlag, IposFlag
        Integer(4),intent(in) :: delta
        Integer(4), intent(in) :: N_objects
        Type(ObjectType), intent(in) :: objects(:)
        Real(8),intent(in) :: OriPosi(3)
        Real(8),intent(in),optional :: detaV(3)
   
        Integer(4) :: iBoundary,i
        Real(8) :: LineEnds(2,2)
        Real(8) :: Trail_2D(2,2)
        !Real(8) :: Trail_2Da(2,3)   !/z,r,theta/
        !Real(8) :: Trail_3D(2,3)    !/x,y,z/
        Real(8) :: InterPoint(2),InterPointTemp(2)
        Integer(4) :: InterFlag
        Real(8) :: x,y,z,r,theta
        
        Real(8) :: TimeRemain, TimeTemp
        Integer(4) :: NBoundary = 8, CrossBoundary
        Integer(4) :: I_Count
        
        double precision :: ranf1, ranf2, ranf3
        REAL(8) :: beta, velocity_tangential,theta_v,VFactor,Kappa,gama
        Real(8) :: KB=1.3807d-23
        REAL(8), PARAMETER	::	pii	= 3.14159265358979D0

        Trail_2D(1,:) = (/OriPosi(1),OriPosi(2)/)
        Trail_2D(2,:) = (/PO%X,PO%Y/)
        beta=PB%Mass/(2*KB*PB%SO%InitTemperature)
        PB%VFactor = v_ref
        Kappa=2.0
        gama=3.0
        
        !> ----------------------------------
        !I_Count = 0
        !TimeRemain = 0.
        !Do iBoundary = 1, NBoundary
        !    Select Case(iBoundary)
        !    Case(1)
        !        LineEnds(1,:) = (/dxmin,dymin/)
        !        LineEnds(2,:) = (/dxmin,dymax/) 
        !    Case(2)
        !        LineEnds(1,:) = (/dxmax,dymin/)
        !        LineEnds(2,:) = (/dxmax,dymax/)
        !    Case(3)
        !        LineEnds(1,:) = (/dxmax,dymin/)
        !        LineEnds(2,:) = (/dxmin,dymin/)
        !    Case(4)
        !        LineEnds(1,:) = (/dxmax,dymax/)
        !        LineEnds(2,:) = (/dxmin,dymax/)
        !    Case(5)
        !        LineEnds(1,:) = (/objects(1)%Locations(2,1),objects(1)%Locations(1,2)/)
        !        LineEnds(2,:) = (/objects(3)%Locations(2,1),objects(3)%Locations(2,2)/)
        !    Case(6)
        !        LineEnds(1,:) = (/objects(5)%Locations(1,1),objects(5)%Locations(2,2)/)
        !        LineEnds(2,:) = (/objects(5)%Locations(1,1),objects(5)%Locations(1,2)/)
        !    Case(7)
        !        LineEnds(1,:) = (/objects(3)%Locations(2,1),objects(3)%Locations(2,2)/)
        !        LineEnds(2,:) = (/objects(4)%Locations(2,1),objects(4)%Locations(1,2)/)
        !    Case(8)
        !        LineEnds(1,:) = (/objects(5)%Locations(2,1),objects(5)%Locations(2,2)/)
        !        LineEnds(2,:) = (/objects(5)%Locations(1,1),objects(5)%Locations(2,2)/)
        !    End Select
        !    Call LineIntersection(LineEnds,Trail_2D,InterPointTemp,InterFlag)
        !    If (InterFlag == 1) Then
        !        I_Count = I_Count + 1
        !        If (present(detaV) .and. detaV(1)/=0.) Then
        !            TimeTemp = (PO%X - InterPointTemp(1)) / detaV(1)
        !        Elseif (PO%Vx /= 0.) Then
        !            TimeTemp = (PO%X - InterPointTemp(1)) / PO%Vx
        !        Else
        !            print*,'Vx = 0 !!'
        !            pause
        !        Endif
        !        If (TimeTemp > TimeRemain) Then ! only record the first boundary the particle intersected
        !            TimeRemain = TimeTemp
        !            CrossBoundary = iBoundary
        !            InterPoint = InterPointTemp
        !        End If
        !    Endif
        !End Do
        !
        !If (I_Count > 0) Then
        !    Select Case(CrossBoundary)
        !    Case(1)
        !        PB%nLoss(1) = PB%nLoss(1) + 1 !> ab.ZWZ
        !        PO%X = -2000 
        !    Case(2)
        !        PB%nLoss(2) = PB%nLoss(2) + 1 !> ab.ZWZ
        !        PO%X = -2000
        !    Case(3)
        !        PB%nLoss(3) = PB%nLoss(3) + 1 !> ab.ZWZ
        !        If (delta == 0) Then
        !            ! Abortion boundary -- delete
        !            PO%X = -2000
        !        Elseif (delta == 1) Then
        !            Print*,'wrong cross'
        !            Print*,PO
        !            pause
        !        Endif
        !    Case(4)
        !        PB%nLoss(4) = PB%nLoss(4) + 1 !> ab.ZWZ
        !        ! Abortion boundary -- delete
        !        PO%X = -2000
        !    Case(5)
        !        PO%X = -2000
        !        If (InterPoint(2) < objects(2)%Locations(1,2) .or. InterPoint(2) > objects(2)%Locations(2,2)) Then
        !            If (PB%SO%SpecyIndex == 0) Then
        !                Call ElectronIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,0)
        !            Else
        !                Call IonIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,0)
        !            End If
        !        End If
        !    Case(6)
        !        PO%X = -2000
        !        If (PB%SO%SpecyIndex == 0) Then
        !            Call ElectronIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,1)
        !        Else
        !            Call IonIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,1)
        !        End If
        !    Case(7)
        !        PO%X = -2000
        !        If (PB%SO%SpecyIndex == 0) Then
        !            Call ElectronIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,3)
        !        Else
        !            Call IonIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,3)
        !        End If
        !    Case(8)
        !        PO%X = -2000
        !        If (PB%SO%SpecyIndex == 0) Then
        !            Call ElectronIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,2)
        !        Else
        !            Call IonIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,2)
        !        End If
        !    End Select
        !End if
        !
        !> ----------------------------------
        
        If (PO%X <= dxmin) Then
            
            !Abortion boundary -- delete
            !PB%nLoss(1) = PB%nLoss(1) + 1 !> ab.ZWZ
            !PO%X = -2000
            !! Reflex boundary -- rebounce
                !Iposflag = 1
                !TimeRemain = (PO%X - dxmin)/PO%Vx
                !TimeMove = TimeRemain
                !PO%Y = PO%Y - PO%Vy * TimeRemain
                !PO%Vx = -PO%Vx
                !PO%X = dxmin + 10E-6
            !diffuse reflection zzj24/10/15 沈青 稀薄气体动力学3.2 Maxwellian
            !Iposflag = 1
            !TimeRemain = (PO%X - dxmin)/PO%Vx
            !TimeMove = TimeRemain
            !PO%Y = PO%Y - PO%Vy * TimeRemain
            !PO%X = dxmin + 10E-6
            !call DRandom(ranf1)
            !call DRandom(ranf2)
            !call DRandom(ranf3)
            !PO%Vx =SQRT ((-DLOG (ranf1))/beta)
            !velocity_tangential=SQRT ((-DLOG (ranf2))/beta)
            !theta_v=2*pii*ranf3
            !PO%Vy=velocity_tangential*COS (theta_v)
            !PO%Vz=velocity_tangential*SIN (theta_v)
            !VFactor = 1.0 / PB%VFactor
            !call PO%VelRes(VFactor)
            
            !kappa漫反射
            Iposflag = 1
            TimeRemain = (PO%X - dxmin)/PO%Vx
            TimeMove = TimeRemain
            PO%Y = PO%Y - PO%Vy * TimeRemain
            PO%X = dxmin + 10E-6
            !call PO%VelKappaInit(PB%SO%Mass,PB%SO%InitTemperature,PB%VFactor)
            call DRandom(ranf1)
            CALL RANDOM_NUMBER(ranf1)
            PO%Vx =SQRT ((1.d0*Kappa-1.5)/beta*((ranf1)**(-1/(Kappa-1))-1))
            VFactor = 1.0 / PB%VFactor
            call PO%VelRes(VFactor)
            !poly漫反射
            !Iposflag = 1
            !TimeRemain = (PO%X - dxmin)/PO%Vx
            !TimeMove = TimeRemain
            !PO%Y = PO%Y - PO%Vy * TimeRemain
            !PO%X = dxmin + 10E-6
            !call RANDOM_NUMBER (ranf1)
            !PO%Vx=SQRT (gama/((gama-1)*beta)*(1-(ranf1)**((2.0*gama-2)/(gama+1))))
            !VFactor = 1.0 / PB%VFactor
            !call PO%VelRes(VFactor)
            
            !LineEnds(1,:) = (/dxmin,dymin/)
            !LineEnds(2,:) = (/dxmin,dymax/)
            !Call LineIntersection(LineEnds,Trail_2D,InterPoint,InterFlag)
            !If (InterFlag == 1) Then
            !    If (PB%SO%SpecyIndex == 0) Then
            !        Call ElectronIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,0)
            !    Else
            !        Call IonIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,0)
            !    End If
            !Endif
        Elseif (PO%X > dxmax) Then
            PB%nLoss(2) = PB%nLoss(2) + 1 !> ab.ZWZ
            PO%X = -2000
            !LineEnds(1,:) = (/dxmax,dymin/)
            !LineEnds(2,:) = (/dxmax,dymax/)
            !Call LineIntersection(LineEnds,Trail_2D,InterPoint,InterFlag)
            !Iposflag = 1
            !    TimeRemain = (PO%X - dxmin)/PO%Vx
            !    TimeMove = TimeRemain
            !    PO%Y = PO%Y - PO%Vy * TimeRemain
            !    PO%Vx = -PO%Vx
            !    PO%X = dxmax - 10E-6
            !If (InterFlag == 1) Then
            !    If (PB%SO%SpecyIndex == 0) Then
            !        Call ElectronIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,1)
            !    Else
            !        Call IonIndcedSEEOne(ParticleGlobal(0),PO,InterPoint,1)
            !    End If
            !Endif
        Elseif (PO%Y < dymin) Then
            !PB%nLoss(3) = PB%nLoss(3) + 1 !> ab.ZWZ
            If (delta == 0) Then
                ! Abortion boundary -- delete
                !PO%X = -2000
                !! Reflex boundary -- rebounce
                Iposflag = 1
                TimeRemain = (PO%Y - dymin)/PO%Vy
                TimeMove = TimeRemain
                PO%X = PO%X - PO%Vx * TimeRemain
                PO%Vy = -PO%Vy
                PO%Y = dymin + 10E-6
            
                If (PO%Y > dymax) Then
                    PO%Y = PO%Y -4.0
                    PRINT*, 'Case5: Check Partition and Code, Stop'
                    Pause
                Endif
            Elseif (delta == 1) Then
                Print*,'wrong cross'
                Print*,PO
                pause
            Endif
        Elseif (PO%Y > dymax) Then
            !PB%nLoss(4) = PB%nLoss(4) + 1 !> ab.ZWZ
            !! Abortion boundary -- delete
            !PO%X = -2000
            
            ! Reflex boundary -- rebounce
            If (delta == 0) Then
                Iposflag = 1
                TimeRemain = (PO%Y - dymax)/PO%Vy
                TimeMove = TimeRemain
                PO%X = PO%X - PO%Vx * TimeRemain
                PO%Vy = -PO%Vy
                PO%Y = dymax - 10E-6
                
                If (PO%Y < dymin) Then
                    PRINT*, 'Case6: Check Partition and Code, Stop'
                    Pause
                Endif
            Endif
        Endif
        
        If (N_objects /= 0) Then
            Do i = 1, N_objects
                ! asorb objects
                If (objects(i)%Shape == 3) Then
                    If (PO%X > objects(i)%Locations(1,1) .And. PO%X < objects(i)%Locations(2,1) .AND. &
                        PO%Y > objects(i)%Locations(1,2) .AND. PO%Y < objects(i)%Locations(2,2) ) THEN  
                        PO%X = -2000
                    End If
                Elseif (objects(i)%Shape == 1) Then
                    If ( (PO%X - objects(i)%Locations(1,1))**2 + (PO%Y - objects(i)%Locations(1,2))**2 < objects(i)%Dimensions(1)**2) Then
                        PO%X = -2000
                    End If
                End If
            End Do
            
            !# 1
            !LineEnds(1,:) = (/objects(1)%Locations(2,1),objects(1)%Locations(2,2)/)
            !LineEnds(2,:) = (/objects(1)%Locations(1,1),objects(1)%Locations(2,2)/)
            !Call LineIntersection(LineEnds,Trail_2D,InterPoint,InterFlag)
            !If (InterFlag == 1) Then
            !    PO%X = -2000
            !Endif
            !
            !!# 2
            !LineEnds(1,:) = (/objects(1)%Locations(2,1),objects(1)%Locations(1,2)/)
            !LineEnds(2,:) = (/objects(1)%Locations(2,1),objects(1)%Locations(2,2)/)
            !Call LineIntersection(LineEnds,Trail_2D,InterPoint,InterFlag)
            !If (InterFlag == 1) Then
            !    PO%X = -2000
            !Endif
            !
            !!# 3
            !LineEnds(1,:) = (/objects(2)%Locations(1,1),objects(2)%Locations(2,2)/)
            !LineEnds(2,:) = (/objects(2)%Locations(1,1),objects(2)%Locations(1,2)/)
            !Call LineIntersection(LineEnds,Trail_2D,InterPoint,InterFlag)
            !If (InterFlag == 1) Then
            !    PO%X = -2000
            !Endif
            !
            !!# 4
            !LineEnds(1,:) = (/objects(2)%Locations(2,1),objects(2)%Locations(2,2)/)
            !LineEnds(2,:) = (/objects(2)%Locations(1,1),objects(2)%Locations(2,2)/)
            !Call LineIntersection(LineEnds,Trail_2D,InterPoint,InterFlag)
            !If (InterFlag == 1) Then
            !    PO%X = -2000
            !Endif
        Endif
            
    End Subroutine AdjustParticleBoundary 
             
    Subroutine ParticleInformationPassingInit(PB,SO,isp)
        Implicit none
        Type(ParticleBundle),intent(inout) :: PB
        Type(SpecyOne), intent(in),Target :: SO 
        Integer(4),intent(in) :: isp            
        
        PB%XFactor = L_ref
        PB%VFactor = v_ref
        
        PB%Charge = qs(isp)
        PB%Mass = xm(isp)
        PB%Weight = affp_bjw(isp)
        
        PB%Dx = hx(1)*L_ref
        PB%Dt = dt*time_ref
        
        !PB%XFactor = PB%dx
        !PB%VFactor = PB%dx / PB%dt
        
        PB%SO => SO
    
    End Subroutine  ParticleInformationPassingInit
     
    !|----------------------------------------------------------------------
     Subroutine MCCEntrance()   
        Implicit none
        print*,'before MCC:'
        print*,'ntot=',sum(ParticleGlobal(0:ControlFlowGlobal%Ns)%Npar),'ns(1)=',ParticleGlobal(0)%Npar,'ns(2)=',ParticleGlobal(1)%Npar
        !Call AllParticleLocal2MCC()
        Call MCC(ControlFlowGlobal%Ns,ControlFlowGlobal%Ng,ParticleGlobal,SpecyGlobal,GasGlobal,MCCBundleGlobal) 
        !Call AllParticleMCC2Local()
        print*,'after MCC:'
        print*,'ntot=',sum(ParticleGlobal(0:ControlFlowGlobal%Ns)%Npar),'ns(1)=',ParticleGlobal(0)%Npar,'ns(2)=',ParticleGlobal(1)%Npar
     End Subroutine  MCCEntrance
    !|----------------------------------------------------------------------
     
     
    !> added by ZWZ -- Move all particles
    Subroutine MoveAll(PB,CF,N_objects,objects,delta,PushFlag)
        Use Domain_2D, Only:dxmin
        Implicit none
        Class(ParticleBundle), intent(inout) :: PB
        Class(ControlFlow), intent(in) :: CF
        Integer(4), intent(in) :: N_objects
        Type(ObjectType), intent(in) :: objects(:)
        integer(4), intent(in) :: delta
        integer(4), intent(in) :: PushFlag
        Integer(4) :: i,isp
        Integer(4) :: IvelFlag = 0, IposFlag = 0
        Real(8) :: TimeMove
        Real(8) :: OriPosi(3)
        Real(8) :: detaV(3)
       
        PB%nLoss(1) =0
        isp = PB%SO%SpecyIndex
        Do i = 1, PB%Npar
            
            ! ========Turner Case=============
            !TimeMove = CF%dt/time_ref
            ! ================================
            
            TimeMove = dt
            
            IvelFlag = 1
            IposFlag = 1
            Do While ( IvelFlag == 1 .or. IposFlag == 1)
                Select Case(PushFlag)
                Case(0) ! explicit move
                    Call MoveOne(PB%PO(i),isp,TimeMove,IvelFlag, IposFlag, delta, OriPosi)
                    Call AdjustParticleBoundary(PB,PB%PO(i),TimeMove,IvelFlag,IposFlag,N_objects,objects,delta,OriPosi)
                Case(1) ! implicit prepush
                    Call PrePushOne(PB%PO(i),isp,TimeMove,IvelFlag, IposFlag, delta, OriPosi)
                    Call AdjustParticleBoundary(PB,PB%PO(i),TimeMove,IvelFlag,IposFlag,N_objects,objects,delta,OriPosi)
                Case(2) ! implicit postpush
                    Call PostPushOne(PB%PO(i),isp,TimeMove,IvelFlag, IposFlag, delta, OriPosi,detaV)
                    Call AdjustParticleBoundary(PB,PB%PO(i),TimeMove,IvelFlag,IposFlag,N_objects,objects,delta,OriPosi,detaV)
                End Select
                !Call AdjustParticleBoundary(PB,PB%PO(i),TimeMove,IvelFlag,IposFlag,N_objects,objects,delta,OriPosi)
            End Do
        End Do
        
        Do i=PB%NPar,1,-1
            If (PB%PO(i)%X<=dxmin) then
                Call PB%DelOne(i)
            End if
            
            !zzj delete0822
            !If (isp > 0 .And. PB%PO(i)%X > 795) Then
            !    write(*,*) 'particle reach 395 2'
            !    Call SYSTEM_CLOCK(Time2)
            !    Write(*,*) 'RunTime = ', Time2 - Time1
            !    pause
            !End If

        End do
        
        !Do i=PB%NPar,1,-1
        !    If (PB%PO(i)%X<=dxmin .or. PB%PO(i)%X>dxmax) then
        !        print*,i
        !        print*, PB%PO(i)
        !        pause
        !    End if
        !End do
                
    End subroutine MoveAll
     
        
    SUBROUTINE	Quasi_cdp(PB_e,PB_i,delta,it)

        !USE PIC_MAIN_PARAM_2D
          !$ ab.ZWZ for checking atom distribution 
        IMPLICIT NONE
        Class(ParticleBundle), intent(inout) :: PB_e,PB_i
        INTEGER,INTENT(IN) :: delta !$ ab.ZWZ
        INTEGER       :: ii, jj, i_part, ntotp
        !REAL(8)       :: ranum
        double precision ::  ranum
        REAL(8)       :: v_r, FS1, FS2, QA, UU, UN, GG, V_x, V_y, V_z
        REAL(8)		  :: ns_local(2)
        INTEGER       :: add_number(2)
        REAL(8)       ::  b_amb(2,2)
        Real(8),parameter :: PI_1=3.141592653589793238D0



        INTEGER, DIMENSION(2,(nx-1)*(ny-1))		::	ic1, IC
        INTEGER									::	n_cell_in
        INTEGER									::	insp,i,nnx,nny
        INTEGER									::	ipre, jpre
        INTEGER, DIMENSION(nx-1,ny-1)			::	net_at_cell
        REAL(8), DIMENSION(nx-1,ny-1)			::	pro_at_cell
        INTEGER									::	net_at_boundary = 0
        INTEGER									::	mc, num_ion, num_elec, num, ne_old


        INTEGER     :: it
        CHARACTER*40	fname
        INTEGER, DIMENSION(nx-1)    :: nab

        print*,'before quasi:'
        print*,'ns(1)=',ns(1)
        print*,'ns(2)=',ns(2)

        !!! ******************************* 方法一（bjw 2019-5-7）: 各个网格准中性条件 ***********************************************
        ne_old = ns(1) 

        !n_cell_in = 10 / hx(2)  !!! 准中性区域的网格厚度
        n_cell_in = 10 / hx(1)  !!! 准中性区域的网格厚度

        nnx = nx - 1
        nny = ny - 1

        DO insp = 1, 2
	        CALL INDEXM_cdp(PB_e,PB_i,insp, IC, nnx, nny)  !!! nnx, nny网格数
	        ic1(insp,:) = IC(2,:)       !$ 每个网格的粒子数
        ENDDO
        !write(*,*) ic1
        !print*,'before quasi11111111111111:'
        DO ipre=1, 1 + n_cell_in  !!!x方向入射(x最小位置处)
	        net_at_cell = 0
	        net_at_boundary = 0
	        pro_at_cell = 0
	        DO jpre = 1, nny
		        mc = (ipre - 1) * nny + jpre
		        net_at_cell(ipre,jpre) = ic1(2,mc) - ic1(1,mc)          !$ 离子数-电子数
                !print*,'mc = ',mc
                !print*,'ic1(2,mc)=',ic1(2,mc)
                !print*,'ic1(1,mc)=',ic1(1,mc)
                !pause
		        net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
	        ENDDO

	        IF(net_at_boundary > 0) THEN
		        net_at_boundary = 0
		        DO jpre = 1, nny
			        IF(net_at_cell(ipre,jpre) > 0) THEN
				        pro_at_cell(ipre,jpre) = REAL(net_at_cell(ipre,jpre))
				        net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
			        ELSE
				        pro_at_cell(ipre,jpre) = 0.
			        ENDIF
		        ENDDO

		        !pro_at_cell(:,jpre) = pro_at_cell(:,jpre) / REAL(net_at_boundary)
                pro_at_cell(ipre,:) = pro_at_cell(ipre,:) / REAL(net_at_boundary)

		        DO jpre = 2, nny
			        pro_at_cell(ipre,jpre) = pro_at_cell(ipre,jpre) + pro_at_cell(ipre,jpre-1)      
		        ENDDO

		        DO i = 1, net_at_boundary
			        num_elec = ns(1)
			        num_ion = ns(2)
			        num = num_elec + num_ion

                    CALL DRandom(ranum)
			        jpre = 1
			        DO WHILE(ranum > pro_at_cell(ipre, jpre) .AND. jpre < nnx)
				        jpre = jpre + 1
			        ENDDO

			        CALL DRandom(ranum)
			        !part(num+1,1) = 0. + (ipre - ranum) * hx(1)			
			        PB_e%PO(num_elec+1)%X = f_left_wall(1) + (ipre - ranum) * hx(1)			                !$ x位置
			        CALL DRandom(ranum)
			        !part(num+1,2) = 0. + (jpre - ranum) * hx(2)
			        !part(num+1,2) = f_left_wall(2) + (jpre - ranum) * hx(2)
                   !$ ========================= mb.ZWZ ================================= \\
                    IF(delta == 0) THEN
                        PB_e%PO(num_elec+1)%Y= f_left_wall(2) + (jpre - ranum) * hx(2)                      !$ y位置
                        PB_e%PO(num_elec+1)%Z= 0.
                    ELSEIF(delta == 1) THEN
                        PB_e%PO(num_elec+1)%Y= f_left_wall(2) + SQRT(hx(2)**2*((jpre -1)*(jpre -1) &
                                        +((jpre)+(jpre -1)) &
                                        *((jpre)-(jpre -1))*ranum))
                        CALL DRandom(ranum)
                        PB_e%PO(num_elec+1)%Z = 2.*PI*ranum
                        !part(num+1,2)= f_left_wall(2) + SQRT(hx(2)**2*(jpre*jpre &
                        !                -((jpre)+(jpre -1)) &
                        !                *((jpre)-(jpre -1))*ranum))
                        !part(num+1,2)= f_left_wall(2) + (jpre - ranum) * hx(2)
                        !rweight(num+1) = part(num+1,2)*dy_inject(ipf(jj))
                        !YINI(num+1) = part(num+1,2)
                    ENDIF
                    !$ ========================= mb.ZWZ ================================= //
			

        1000        CALL Loadv(V_x, tmpj(ipf(1)), 1)   !!! 只补电子，所以用 temp0(ipf(1))
		            CALL Loadv(V_y, tmpj(ipf(1)), 1)
		            CALL Loadv(V_z, tmpj(ipf(1)), 1)
			        IF(V_x <= 0) THEN
                        !WRITE(*,*) '++++++++++++'
				        goto 1000
			        ENDIF
			        PB_e%PO(num_elec+1)%Vx = V_x
			        PB_e%PO(num_elec+1)%Vy = V_y
			        PB_e%PO(num_elec+1)%Vz = V_z

                    PB_e%PO(num_elec+1)%ax = 0.
			        PB_e%PO(num_elec+1)%ay = 0.
			        PB_e%PO(num_elec+1)%az = 0.
            
			        !part(num+1,7) = 1

			        ns(1) = ns(1) + 1
                    PB_e%NPar=PB_e%NPar+1
            
                    !!$----------- ab.ZWZ for checking atom distribution ----------
                    !IF (DEBUG==1) THEN
                    !    PART_STAT_COUNT = PART_STAT_COUNT + 1
                    !    IF (PART_STAT_COUNT<=PART_STAT_NUM) THEN
                    !        !$----------------------------------------------- 
                    !        PART_STAT_XP(PART_STAT_COUNT) = part(num+1,1)
                    !        PART_STAT_YP(PART_STAT_COUNT) = part(num+1,2)
                    !        PART_STAT_ZP(PART_STAT_COUNT) = part(num+1,3)
                    !        !$-----------------------------------------------
                    !        PART_STAT_VX(PART_STAT_COUNT) = part(num+1,4)
                    !        PART_STAT_VY(PART_STAT_COUNT) = part(num+1,5)
                    !        PART_STAT_VZ(PART_STAT_COUNT) = part(num+1,6)
                    !    ELSE
                    !        PART_STAT_COUNT = PART_STAT_NUM
                    !    END IF
                    !END IF
                    !!$---------- ab.ZWZ for checking atom distribution ----------    
		        ENDDO
	        ENDIF
        ENDDO
        WRITE(6,*) '# Quasi1: to inject on one side',ns(1)-ne_old
        DO ipre=nnx - n_cell_in, nnx  !!!x方向入射(x最大位置处)
	        net_at_cell = 0
	        net_at_boundary = 0
	        pro_at_cell = 0
	        DO jpre = 1, nny
		        mc = (ipre - 1) * nny + jpre
		        net_at_cell(ipre,jpre) = ic1(2,mc) - ic1(1,mc)
		        net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
	        ENDDO
        !write(*,*) ic1(2,1) - ic1(1,1)
	        IF(net_at_boundary > 0) THEN
		        net_at_boundary = 0
		        DO jpre = 1, nny
			        IF(net_at_cell(ipre,jpre) > 0) THEN
				        pro_at_cell(ipre,jpre) = REAL(net_at_cell(ipre,jpre))
				        net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
			        ELSE
				        pro_at_cell(ipre,jpre) = 0.
			        ENDIF
		        ENDDO

		        !pro_at_cell(:,jpre) = pro_at_cell(:,jpre) / REAL(net_at_boundary)
                pro_at_cell(ipre,:) = pro_at_cell(ipre,:) / REAL(net_at_boundary)

		        DO jpre = 2, nny
			        pro_at_cell(ipre,jpre) = pro_at_cell(ipre,jpre) + pro_at_cell(ipre,jpre-1)
		        ENDDO

		        DO i = 1, net_at_boundary
			        num_elec = ns(1)
			        num_ion = ns(2)
			        num = num_elec + num_ion
            
                    CALL DRandom(ranum)
			        jpre = 1
			        DO WHILE(ranum > pro_at_cell(ipre, jpre) .AND. jpre < nnx)
				        jpre = jpre + 1
			        ENDDO

			        CALL DRandom(ranum)
                    PB_e%PO(num_elec+1)%X=0. + (ipre - ranum) * hx(1)
			        !part(num+1,1) = 0. + (ipre - ranum) * hx(1)			
			        CALL DRandom(ranum)
                    !$ ========================= mb.ZWZ ================================= \\
                    IF(delta == 0) THEN
                        PB_e%PO(num_elec+1)%Y=f_left_wall(2) + (jpre - ranum) * hx(2)
                        PB_e%PO(num_elec+1)%Z=0.
                        !part(num+1,2)= f_left_wall(2) + (jpre - ranum) * hx(2)                      !$ y位置
                        !part(num+1,3) = 0.
                    ELSEIF(delta == 1) THEN
                
                        PB_e%PO(num_elec+1)%Y=f_left_wall(2) + SQRT(hx(2)**2*((jpre -1)*(jpre -1) &
                                        +((jpre)+(jpre -1)) &
                                        *((jpre)-(jpre -1))*ranum))
                        !part(num+1,2)= f_left_wall(2) + SQRT(hx(2)**2*((jpre -1)*(jpre -1) &
                        !                +((jpre)+(jpre -1)) &
                        !                *((jpre)-(jpre -1))*ranum))
                        CALL DRandom(ranum)
                        PB_e%PO(num_elec+1)%Z = 2.*PI_1*ranum
                    ENDIF
                    !$ ========================= mb.ZWZ ================================= //

        2000        CALL Loadv(V_x, tmpj(ipf(1)), 1)   !!! 只补电子，所以用 temp0(ipf(1))
		            CALL Loadv(V_y, tmpj(ipf(1)), 1)
                    CALL Loadv(V_z, tmpj(ipf(1)), 1)
			        IF(V_x >= 0) THEN
                        !WRITE(*,*) '++++++++++++'
				        goto 2000
			        ENDIF
			        PB_e%PO(num_elec+1)%Vx = V_x
			        PB_e%PO(num_elec+1)%Vy = V_y
			        PB_e%PO(num_elec+1)%Vz = V_z

                    PB_e%PO(num_elec+1)%ax = 0.
			        PB_e%PO(num_elec+1)%ay = 0.
			        PB_e%PO(num_elec+1)%az = 0.
            
			        !part(num+1,7) = 1

			        ns(1) = ns(1) + 1
                    PB_e%NPar=PB_e%NPar+1
		        ENDDO
            ENDIF
        ENDDO





        WRITE(6,*) '# Quasi: to inject on one side',ns(1)-ne_old

        ntot = ns(1) + ns(2)

        print*,'after quasi:'
        print*,'ns(1)=',ns(1)
        print*,'ns(2)=',ns(2)

    
    End subroutine   Quasi_cdp
    
    
    SUBROUTINE	INDEXM_cdp(PB_e,PB_i,insp, IC, nnx, nny)


        USE Domain_2D

        IMPLICIT NONE
        Class(ParticleBundle), intent(inout) :: PB_e,PB_i
        
        INTEGER								::	insp, nsp, nnx, nny
        INTEGER, DIMENSION(2,nnx * nny)		::	IC
        INTEGER								::	i, j, k, MC, ipret, jpret, num, ip,isp
        REAL(8)								::	x, y
        INTEGER								::	npt

        !IC = 0
        !
        !j = 1
        !
        !num = SIZE(part, 1)
        !npt = ntot
        !
        !DO i=1, ntot
        !	IF(part(i,7) .NE. 0) THEN
        !		x = part(i,1)
        !		y = part(i,2)
        !		nsp = part(i,7)
        !		IF(nsp == insp) THEN
        !			ipret = INT((x - 0) / hx(1)) + 1
        !			jpret = INT((y - 0) / hx(2)) + 1
        !			!IF(ipret > nnx .OR. ipret <= 0 .OR. jpret > nny .OR. jpret <= 0) THEN
        !			IF(ipret > nnx .OR. ipret <= 0 .OR. jpret > nny .OR. jpret <= 0 .OR. x <0 .OR. y<0) THEN
        !				DO ip=1,SIZE(part,2)
        !					part(i, ip) = part(npt, ip)
        !				END DO
        !				ns(nsp) = ns(nsp) - 1
        !				npt=npt-1
        !				goto 100
        !			ENDIF
        !			MC = (ipret - 1) * (ny - 1) + jpret
        !			IC(2,MC) = IC(2,MC) + 1
        !		ENDIF
        !100	ENDIF
        !ENDDO
        !
        !ntot = npt




        IC = 0

        j = 1

        num = 10000000
        !npt = ntot


         IF(1 == insp) THEN
            DO i = 1, PB_e%Npar
	        !IF(part(i,7) .NE. 0) THEN
                !IF(isp .NE. 0) THEN
        
		        x = PB_e%PO(i)%X
		        y = PB_e%PO(i)%Y
		        !nsp = isp+1
		        !IF(1 == insp) THEN
			        ipret = INT((x - 0) / hx(1)) + 1
			        jpret = INT((y - 0) / hx(2)) + 1
			        !IF(ipret > nnx .OR. ipret <= 0 .OR. jpret > nny .OR. jpret <= 0) THEN
			        IF(ipret > nnx .OR. ipret <= 0 .OR. jpret > nny .OR. jpret <= 0 .OR. x <0 .OR. y<0) THEN
				        !DO ip=1,SIZE(part,2)
				        !	part(i, ip) = part(npt, ip)
				        !END DO
				        !ns(nsp) = ns(nsp) - 1
				        !npt=npt-1
                         Call PB_e%DelOne(i) 
			        else
			        MC = (ipret - 1) * (ny - 1) + jpret
			        IC(2,MC) = IC(2,MC) + 1
                    ENDIF
	        !100	ENDIF
	        !ENDIF
           
            ENDDO
         endif
 
    
         IF(2 == insp) THEN   
             DO i = 1, PB_i%Npar
	        !IF(part(i,7) .NE. 0) THEN
                !IF(isp .NE. 0) THEN
        
		        x = PB_i%PO(i)%X
		        y = PB_i%PO(i)%Y
		
		
			        ipret = INT((x - 0) / hx(1)) + 1
			        jpret = INT((y - 0) / hx(2)) + 1
			        !IF(ipret > nnx .OR. ipret <= 0 .OR. jpret > nny .OR. jpret <= 0) THEN
			        IF(ipret > nnx .OR. ipret <= 0 .OR. jpret > nny .OR. jpret <= 0 .OR. x <0 .OR. y<0) THEN
				        !DO ip=1,SIZE(part,2)
				        !	part(i, ip) = part(npt, ip)
				        !END DO
				        !ns(nsp) = ns(nsp) - 1
				        !npt=npt-1
                         Call PB_i%DelOne(i) 
				
			        else
			        MC = (ipret - 1) * (ny - 1) + jpret
			        IC(2,MC) = IC(2,MC) + 1
                    ENDIF
		
	        !ENDIF
             ENDDO
         endif
 

        !ntot = npt



    END SUBROUTINE INDEXM_cdp
    
    
    
    
    
    
    
    
     
End Module ModuleMCCInterface