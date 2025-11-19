Module DiagnosticsTestParticle
    Use ModuleGrid
    Use ModuleField
    Use ModuleParticleBundle
    Use ModuleControlFlow
    Use ModuleParticleBoundary
    Implicit none
    contains  
  
        SUBROUTINE Initial_Diag_Particle_Test(PB)
            implicit none 
            Type(ParticleBundle),intent(inout) :: PB
                real(8) :: Initial_Position
                Initial_Position = 0.05         !Real first relative position whose value should be less than 1
                PB%NPar     = 1
                PB%PO(1)%X  = ZLength * Initial_Position / PB%Xfactor
                PB%PO(1)%Vx = 0.
                PB%PO(1)%Vy = 0.
                PB%PO(1)%Vz = 0.
                PB%PO(1)%Ax = 0.
                PB%PO(1)%Ay = 0.
                PB%PO(1)%Ax = 0.
                PB%PO(1)%Ay = 0.
                PB%PO(1)%Az = 0.
            return
        end subroutine Initial_Diag_Particle_Test
        
        SUBROUTINE Initial_Diag_Particle_Test2(PB)
            implicit none 
            Type(ParticleBundle),intent(inout) :: PB
                real(8) :: Reset_Position
                Reset_Position = 0.5         !Real Reset relative position that should be less than 1
                PB%NPar     = 1
                PB%PO(1)%X  = ZLength  * Reset_Position / PB%Xfactor
                PB%PO(1)%Vx = 0.
                PB%PO(1)%Vy = 0.
                PB%PO(1)%Vz = 0.
                PB%PO(1)%Ax = 0.
                PB%PO(1)%Ay = 0.
                PB%PO(1)%Ax = 0.
                PB%PO(1)%Ay = 0.
                PB%PO(1)%Az = 0.
            return
        end subroutine Initial_Diag_Particle_Test2

    
        subroutine Diag_Particle_Test_In_Field(PB,PBDO,CF,Ni,Nj,FG)
            Implicit none
            Type(Field),intent(in) :: FG
            Type(ControlFlow),intent(in):: CF
            Type(ParticleBundle),intent(in) :: PB
            Type(ParticleBoundaryOne),intent(in) :: PBDO
            Type(ParticleBundle),save :: TestPB,TempPB
            Type(ParticleOne) :: TempPO
            Integer(4) :: Ni,Nj
            integer(8) :: NTime,Control_Parameter
            Real(8) :: Energy,Real_Time

            NTime = Ni + (Nj - 1) * CF%Period
            Real_Time = NTime * PB%Dt
            Control_Parameter = MOD(NTime,PB%NParNormal)
            if(Ni==1) then 
                if(NTime == 1)then
                write(*,*) "Diag_Particle_Test_In_Field Begin:"
                end if
                write(*,*)  "Period ",Nj,"  total step :",NTime
            end if
            IF (NTime == 1) then  
                TempPB = PB
                TestPB = PB
                call Initial_Diag_Particle_Test(TempPB)
                call Initial_Diag_Particle_Test(TestPB)
                TestPB = TempPB
            end if
            TempPB%PO(1)%X = TestPB%PO(Control_Parameter)%X
            TempPO = TempPB%PO(1)
            If (TempPO%X<=PBDO%XMin.or.TempPO%X>=PBDO%XMax) Then
                call Initial_Diag_Particle_Test2(TempPB)
                TempPO = TempPB%PO(1) 
                Energy = TempPO%Energy(TempPB%Mass,TempPB%Vfactor)
                TempPO%Az = Energy / JtoeV       
                TempPO%Ay = Real_Time
                Call TestPB%Addone(TempPO)
                !Stop
            ELSE
                call MoveElectromagneticParticleBundle(TempPB,FG)
                TempPO = TempPB%PO(1)
                Energy = TempPO%Energy(TempPB%Mass,TempPB%Vfactor)
                TempPO%Az = Energy / JtoeV       
                TempPO%Ay = Real_Time
                Call TestPB%Addone(TempPO)
                !Write(*,*) TestPB%Npar
            END IF
            if((Control_Parameter==0) .OR. (NTime == CF%NDiagShort * CF%Period)) then 
                call Dump_Diag_Particle_Test(TestPB)
            end if  
            return
        end subroutine Diag_Particle_Test_In_Field 
        
        subroutine Dump_Diag_Particle_Test(PB)
            implicit none
            Type(ParticleBundle),intent(inout) :: PB
            Type(ParticleOne) :: TempPO
            character(len=40) :: filename = "Diag Particle In Field.dat"
            Integer(4) :: i,j
            open(Unit = 100,file=filename,Position='APPEND')
            !open(Unit = 100,file=filename)
            do i=1,PB%NPar
                write(100,FMt="(7(es21.14,1x))") PB%PO(i)%Ay,PB%PO(i)%X*PB%Xfactor,&
                                                 PB%PO(i)%Vx*PB%Vfactor,PB%PO(i)%Vy*PB%Vfactor,PB%PO(i)%Vz*PB%Vfactor,&
                                                 PB%PO(i)%Ax,PB%PO(i)%Az
            end do
            close(100)
            TempPO = PB%PO(PB%NPar)
            call Initial_Diag_Particle_Test(PB)
            PB%PO(1) = TempPO
            !Control_Parameter = 1
            write(*,*) "write data successful"
            return
        end subroutine Dump_Diag_Particle_Test
  End Module DiagnosticsTestParticle