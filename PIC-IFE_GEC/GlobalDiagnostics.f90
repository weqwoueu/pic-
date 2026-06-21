!===============================================================================
!  Module:  ModuleGlobalDiagnostics
!  Purpose: Global energy/particle conservation diagnostics.
!           Outputs global_diagnostics.csv per the verification protocol.
!
!  Output columns:
!    t, Ne_domain, Ni_domain, Ne_lost, Ni_lost, Ne_injected, Ni_injected,
!    Ke, Ki, Efield, Einjected, Elost, Ebalance_error
!===============================================================================

Module ModuleGlobalDiagnostics
    Use ModuleMCCInterface, Only: ParticleGlobal, ControlFlowGlobal
    Use Field_2D,            Only: efx, efy
    Use IFE_Data,            Only: Field_Size
    Use Domain_2D,           Only: hx
    Use Constant_Variable_2D,Only: EPSILON0, Efield_ref, Phi_ref, L_ref
    Use TimeControl,         Only: dt, nt
    Use Cell_Volume_Data,    Only: Cell_Volume_zwz
    Implicit None

    ! ---- private, save, cumulative trackers ----
    Real(8), Save, Private :: E_lost_cum  = 0.0d0    ! total escaped kinetic energy (J)
    Real(8), Save, Private :: E_inj_cum   = 0.0d0    ! total injected kinetic energy (J)

    Real(8), Save, Private :: Ne_lost_cum = 0.0d0    ! weighted electron count lost
    Real(8), Save, Private :: Ni_lost_cum = 0.0d0    ! weighted ion count lost
    Real(8), Save, Private :: Ne_inj_cum  = 0.0d0    ! weighted electron count injected
    Real(8), Save, Private :: Ni_inj_cum  = 0.0d0    ! weighted ion count injected

    ! ---- initial values for balance check ----
    Real(8), Save, Private :: E_dom_0     = 0.0d0
    Real(8), Save, Private :: Ne_0        = 0.0d0
    Real(8), Save, Private :: Ni_0        = 0.0d0

    ! ---- file unit ----
    Integer,  Save, Private :: diag_unit   = 545
    Logical,  Save, Private :: is_open     = .False.

    ! ---- public interface ----
    Public :: InitGlobalDiagnostics
    Public :: RecordParticleLost
    Public :: RecordParticleInjected
    Public :: WriteGlobalDiagnostics

Contains

    !---------------------------------------------------------------------------
    Subroutine InitGlobalDiagnostics()
        Implicit None
        Integer  :: isp, i
        Real(8)  :: ek, mass_j, vfactor, weight

        ! Open CSV file (overwrite)
        Open(diag_unit, File='./OUTPUT/global_diagnostics.csv', &
             Action='Write', Status='Replace')
        Write(diag_unit, '(A)') &
            't,Ne_domain,Ni_domain,Ne_lost,Ni_lost,Ne_injected,Ni_injected,' // &
            'Ke,Ki,Efield,Einjected,Elost,Ebalance_error'
        Close(diag_unit)
        is_open = .True.

        ! ---- record initial E_dom, Ne, Ni ----
        ! Kinetic energy of all particles
        E_dom_0 = 0.0d0
        Ne_0    = 0.0d0
        Ni_0    = 0.0d0
        Do isp = 0, ControlFlowGlobal%Ns
            mass_j  = ParticleGlobal(isp)%Mass
            vfactor = ParticleGlobal(isp)%VFactor
            weight  = ParticleGlobal(isp)%Weight
            Do i = 1, ParticleGlobal(isp)%Npar
                ek = ParticleGlobal(isp)%PO(i)%Energy(mass_j, vfactor)  ! Joules
                If (isp == 0) Then
                    Ne_0 = Ne_0 + weight
                Else
                    Ni_0 = Ni_0 + weight
                End If
                E_dom_0 = E_dom_0 + ek * weight
            End Do
        End Do
        ! Add initial field energy
        E_dom_0 = E_dom_0 + ComputeFieldEnergy()

        ! Reset cumulative counters
        E_lost_cum  = 0.0d0
        E_inj_cum   = 0.0d0
        Ne_lost_cum = 0.0d0
        Ni_lost_cum = 0.0d0
        Ne_inj_cum  = 0.0d0
        Ni_inj_cum  = 0.0d0

        ! Write initial row
        Call AppendRow(0.0d0, 0)
    End Subroutine InitGlobalDiagnostics

    !---------------------------------------------------------------------------
    ! Called when a particle is about to be deleted (escaped domain)
    ! ispe : species index (0 = electron, >0 = ion)
    !---------------------------------------------------------------------------
    Subroutine RecordParticleLost(ispe, mass_kg, vfactor, vx, vy, vz, weight)
        Implicit None
        Integer, Intent(In) :: ispe
        Real(8), Intent(In) :: mass_kg, vfactor, vx, vy, vz, weight
        Real(8) :: ek

        ek = 0.5d0 * mass_kg * (vx*vx + vy*vy + vz*vz) * vfactor * vfactor
        E_lost_cum = E_lost_cum + ek * weight
        If (ispe == 0) Then
            Ne_lost_cum = Ne_lost_cum + weight
        Else
            Ni_lost_cum = Ni_lost_cum + weight
        End If
    End Subroutine RecordParticleLost

    !---------------------------------------------------------------------------
    ! Called when a particle is injected at a boundary
    !---------------------------------------------------------------------------
    Subroutine RecordParticleInjected(ispe, mass_kg, vfactor, vx, vy, vz, weight)
        Implicit None
        Integer, Intent(In) :: ispe
        Real(8), Intent(In) :: mass_kg, vfactor, vx, vy, vz, weight
        Real(8) :: ek

        ek = 0.5d0 * mass_kg * (vx*vx + vy*vy + vz*vz) * vfactor * vfactor
        E_inj_cum = E_inj_cum + ek * weight
        If (ispe == 0) Then
            Ne_inj_cum = Ne_inj_cum + weight
        Else
            Ni_inj_cum = Ni_inj_cum + weight
        End If
    End Subroutine RecordParticleInjected

    !---------------------------------------------------------------------------
    ! Compute current electrostatic field energy in Joules.
    ! E_field = 0.5 * epsilon0 * integral (E_x^2 + E_y^2) dV
    !   efx, efy  : normalized E-field
    !   Efield_ref = Phi_ref / L_ref   (physical E per normalized unit)
    !   Cell_Volume_zwz(2,i) : physical cell volume in m^3 (already includes L_ref)
    !---------------------------------------------------------------------------
    Function ComputeFieldEnergy() Result(efield_j)
        Implicit None
        Real(8) :: efield_j
        Integer :: i
        Real(8) :: ex_phys, ey_phys, dV

        efield_j = 0.0d0
        Do i = 1, Field_Size
            ex_phys = efx(i, 1) * Efield_ref
            ey_phys = efy(i, 1) * Efield_ref
            dV      = Cell_Volume_zwz(2, i)
            efield_j = efield_j + 0.5d0 * EPSILON0 * (ex_phys*ex_phys + ey_phys*ey_phys) * dV
        End Do
    End Function ComputeFieldEnergy

    !---------------------------------------------------------------------------
    ! Compute current domain energy and particle counts, then append CSV row.
    !---------------------------------------------------------------------------
    Subroutine WriteGlobalDiagnostics(it, xt)
        Implicit None
        Integer, Intent(In) :: it
        Real(8), Intent(In) :: xt
        Integer  :: isp, i
        Real(8)  :: Ke, Ki, Efield, E_dom_t
        Real(8)  :: Ne_dom, Ni_dom
        Real(8)  :: mass_j, vfactor, weight, ek
        Real(8)  :: Ebalance

        Ke = 0.0d0
        Ki = 0.0d0
        Ne_dom = 0.0d0
        Ni_dom = 0.0d0

        Do isp = 0, ControlFlowGlobal%Ns
            mass_j  = ParticleGlobal(isp)%Mass
            vfactor = ParticleGlobal(isp)%VFactor
            weight  = ParticleGlobal(isp)%Weight
            Do i = 1, ParticleGlobal(isp)%Npar
                ek = ParticleGlobal(isp)%PO(i)%Energy(mass_j, vfactor)
                If (isp == 0) Then
                    Ke = Ke + ek * weight
                    Ne_dom = Ne_dom + weight
                Else
                    Ki = Ki + ek * weight
                    Ni_dom = Ni_dom + weight
                End If
            End Do
        End Do

        Efield  = ComputeFieldEnergy()
        E_dom_t = Ke + Ki + Efield

        Ebalance = (E_dom_t + E_lost_cum - E_inj_cum - E_dom_0) / E_dom_0

        Call AppendRow(xt, it)

        Contains
            Subroutine AppendRow(tval, itval)
                Real(8), Intent(In) :: tval
                Integer, Intent(In) :: itval
                Open(diag_unit, File='./OUTPUT/global_diagnostics.csv', &
                     Action='Write', Position='Append', Status='Old')
                Write(diag_unit, '(E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8,",", &
                    &E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8)') &
                    tval, Ne_dom, Ni_dom, Ne_lost_cum, Ni_lost_cum, &
                    Ne_inj_cum, Ni_inj_cum, Ke, Ki, Efield, E_inj_cum, E_lost_cum, Ebalance
                Close(diag_unit)
            End Subroutine AppendRow
    End Subroutine WriteGlobalDiagnostics

End Module ModuleGlobalDiagnostics
