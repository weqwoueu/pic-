!===============================================================================
!  Module:  ModuleGlobalDiagnostics
!  Purpose: Global energy/particle conservation diagnostics.
!           Outputs global_diagnostics.csv per the verification protocol.
!
!           Uses ModuleDiagCounters for shared loss/injection accumulators
!           (updated directly by MCCInterface, read here).  No circular deps.
!
!  Output columns:
!    t, Ne_domain, Ni_domain, Ne_lost, Ni_lost, Ne_injected, Ni_injected,
!    Ke, Ki, Efield, Einjected, Elost, Ebalance_error
!===============================================================================

Module ModuleGlobalDiagnostics
    Use ModuleMCCInterface, Only: ParticleGlobal, ControlFlowGlobal
    Use ModuleDiagCounters, Only: E_lost_cum, E_inj_cum, &
                                  Ne_lost_cum, Ni_lost_cum, &
                                  Ne_inj_cum, Ni_inj_cum
    Use Field_2D,            Only: efx, efy
    Use IFE_Data,            Only: Field_Size
    Use Constant_Variable_2D,Only: EPSILON0, Efield_ref
    Use Cell_Volume_Data,    Only: Cell_Volume_zwz
    Implicit None

    ! ---- initial values for balance check ----
    Real(8), Save, Private :: E_dom_0     = 0.0d0
    Real(8), Save, Private :: Ne_0        = 0.0d0
    Real(8), Save, Private :: Ni_0        = 0.0d0

    ! ---- file unit ----
    Integer, Save, Private :: diag_unit = 545

    ! ---- public interface ----
    Public :: InitGlobalDiagnostics
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

        ! ---- record initial E_dom, Ne, Ni ----
        E_dom_0 = 0.0d0
        Ne_0    = 0.0d0
        Ni_0    = 0.0d0
        Do isp = 0, ControlFlowGlobal%Ns
            mass_j  = ParticleGlobal(isp)%Mass
            vfactor = ParticleGlobal(isp)%VFactor
            weight  = ParticleGlobal(isp)%Weight
            Do i = 1, ParticleGlobal(isp)%Npar
                ek = ParticleGlobal(isp)%PO(i)%Energy(mass_j, vfactor)
                If (isp == 0) Then
                    Ne_0 = Ne_0 + weight
                Else
                    Ni_0 = Ni_0 + weight
                End If
                E_dom_0 = E_dom_0 + ek * weight
            End Do
        End Do
        E_dom_0 = E_dom_0 + ComputeFieldEnergy()

        ! Reset shared counters
        E_lost_cum  = 0.0d0
        E_inj_cum   = 0.0d0
        Ne_lost_cum = 0.0d0
        Ni_lost_cum = 0.0d0
        Ne_inj_cum  = 0.0d0
        Ni_inj_cum  = 0.0d0

        ! Write initial row (t=0)
        Call AppendRow(0.0d0, Ne_0, Ni_0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, E_dom_0)
    End Subroutine InitGlobalDiagnostics

    !---------------------------------------------------------------------------
    ! Append a row to the CSV file with the given pre-computed values.
    ! (Called by Init for t=0, by WriteGlobalDiagnostics for each step.)
    !---------------------------------------------------------------------------
    Subroutine AppendRow(tval, Ne_dom, Ni_dom, Ke, Ki, Efield, Ebalance)
        Implicit None
        Real(8), Intent(In) :: tval, Ne_dom, Ni_dom, Ke, Ki, Efield, Ebalance
        Open(diag_unit, File='./OUTPUT/global_diagnostics.csv', &
             Action='Write', Position='Append', Status='Old')
        Write(diag_unit, '(E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8,",", &
            &E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8,",",E15.8)') &
            tval, Ne_dom, Ni_dom, Ne_lost_cum, Ni_lost_cum, &
            Ne_inj_cum, Ni_inj_cum, Ke, Ki, Efield, E_inj_cum, E_lost_cum, Ebalance
        Close(diag_unit)
    End Subroutine AppendRow

    !---------------------------------------------------------------------------
    ! Compute current electrostatic field energy in Joules.
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
    ! Write a diagnostic row for the current time step.
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

        Call AppendRow(xt, Ne_dom, Ni_dom, Ke, Ki, Efield, Ebalance)
    End Subroutine WriteGlobalDiagnostics

End Module ModuleGlobalDiagnostics
