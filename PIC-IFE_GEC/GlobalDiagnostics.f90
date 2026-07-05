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
                                  Ne_inj_cum, Ni_inj_cum, &
                                  ResetDiagCounters
    Use Field_2D,            Only: efx, efy
    Use IFE_Data,            Only: Field_Size
    Use Constant_Variable_2D,Only: EPSILON0, Efield_ref, L_ref
    Use Cell_Volume_Data,    Only: Cell_Volume_zwz
    Use Domain_2D,           Only: delta_global
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
    Subroutine InitGlobalDiagnostics(t0)
        Implicit None
        Real(8), Intent(In), Optional :: t0
        Real(8)  :: t_init
        Real(8)  :: Ke0, Ki0, Efield0

        ! Open CSV file (overwrite)
        Open(diag_unit, File='./OUTPUT/global_diagnostics.csv', &
             Action='Write', Status='Replace')
        Write(diag_unit, '(A)') &
            't,Ne_domain,Ni_domain,Ne_lost,Ni_lost,Ne_injected,Ni_injected,' // &
            'Ke,Ki,Efield,Einjected,Elost,Ebalance_error'
        Close(diag_unit)

        ! ---- record initial E_dom, Ne, Ni ----
        t_init = 0.0d0
        If (Present(t0)) t_init = t0

        Call ComputeDomainTotals(Ke0, Ki0, Ne_0, Ni_0)
        Efield0 = ComputeFieldEnergy()
        E_dom_0 = Ke0 + Ki0 + Efield0

        ! Reset shared counters
        Call ResetDiagCounters()

        ! Write the actual initial state after particles and fields are ready.
        Call AppendRow(t_init, Ne_0, Ni_0, Ke0, Ki0, Efield0, 0.0d0)
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
    Function ParticleWeight(isp, i) Result(weight)
        Implicit None
        Integer, Intent(In) :: isp, i
        Real(8) :: weight

        If (ParticleGlobal(isp)%UnequalWeightFlag) Then
            weight = ParticleGlobal(isp)%PO(i)%WQ
        Else
            weight = ParticleGlobal(isp)%Weight
        End If
    End Function ParticleWeight

    !---------------------------------------------------------------------------
    Subroutine ComputeDomainTotals(Ke, Ki, Ne_dom, Ni_dom)
        Implicit None
        Real(8), Intent(Out) :: Ke, Ki, Ne_dom, Ni_dom
        Integer  :: isp, i
        Real(8)  :: mass_j, vfactor, weight, ek

        Ke = 0.0d0
        Ki = 0.0d0
        Ne_dom = 0.0d0
        Ni_dom = 0.0d0

        Do isp = 0, ControlFlowGlobal%Ns
            mass_j  = ParticleGlobal(isp)%Mass
            vfactor = ParticleGlobal(isp)%VFactor
            Do i = 1, ParticleGlobal(isp)%Npar
                weight = ParticleWeight(isp, i)
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
    End Subroutine ComputeDomainTotals

    !---------------------------------------------------------------------------
    ! Compute current electrostatic field energy in Joules.
    !---------------------------------------------------------------------------
    Function ComputeFieldEnergy() Result(efield_j)
        Implicit None
        Real(8) :: efield_j
        Integer :: i
        Real(8) :: ex_phys, ey_phys, dV, geom_scale

        efield_j = 0.0d0
        ! Cell_Volume_zwz is stored in normalized coordinates.
        geom_scale = L_ref * L_ref
        If (delta_global == 1) geom_scale = L_ref * L_ref * L_ref

        Do i = 1, Field_Size
            ex_phys = efx(i, 1) * Efield_ref
            ey_phys = efy(i, 1) * Efield_ref
            dV      = Cell_Volume_zwz(2, i) * geom_scale
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
        Real(8)  :: Ke, Ki, Efield, E_dom_t
        Real(8)  :: Ne_dom, Ni_dom
        Real(8)  :: Ebalance

        Call ComputeDomainTotals(Ke, Ki, Ne_dom, Ni_dom)

        Efield  = ComputeFieldEnergy()
        E_dom_t = Ke + Ki + Efield
        If (Abs(E_dom_0) > 1.0d-300) Then
            Ebalance = (E_dom_t + E_lost_cum - E_inj_cum - E_dom_0) / E_dom_0
        Else
            Ebalance = 0.0d0
        End If

        Call AppendRow(xt, Ne_dom, Ni_dom, Ke, Ki, Efield, Ebalance)
    End Subroutine WriteGlobalDiagnostics

End Module ModuleGlobalDiagnostics
