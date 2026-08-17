!===============================================================================
!  Module:  ModuleDiagCounters
!  Purpose: Shared loss/injection counters, used by both MCCInterface and
!           GlobalDiagnostics.  NO circular module dependencies.
!===============================================================================

Module ModuleDiagCounters
    Implicit None
    Real(8), Save, Public :: E_lost_cum  = 0.0d0
    Real(8), Save, Public :: E_inj_cum   = 0.0d0
    Real(8), Save, Public :: Ne_lost_cum = 0.0d0
    Real(8), Save, Public :: Ni_lost_cum = 0.0d0
    Real(8), Save, Public :: Ne_inj_cum  = 0.0d0
    Real(8), Save, Public :: Ni_inj_cum  = 0.0d0

    Public :: ResetDiagCounters
    Public :: AddDiagLost
    Public :: AddDiagInjected
    Public :: AddDiagThermalExchange

Contains

    Subroutine ResetDiagCounters()
        Implicit None
        E_lost_cum  = 0.0d0
        E_inj_cum   = 0.0d0
        Ne_lost_cum = 0.0d0
        Ni_lost_cum = 0.0d0
        Ne_inj_cum  = 0.0d0
        Ni_inj_cum  = 0.0d0
    End Subroutine ResetDiagCounters

    Subroutine AddDiagLost(isp, weight, energy)
        Implicit None
        Integer(4), Intent(In) :: isp
        Real(8), Intent(In)    :: weight, energy

        E_lost_cum = E_lost_cum + energy * weight
        If (isp == 0) Then
            Ne_lost_cum = Ne_lost_cum + weight
        Else
            Ni_lost_cum = Ni_lost_cum + weight
        End If
    End Subroutine AddDiagLost

    Subroutine AddDiagInjected(isp, weight, energy)
        Implicit None
        Integer(4), Intent(In) :: isp
        Real(8), Intent(In)    :: weight, energy

        E_inj_cum = E_inj_cum + energy * weight
        If (isp == 0) Then
            Ne_inj_cum = Ne_inj_cum + weight
        Else
            Ni_inj_cum = Ni_inj_cum + weight
        End If
    End Subroutine AddDiagInjected

    Subroutine AddDiagThermalExchange(isp, weight, energy_lost, energy_injected)
        Implicit None
        Integer(4), Intent(In) :: isp
        Real(8), Intent(In) :: weight, energy_lost, energy_injected

        E_lost_cum = E_lost_cum + energy_lost * weight
        E_inj_cum  = E_inj_cum  + energy_injected * weight
        If (isp == 0) Then
            Ne_lost_cum = Ne_lost_cum + weight
            Ne_inj_cum  = Ne_inj_cum  + weight
        Else
            Ni_lost_cum = Ni_lost_cum + weight
            Ni_inj_cum  = Ni_inj_cum  + weight
        End If
    End Subroutine AddDiagThermalExchange
End Module ModuleDiagCounters
