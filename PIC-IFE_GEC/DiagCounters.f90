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
End Module ModuleDiagCounters
