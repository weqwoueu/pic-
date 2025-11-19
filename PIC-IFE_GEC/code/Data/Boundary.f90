Module ModuleBoundary
    Use ModuleShape
    Implicit none

    
    Type Boundary
        Integer(4) :: Shape=1 ! 1.Line  2.Circle
	    Class(GeneralLine),Allocatable :: BoundaryLine
        Logical :: SEEFlag = .True.
    contains
        procedure :: Alloc=>AllocateBoundaryLine
        procedure :: AdjustPar=>AdjustParticleBoundary
    End Type Boundary
    
    Type,extends(Boundary) :: AxisBoundary
        
    End Type AxisBoundary
    
    Type,extends(Boundary) :: ReflexBoundary
        
    End Type ReflexBoundary
    

    Type,extends(Boundary) :: AbortionBoundary
        
    End Type AbortionBoundary
    
    
    Type,extends(AbortionBoundary) :: ConductorBoundary
        
    End Type ConductorBoundary
    
    
    Type,extends(AbortionBoundary) :: DielectricBoundary
        
    End Type DielectricBoundary
    
    Type BoundaryAbstract
        Integer(4) :: BoundaryType=3 ! 1.Axis  2.Reflex  3.Abortion  31.Conductor  32.Dielectric
        Class(Boundary),Allocatable :: B
    contains    
        procedure :: Alloc=>AllocateBoundary
    End Type BoundaryAbstract
    
    Contains
    
    Subroutine AllocateBoundaryLine(B)
        Class(Boundary),intent(inout) :: B
        Select Case (B%Shape)
        Case(1)
            Allocate(Line::B%BoundaryLine)
        Case(2)
            Allocate(Circle::B%BoundaryLine)
        End Select
        return
    End Subroutine AllocateBoundaryLine
    
    
    Subroutine AdjustParticleBoundary(B,PO)
        Use ModuleParticleOne
        Class(Boundary),intent(in) :: B
        Class(ParticleOne),intent(inout) :: PO
        Select Type (B)
        Type is (AxisBoundary)
            
        Type is (ReflexBoundary)
            
        Type is (AbortionBoundary)
            
        End Select
        return
    End Subroutine AdjustParticleBoundary
    
    
    Subroutine AllocateBoundary(BA,S)
        Class(BoundaryAbstract),intent(inout) :: BA
        Integer(4),intent(in),optional :: S
        Select Case (BA%BoundaryType)
        Case(1)
            Allocate(AxisBoundary::BA%B)
        Case(2)
            Allocate(ReflexBoundary::BA%B)
        Case(3)
            Allocate(AbortionBoundary::BA%B)
        Case(31)
            Allocate(ConductorBoundary::BA%B)
        Case(32)
            Allocate(DielectricBoundary::BA%B)
        End Select
        If (present(S)) Then
            BA%B%Shape = S
        Endif
        Call BA%B%Alloc
    End Subroutine AllocateBoundary
    
End Module ModuleBoundary