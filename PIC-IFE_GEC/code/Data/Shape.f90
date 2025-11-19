Module ModuleShape
    !Use Constants, Only:Pi
    Implicit none
    
    Real(8), PARAMETER :: Unit = 1.
    
    Type Point
        Real(8) :: X=0, Y=0, Z=0
    contains
        procedure :: Init=>InitializePoint
        procedure :: Copy=>CopyPoint
        procedure :: DistanceP=>DistanceFromPoint
        procedure :: InLine
    End Type Point
    
    Type(Point) :: Origin !> origin point, default: (x,y,z) = (0,0,0)
        
    Type,extends(Point) :: Vector
    contains
        procedure :: Modulus=>ModulusVector
        procedure :: Scale=>ScaleVector
        procedure :: PerpenV=>GetPerpendicularDirectionVectorFormVector
    End Type Vector
    
    Type GeneralLine
        
    contains
        procedure :: Init=>InitializeGeneralLine
        procedure :: GLInters => GeneralLineIntersection
    End Type GeneralLine
    
    Type,extends(GeneralLine) :: Line !> default direction: P1 -> P2
        Type(Point) :: P1
        Type(Point) :: P2
    contains
        procedure :: Init=>InitializeLine
        procedure :: Copy=>CopyLine
        procedure :: Length=>LengthLine
        procedure :: LInters=>LineIntersection
        procedure :: LCInters=>LineCircleIntersection
        procedure :: DirectionV=>GetDirectionVector
        procedure :: PerpenV=>GetPerpendicularDirectionVectorFromLine
    End Type Line
    
    
    Type,extends(GeneralLine) :: Circle
        Type(Point) :: Center
        Real(8) :: R
    contains
        procedure :: Perimeter=>PerimeterCircle
    End Type Circle
    
    !Type Shape
    !    
    !End Type Shape
    !
    !Type,extends(Shape) :: Rectangle
    !    Type(DLine) :: Edges(4)
    !End Type Rectangle
    
    Contains

    Logical Function Eq(a,b)
        Real(8),intent(in) :: a,b
        Eq = Abs(a-b) < 1E-6
    End Function Eq
    
    !>--------------------------------------------------------------------------
    
    Subroutine InitializePoint(P,X,Y,Z)
        Class(Point),intent(inout) :: P
        Real(8),intent(in),optional :: X,Y,Z
        If (present(X)) P%X = X
        If (present(Y)) P%Y = Y
        If (present(Z)) P%Z = Z
        return
    End Subroutine InitializePoint
    
    Subroutine CopyPoint(PD,PC)
        Class(Point),intent(inout) :: PD
        Class(Point),intent(in) :: PC
        Select Type (PD)
        Type is (Point)
            Select Type (PC)
            Type is (Point)
                PD = PC
            End Select
        End Select
    End Subroutine CopyPoint
    
    
    Real(8) Function DistanceFromPoint(P,PTarget) 
        Implicit none
        Class(Point),intent(in) :: P,PTarget
        DistanceFromPoint = Sqrt( Abs(P%X-PTarget%X)**2 + Abs(P%Y-PTarget%Y)**2 )
        return 
    End Function DistanceFromPoint
    
    
    Logical Function InLine(P,L)
        Implicit none
        Class(Point),intent(in) :: P
        Class(GeneralLine),intent(in) :: L
        InLine = .False.
        Select Type (L)
        Type is (Line)
            Associate(X1=>L%P1%X, X2=>L%P2%X, Y1=>L%P1%Y, Y2=>L%P2%Y, X=>P%X, Y=>P%Y)
                If ( Eq( (Y1 - Y2) * (X - X1), (X1 - X2) * (Y - Y1) ) ) Then    ! in line
                    If ((X2-X)*(X1-X)<=0 .and. (Y1-Y)*(Y2-Y)<0) Then            ! in line segment
                        InLine = .True.
                    End If
                End If
            End Associate
        Type is (Circle)
            If ( Eq(P%DistanceP(L%Center),L%R) ) Then
                InLine = .True.
            Endif
        End Select
    End Function InLine 
    
    !> ----------------------------------------------------------------------------------------------
    
    Real(8) Function ModulusVector(V)
        Class(Vector),intent(in) :: V
        ModulusVector = V%DistanceP(Origin)
    End Function ModulusVector
    
    Subroutine ScaleVector(V,Factor)
        Class(Vector),intent(inout) :: V
        Real(8),intent(in) :: Factor
        V%X = V%X * Factor
        V%Y = V%Y * Factor
        V%Z = V%Z * Factor
    End Subroutine ScaleVector
    
    !> Side: -1 -- left || 1 -- right || default -- right
    Type(Vector) Function GetPerpendicularDirectionVectorFormVector(V,Side)
        Class(Vector),intent(in) :: V
        Integer(4),intent(in),optional :: Side 
        Associate(Output=>GetPerpendicularDirectionVectorFormVector)
            Call Output%init(-V%Y, V%X)
            If (.Not.Eq(Output%Modulus(), Unit)) Call Output%Scale(1/Output%Modulus())
            If (present(Side)) Call Output%Scale(Real(Side,8))
        End Associate
    End Function GetPerpendicularDirectionVectorFormVector
    !> ----------------------------------------------------------------------------------------------
    
    ! when the GL is circle x2 represent R
    Subroutine InitializeGeneralLine(L,X1,Y1,X2,Y2)
        Class(GeneralLine),intent(inout) :: L
        Real(8),intent(in) :: X1,Y1,X2,Y2
        Select Type (L)
        Type is (Line)
            Call L%init(X1,Y1,X2,Y2)
        Type is (Circle)
            !TODO
        End Select
        return
    End Subroutine InitializeGeneralLine
    
    Integer(4) Function GeneralLineIntersection(GL1,GL2,InterPoint1,InterPoint2)
        Class(GeneralLine),intent(in) :: GL1, GL2
        Class(Point),intent(out),optional :: InterPoint1,InterPoint2
        Select Type (GL1)
        Type is (Line)
            Select Type (GL2)
            Type is (Line)
                GeneralLineIntersection = GL1%LInters(GL2,InterPoint1)
            Type is (Circle)
                GeneralLineIntersection = GL1%LCInters(GL2,InterPoint1,InterPoint2)
            End Select
        Type is (Circle)
            !TODO
        End Select
    End Function GeneralLineIntersection
    
    !> ----------------------------------------------------------------------------------------------
    
    Subroutine InitializeLine(L,x1,y1,x2,y2)
        Class(Line),intent(inout) :: L
        Real(8),intent(in) :: x1,y1,x2,y2
        Call L%P1%init(x1,y1)
        Call L%P1%init(x2,y2)
    End Subroutine InitializeLine
    
    !Subroutine InitilizeLine(L,P1,P2)
    !    Class(Line),intent(inout) :: L
    !    Type(Point),intent(in) :: P1,P2
    !    Call L%P1%Copy(P1)
    !    Call L%P2%Copy(P2)
    !End Subroutine InitilizeLine
    
    
    Subroutine CopyLine(LD,LC)
        Class(Line),intent(inout) :: LD
        Class(Line),intent(in) :: LC
        Select Type (LD)
        Type is (Line)
            Select Type (LC)
            Type is (Line)
                LD = LC
            End Select
        End Select
    End Subroutine CopyLine
    
    Real(8) Function LengthLine(L) 
        Implicit none
        Class(Line),intent(in) :: L
        LengthLine = L%P1%DistanceP(L%P2)
        return 
    End Function LengthLine
    
    
    !> Return Value: -1 -- Cross Out of Segment || 0 -- Parallel || 1 -- Cross in Segment || 2 -- Coline
    Integer(4) Function LineIntersection(L1,L2,InterPoint)
        Class(Line),intent(in) :: L1, L2
        Class(Point),intent(out),optional :: InterPoint
        Real(8) :: A1,B1,C1,A2,B2,C2
        
        If (Present(InterPoint)) Call InterPoint%Init
        
        A1 = L1%P2%Y - L1%P1%Y
        B1 = L1%P1%X - L1%P2%X
        C1 = L1%P2%X * L1%P1%Y - L1%P1%X * L1%P2%Y
        A2 = L2%P2%Y - L2%P1%Y
        B2 = L2%P1%X - L2%P2%X
        C2 = L2%P2%X * L2%P1%Y - L2%P1%X * L2%P2%Y
        
        If ( Eq(A1 * B2, B1 * A2) ) Then
            If ( EQ((A1 + B1) * C2, (A2 + B2) * C1) ) Then
                ! coline
                LineIntersection = 2
            Else
                ! parallel
                LineIntersection = 0
            End if
        Else
            ! cross
            InterPoint%X = (B2 * C1 - B1 * C2) / (A2 * B1 - A1 * B2)
            InterPoint%Y = (A1 * C2 - A2 * C1) / (A2 * B1 - A1 * B2)
            IF (InterPoint%inLine(L1) .and. InterPoint%inLine(L2)) Then
                LineIntersection = 1
            Else
                LineIntersection = -1
            Endif
        End if
    End Function LineIntersection
    
    
    !> Return Value: -2 -- cross out of segment || -1 -- tangent out of segment|| 0 -- apart || 1 -- tangent in segment 
    !>                2 -- cross two points in segment || 10 -- cross one point in segment
    Integer(4) Function LineCircleIntersection(L,C,InterPoint1,InterPoint2)
        Class(Line),intent(in) :: L
        Class(Circle),intent(in) :: C
        Class(Point),intent(out),optional :: InterPoint1,InterPoint2
        LineCircleIntersection = 0
        !TODO
    End Function LineCircleIntersection
    
    
    Type(Vector) Function GetDirectionVector(L)
        Class(Line),intent(in) :: L
        GetDirectionVector%X = (L%P2%X - L%P1%X) / L%Length()
        GetDirectionVector%Y = (L%P2%Y - L%P1%Y) / L%Length()
        return      
    End Function GetDirectionVector
    
    
    !> Side: -1 -- left || 1 -- right || default -- right
    Type(Vector) Function GetPerpendicularDirectionVectorFromLine(L,Side)
        Class(Line),intent(in) :: L
        Integer(4),intent(in),optional :: Side
        Type(Vector) :: DV
        DV = L%DirectionV()
        If (present(Side)) Then
            GetPerpendicularDirectionVectorFromLine = DV%PerpenV(Side)
        Else
            GetPerpendicularDirectionVectorFromLine = DV%PerpenV()
        Endif
    End Function GetPerpendicularDirectionVectorFromLine
    
    !> ----------------------------------------------------------------------------------------------
    
    Real(8) Function PerimeterCircle(C) 
        !Use Constants, Only:Pi
        Implicit none
        Class(Circle),intent(in) :: C
        !PerimeterCircle = 2*Pi*(C%R)**2
        PerimeterCircle = 2*(C%R)**2
        return 
    End Function PerimeterCircle
    
End Module ModuleShape