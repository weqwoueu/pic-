SUBROUTINE Check_Sub_Nodes_2D(Subvert1, Subvert2, BCvert1, BCvert2, IN1, IN2 )

IMPLICIT NONE

REAL(8), INTENT(IN)		::	Subvert1(2), Subvert2(2), BCvert1(2), BCvert2(2)
INTEGER, INTENT(OUT)	::	IN1, IN2

REAL(8)  :: P1(2), P2(2), Dist1, Dist2, SmallValue 

SmallValue = 1.0D-5

P1 = BCvert1 - Subvert1
P2 = BCvert2 - Subvert1

Dist1 = DSQRT(P1(1)*P1(1) + P1(2)*P1(2))
Dist2 = DSQRT(P2(1)*P2(1) + P2(2)*P2(2))

IN1 = -1

IF (Dist1<SmallValue) IN1 = 1
IF (Dist2<SmallValue) IN1 = 2


IF (Dist1>=SmallValue .AND. Dist2>=SmallValue) THEN
   Dist1 = P1(1)*P2(2) - P1(2)*P2(1)
   IF ( Dist1 == 0.0) THEN ! Point on line
      IN1 = 0
   ENDIF
ENDIF 

P1 = BCvert1 - Subvert2
P2 = BCvert2 - Subvert2

Dist1 = DSQRT(P1(1)*P1(1) + P1(2)*P1(2))
Dist2 = DSQRT(P2(1)*P2(1) + P2(2)*P2(2))

IN2 = -1

IF (Dist1<SmallValue) IN2 = 1
IF (Dist2<SmallValue) IN2 = 2


IF (Dist1>=SmallValue .AND. Dist2>=SmallValue) THEN
   Dist1 = P1(1)*P2(2) - P1(2)*P2(1)
   IF ( Dist1 == 0.0) THEN ! Point on line
      IN2 = 0
   ENDIF
ENDIF 

END