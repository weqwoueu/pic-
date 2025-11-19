SUBROUTINE Node_in_Boundary_2D(point,line_begin, line_end, IN)

IMPLICIT NONE

REAL(8), INTENT(IN)		::	point(2),line_begin(2), line_end(2)
LOGICAL, INTENT(OUT)	::	IN

REAL(8) :: P1(2), P2(2), Dist1, Dist2, SmallValue 

SmallValue = 1.0D-5

P1 = line_begin - point
P2 = line_end - point

Dist1 = DSQRT(P1(1)*P1(1) + P1(2)*P1(2))
Dist2 = DSQRT(P2(1)*P2(1) + P2(2)*P2(2))

IN = .FALSE.

IF (Dist1<SmallValue .OR. Dist2<SmallValue) THEN
  IN =.TRUE.
ELSE
   Dist1 = P1(1)*P2(2) - P1(2)*P2(1)
   IF ( ABS(Dist1) < SmallValue) THEN ! Point on line
      Dist2 = P1(1)*P2(1) + P1(2)*P2(2)
	  IF ( Dist2 < 0.0) THEN  ! Point on line segment
		  IN =.TRUE.
	  ENDIF
   ENDIF
ENDIF 

!if (IN) then
!print*,  'point=',point
!print*,  'line_begin=', line_begin
!print*,  'line_end=', line_end
!print*,  'IN=',IN
!print*, P1, P2
!!pause
!endif

END