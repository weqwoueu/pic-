SUBROUTINE		GAM(X,Z)

IMPLICIT NONE

REAL(8), INTENT(IN)		::	X
REAL(8), INTENT(OUT)		::	Z
REAL(8)				    	::	Y,A

!--calculates the Gamma function of X.
      A=1.
      Y=X
      IF (Y.LT.1.) THEN
        A=A/Y
      ELSE
50      Y=Y-1
        IF (Y.GE.1.) THEN
          A=A*Y
          GO TO 50
        END IF
      END IF
      Z=A*(1.-0.5748646*Y+0.9512363*Y**2-0.6998588*Y**3+ &
         0.4245549*Y**4-0.1010678*Y**5)


END
