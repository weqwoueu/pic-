!SUBROUTINE Weighting(P1,P2,P3,P4)
!    
!    IMPLICIT NONE
!
!    INTEGER :: delta
!    REAL(8) :: R1, R2, R , DEN
!    REAL(8),INTENT(OUT) :: P1, P2, P3, P4
!
!    
!    IF (delta == 0) THEN
!        P1 = xcellmdx*ycellmdy
!        P2 = dx*ycellmdy
!        P3 = xcellmdx*dy
!        P4 = dx*dy      
!    ELSEIF (delta == 1) THEN
!        R1 = YSTART + FLOAT(J-1)*hx(2)
!        R2 = YSTART + float(J)  *hx(2)
!        R = 
!        DEN = R2*R2 - R1*R1
!        
!        P1 = (1.-XTEM) * (R2*R2 - R * R)/DEN
!        P2 = XTEM      * (R2*R2 - R * R)/DEN
!        P3 = (1.-XTEM) * (R * R - R1*R1)/DEN
!        P4 = XTEM      * (R * R - R1*R1)/DEN
!    ENDIF
!    
!END SUBROUTINE
!    
!SUBROUTINE PRE_INTERP(xp,x0,dx,nx,yp,y0,dy,ny,XTEM,YTEM,IPRE,JPRE)
!
!    IMPLICIT NONE
!
!    REAL(8), INTENT(IN) :: xp, dx, yp, dy
!    REAL(8), INTENT(OUT) :: XTEM, YTEM
!    INTEGER, INTENT(OUT) :: IPRE, JPRE
!    REAL(8) :: Ferror = 1.0E-5
!
!    XTEM = (xp - x0) /dx
!	IPRE = FLOOR(xp)        !!! 粒子在第i个单元里
!	XTEM = XTEM - IPRE
!    
!    YTEM = (yp - y0) /dy
!	JPRE = FLOOR(yp)        !!! 粒子在第j个单元里
!	YTEM = YTEM - JPRE
!    
!    IF (IPRE == nx) THEN
!        IPRE = IPRE-1
!        XTEM = 1.-FERROR
!    END IF
!    IF (JPRE == ny) THEN
!        JPRE = JPRE-1
!        YTEM= 1.-FERROR
!    END IF
!
!
!END SUBROUTINE