SUBROUTINE Load_Charged_Particle(delta, xmin, xmax, ymin, ymax)

USE Field_2D
!USE MCC_Main_Param
USE Constant_Variable_2D
USE Particle_2D
USE PIC_MAIN_PARAM_2D

IMPLICIT NONE

INTEGER :: delta

REAL(8)		xmin, xmax, ymin, ymax, RANUM

REAL(8) :: Vol_Region_M3 = 0.0

INTEGER(8) :: Num_Load

Integer(4) :: nnp

INTEGER :: INSP, J



!Vol_Region_M3 = PI*(ymax*ymax-ymin*ymin)*(xmax-32.5)*L_ref**3
Vol_Region_M3 = (xmax-xmin)*L_ref * PI*(ymax-ymin)**2 *L_ref**2

print*,dens0(ipf(1))
print*,n_ref
print*,affp_bjw(ipf(1))
print*,Vol_Region_M3

Num_Load = dens0(ipf(1))*n_ref * Vol_Region_M3 / affp_bjw(ipf(1))
!Num_Load = dens0(ipf(1))*n_ref * Vol_Region_M3
!Num_Load = dens0(ipf(1)) * n_ref * Vol_Region_M3/10000
!print*,dens0(ipf(1)) * n_ref * Vol_Region_M3
print*,Num_Load
print*,REAL(Num_Load/80.0/160.0)

!Num_load = 128*(xmax-xmin)/0.25*(ymax-ymin)/0.25
!Num_load = 16*(xmax-xmin)*(ymax-ymin)
!
!print*,Num_Load



INSP = ipf(1)
DO J=ntot+1, ntot+Num_Load
    CALL RANDOM_NUMBER(RANUM)
    part(J,1) = xmin + (xmax-xmin)*RANUM
    CALL RANDOM_NUMBER(RANUM)
    IF (delta == 0) THEN
        part(J,2) = ymin + (ymax-ymin)*RANUM
        part(J,3) = 0.
    ELSE IF(delta == 1) THEN
        part(J,2) = SQRT(ymin*ymin+(ymax+ymin)*(ymax-ymin)*RANUM)
        CALL RANDOM_NUMBER(RANUM)
        part(J,3) = 2*PI*RANUM
    ENDIF
    CALL LOADV(part(J,4), tmpj(INSP), INSP)
    CALL LOADV(part(J,5), tmpj(INSP), INSP)
    CALL LOADV(part(J,6), tmpj(INSP), INSP)
    part(J,7) = 1
ENDDO
ns(insp) = ns(insp) + Num_Load
ntot = ntot + Num_Load

INSP = ipf(2)
DO J=ntot+1, ntot+Num_Load
    CALL RANDOM_NUMBER(RANUM)
    part(J,1) = xmin + (xmax-xmin)*RANUM
    CALL RANDOM_NUMBER(RANUM)
    IF (delta == 0) THEN
        part(J,2) = ymin + (ymax-ymin)*RANUM
        part(J,3) = 0.
    ELSE IF(delta == 1) THEN
        part(J,2) = SQRT(ymin*ymin+(ymax+ymin)*(ymax-ymin)*RANUM)
        CALL RANDOM_NUMBER(RANUM)
        part(J,3) = 2*PI*RANUM
    ENDIF
    CALL LOADV(part(J,4), tmpj(INSP), INSP)
    CALL LOADV(part(J,5), tmpj(INSP), INSP)
    CALL LOADV(part(J,6), tmpj(INSP), INSP)
    part(J,7) = 2
ENDDO
ns(insp) = ns(insp) + Num_Load
ntot = ntot + Num_Load

!DO j=1, ns(INSP)
!    IF( part(j,1) .LE. 32.5 ) THEN
!        part(j,1) = -2000
!    ENDIF
!ENDDO

!$ === ab. ZWZ 2020/10/29 m====
ntot =0
DO insp=1, ispe_tot
    ntot = ntot + ns(insp)
ENDDO
!$ === ab. ZWZ 2020/10/29 m====
END SUBROUTINE
