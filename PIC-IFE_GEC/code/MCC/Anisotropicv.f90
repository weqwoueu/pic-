SUBROUTINE Anisotropicv(vel,vz,vr,vphi)

USE IFE_MAIN_PARAM

IMPLICIT NONE

!double precision          randum
!external                   randum

REAL(8)					    RANUM
    
REAL(8)					::	vel,vz,vr,vphi,coschi,sinchi,phi1
	

CALL DRandom(RANUM)
coschi=2*RANUM-1
sinchi=sqrt(1-coschi**2)
CALL DRandom(RANUM)

phi1=2.0*pi*RANUM
vz=vel*sinchi*cos(phi1)
vr=vel*sinchi*sin(phi1)
vphi=vel*coschi

!write(*,*) vz,vr,vphi
!PAUSE
END SUBROUTINE