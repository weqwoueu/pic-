SUBROUTINE Loadv(v, energy, insp)

USE Particle_2D
IMPLICIT NONE


REAL(8)		v, energy
INTEGER		insp
REAL(8)		VMP_Sheath, pacc
!REAL(8)		ranum, ranum1
double precision ::  ranum, ranum1
REAL(8)        RF

VMP_Sheath = SQRT(2.* energy /xm(insp))

100	CALL DRandom(ranum)
!100     ranum = randum()
	v = (6 * ranum - 3) * VMP_Sheath
	pacc = exp(-v**2/(VMP_Sheath**2))
    CALL DRandom(ranum1)
	IF(ranum1 > pacc) THEN
		goto 100
	ENDIF
END SUBROUTINE