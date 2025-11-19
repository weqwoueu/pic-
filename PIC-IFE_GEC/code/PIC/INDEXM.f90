SUBROUTINE	INDEXM(insp, IC, nnx, nny)

USE Particle_2D
USE Domain_2D
IMPLICIT NONE

INTEGER								::	insp, nsp, nnx, nny
INTEGER, DIMENSION(2,nnx * nny)		::	IC
INTEGER								::	i, j, k, MC, ipret, jpret, num, ip
REAL(8)								::	x, y
INTEGER								::	npt

IC = 0

j = 1

num = SIZE(part, 1)
npt = ntot

DO i=1, ntot
	IF(part(i,7) .NE. 0) THEN
		x = part(i,1)
		y = part(i,2)
		nsp = part(i,7)
		IF(nsp == insp) THEN
			ipret = INT((x - 0) / hx(1)) + 1
			jpret = INT((y - 0) / hx(2)) + 1
			!IF(ipret > nnx .OR. ipret <= 0 .OR. jpret > nny .OR. jpret <= 0) THEN
			IF(ipret > nnx .OR. ipret <= 0 .OR. jpret > nny .OR. jpret <= 0 .OR. x <0 .OR. y<0) THEN
				DO ip=1,SIZE(part,2)
					part(i, ip) = part(npt, ip)
				END DO
				ns(nsp) = ns(nsp) - 1
				npt=npt-1
				goto 100
			ENDIF
			MC = (ipret - 1) * (ny - 1) + jpret
			IC(2,MC) = IC(2,MC) + 1
		ENDIF
100	ENDIF
ENDDO

ntot = npt



END SUBROUTINE