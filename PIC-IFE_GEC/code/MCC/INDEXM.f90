SUBROUTINE	INDEXM(insp, IC, nnx, nny) 

USE PIC_Data
USE DSMC_Data_2D
USE IFE_Data_2D
USE Particle_2D     !$ ab.ZWZ 20210709
USE Domain_2D       !$ ab.ZWZ 20210709

IMPLICIT NONE

INTEGER								::	insp, nsp, nnx, nny
INTEGER, DIMENSION(2,nnx * nny)		::	IC
INTEGER								::	i, j, k, MC, ipret, jpret
REAL(8)								::	x, y, Vert_o(2)
INTEGER, DIMENSION(:), ALLOCATABLE		::	ISC


Vert_o(1) = Domain%xmin - hx
Vert_o(2) = Domain%ymin - hy


IC = 0

ALLOCATE(ISC(ntot))
ISC=0

DO i=1, ntot
	IF(part(i,7) .NE. 0) THEN
		x = part(i,1)
		y = part(i,2)
		nsp = part(i,7)
		IF(nsp == insp) THEN
			ipret = INT((part(i,1) - Vert_o(1)) / hx)
			jpret = INT((part(i,2) - Vert_o(2)) / hy)
!			IF(ipret > nnx .OR. ipret <= 0 .OR. jpret > nny .OR. jpret <= 0) THEN
!				DO ip=1,SIZE(part,2)
!					part(i, ip) = part(npt, ip)
!				END DO
!					ns(nsp) = ns(nsp) - 1
!					npt=npt-1
!					goto 100
!			ENDIF
			MC = (ipret - 1) * nny + jpret
			ISC(i)=MC
			IC(2,MC) = IC(2,MC) + 1     !MC号网格里的粒子个数
		ENDIF
100	ENDIF
ENDDO


IF(insp>=3)THEN
	j=0
	DO i=1,nnx * nny
		IC(1,i)=j !IC(1,N),第N个网格前的粒子总数-1
		j=j+IC(2,i)
	ENDDO

	IC(2,:)=0

	DO i=1, ntot
	    IF(part(i,7) == insp) THEN
		    MC=ISC(i)
		    IC(2,MC)=IC(2,MC)+1
		    k=IC(1,MC)+IC(2,MC)
    		IR(k)=i !IR,得到原始号码为N的分子根据其新位置的编号K
    	ENDIF
	ENDDO
ENDIF

DEALLOCATE(ISC)
END SUBROUTINE