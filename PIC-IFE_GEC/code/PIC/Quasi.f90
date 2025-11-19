SUBROUTINE	Quasi

USE Domain_2D
!use ifport
USE Particle_2D
IMPLICIT NONE

double precision randum
external randum

INTEGER									::	n_cell_in = 2
INTEGER, DIMENSION(2,(nx-1)*(ny-1))		::	ic1, IC
INTEGER									::	insp
INTEGER									::	nnx, nny
INTEGER									::	ipre, jpre
INTEGER, DIMENSION(nx-1,ny-1)			::	net_at_cell
REAL(8), DIMENSION(nx-1,ny-1)			::	pro_at_cell
INTEGER									::	net_at_boundary = 0
INTEGER									::	mc, num_ion, num_elec, num
INTEGER									::	i, j
REAL(8)									::	tmpj_elec
REAL									::	ranum
REAL(8)									::	x, y, u, v
REAL(8)                                    ::  RF

nnx = nx - 1
nny = ny - 1
tmpj_elec = TMPJ(ipf(1))

DO insp = 1, 2
	CALL INDEXM(insp, IC, nnx, nny)
	ic1(insp,:) = IC(2,:)
ENDDO

!	OPEN(1, ACTION = 'WRITE', FILE = 'IC1_1.dat')
!    DO i=1, nnx*nny
!        WRITE(1,50) IC1(1,:)
!    ENDDO
!    50 FORMAT (I2)
!    CLOSE(1)
!
!
!	OPEN(1, ACTION = 'WRITE', FILE = 'IC1_2.dat')
!    DO i=1, nnx*nny
!        WRITE(1,50) IC1(2,:)
!    ENDDO
!    CLOSE(1)


DO jpre=nnx - n_cell_in, nny
	net_at_cell = 0
	net_at_boundary = 0
	pro_at_cell = 0
	DO ipre = 1, nnx
		mc = (ipre - 1) * nny + jpre
		net_at_cell(ipre,jpre) = ic1(2,mc) - ic1(1,mc)
		net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
	ENDDO

	IF(net_at_boundary > 0) THEN
		net_at_boundary = 0
		DO ipre = 1, nnx
			IF(net_at_cell(ipre,jpre) > 0) THEN
				pro_at_cell(ipre,jpre) = REAL(net_at_cell(ipre,jpre))
				net_at_boundary = net_at_boundary + net_at_cell(ipre,jpre)
			else
				pro_at_cell(ipre,jpre) = 0.
			ENDIF
		ENDDO

		pro_at_cell(:,jpre) = pro_at_cell(:,jpre) / REAL(net_at_boundary)

		DO ipre = 2, nnx
			pro_at_cell(ipre,jpre) = pro_at_cell(ipre,jpre) + pro_at_cell(ipre-1,jpre)
		ENDDO

		DO i = 1, net_at_boundary
			num_elec = ns(1)
			num_ion = ns(2)
			num = num_elec + num_ion
!			CALL random(ranum)
            ranum = randum()
			ipre = 1
			DO WHILE(ranum > pro_at_cell(ipre, jpre) .AND. ipre < nny)
				ipre = ipre + 1
			ENDDO

!			CALL random(ranum)
            ranum = randum()
			part(num+1,2) = 0 + (jpre - ranum) * hx(1)
			
!			CALL random(ranum)
            ranum = randum()
			part(num+1,1) = 0 + (ipre - ranum) * hx(2)

			part(num+1,3) = 50

			CALL Loadv(u, tmpj_elec, 1)
			part(num+1,4) = u

100			CALL Loadv(v, tmpj_elec, 1)
			IF(v >= 0) THEN
				goto 100
			ENDIF
			
			part(num+1,5) = v

			part(num+1,6) = 0

			part(num+1,7) = 1

			ns(1) = ns(1) + 1
		ENDDO
	ENDIF
ENDDO

ntot = ns(1) + ns(2)

!ntot = 9739
!
!part(1:9739,:) = 0
!
!	OPEN(1, ACTION = 'READ', FILE = 'Phase236.dat')
!    READ(1,*) 
!    READ(1,*) 
!    DO i=1, ntot
!        READ(1,50) part(i,1), part(i,2), part(i,4), part(i,5), part(i,7)
!    ENDDO
!    50 FORMAT (E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',F15.0)
!    CLOSE(1)
!
!    DO i=1, ntot
!        part(i,2)=part(i,2)+20
!    ENDDO
!   

END SUBROUTINE