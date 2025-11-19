SUBROUTINE GetParFlux_2D_Temp

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D
USE Particle_2D

IMPLICIT NONE

INTEGER		i, j, isp


IF (.NOT.ALLOCATED(vj_s1)) THEN
	ALLOCATE(vj_s1(2,0:nx+1,0:ny+1,ispe_tot), vj_s2(2,0:nx+1,0:ny+1,ispe_tot))
END IF

DO isp=1,ispe_tot
	DO j=0, ny+1
		DO i=0, nx+1
			vj_s1(1,i,j,isp) = vj_sx(i,j,isp)
			vj_s1(2,i,j,isp) = vj_sy(i,j,isp)
			vj_s2(1,i,j,isp) = vj_sx(i,j,isp)**2
			vj_s2(2,i,j,isp) = vj_sy(i,j,isp)**2
		END DO          
	END DO
END DO
      
END