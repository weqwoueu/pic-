SUBROUTINE	PTOG(para,i_part)

USE MCC_Data_2D
USE IFE_Data
USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Particle_2D
USE Field_2D
IMPLICIT NONE

REAL(8)		            ::	para(0:nx+1,0:ny+1)
INTEGER					::	i_part,i,j,delta
REAL(8)					::	rxp, ryp, dx, dy
REAL(8)		            xcellmdx, ycellmdy, R1, R2, den

!delta = 0

	rxp = (part(i_part,1) - Vert_o(1))*hxi(1)
	i = rxp
	dx = rxp-i

	ryp = (part(i_part,2) - Vert_o(2))*hxi(2)
	j = ryp
	dy = ryp-j

	xcellmdx = 1. -dx
	ycellmdy = 1. -dy

	!IF( delta == 0 ) THEN
	IF( delta_global == 0 ) THEN           !$ mb.ZWZ 2021/7/11
		para(i  ,j  ) = para(i  ,j  ) + part(i_part,11)*xcellmdx*ycellmdy
		para(i+1,j  ) = para(i+1,j  ) + part(i_part,11)*dx*ycellmdy
		para(i  ,j+1) = para(i  ,j+1) + part(i_part,11)*xcellmdx*dy
		para(i+1,j+1) = para(i+1,j+1) + part(i_part,11)*dx*dy

	!ELSEIF( delta == 1 ) THEN
	ELSEIF( delta_global == 1 ) THEN       !$ mb.ZWZ 2021/7/11

        R1=float(j - 1)*hx(2)
        R2=float(j)*hx(2)
        den=R2*R2-R1*R1

		para(i  ,j  ) = para(i  ,j  ) + part(i_part,11)*xcellmdx*(R2*R2-part(i_part,2)*part(i_part,2))/den
		para(i+1,j  ) = para(i+1,j  ) + part(i_part,11)*dx*(R2*R2-part(i_part,2)*part(i_part,2))/den
		para(i  ,j+1) = para(i  ,j+1) + part(i_part,11)*xcellmdx*(part(i_part,2)*part(i_part,2)-R1*R1)/den
		para(i+1,j+1) = para(i+1,j+1) + part(i_part,11)*dx*(part(i_part,2)*part(i_part,2)-R1*R1)/den

	ENDIF
	

END SUBROUTINE
