
SUBROUTINE	GTOP(i_part,parapre,para)

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Particle_2D
USE Field_2D
IMPLICIT NONE

REAL(8)		            ::	para(0:nx+1,0:ny+1)
INTEGER					::	i_part,i,j,delta
REAL(8)					::	parapre, xp, yp, dx, dy
REAL(8)		            xcellmdx, ycellmdy, R1, R2, den

!delta = 0 


    parapre=0 
   
	xp = (part(i_part,1) - Vert_o(1))*hxi(1)
	i = xp
	dx = xp-i

	yp = (part(i_part,2) - Vert_o(2))*hxi(2)
	j = yp
	dy = yp-j

	xcellmdx = 1. -dx
	ycellmdy = 1. -dy
		
	!IF (delta == 0) THEN					!Cartesian coordinates
	IF (delta_global == 0) THEN					!$ mb.ZWZ 2021/7/11
		parapre = parapre + xcellmdx*ycellmdy*para(i,j)
		parapre = parapre + dx*ycellmdy*para(i+1,j)
		parapre = parapre + xcellmdx*dy*para(i,j+1)
		parapre = parapre + dx*dy*para(i+1,j+1)

	!ELSEIF( delta == 1 ) THEN
	ELSEIF( delta_global == 1 ) THEN               !$ mb.ZWZ 2021/7/11

        R1=float(j - 1)*hx(2)
        R2=float(j)*hx(2)
        den=R2*R2-R1*R1

		parapre = parapre + para(i,j)*xcellmdx*(R2*R2-part(i_part,2)*part(i_part,2))/den
		parapre = parapre + para(i+1,j)*dx*(R2*R2-part(i_part,2)*part(i_part,2))/den
		parapre = parapre + para(i,j+1)*xcellmdx*(part(i_part,2)*part(i_part,2)-R1*R1)/den
		parapre = parapre + para(i+1,j+1)*dx*(part(i_part,2)*part(i_part,2)-R1*R1)/den

	ENDIF

END SUBROUTINE
