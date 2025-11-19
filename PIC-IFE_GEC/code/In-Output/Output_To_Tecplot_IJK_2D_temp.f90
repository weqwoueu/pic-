SUBROUTINE Output_To_Tecplot_IJK_2D_temp( nnx, nny, VertX, Phi, Rho, Output_Option)


IMPLICIT NONE

INTEGER						::	nnx, nny

REAL(8)		::	VertX(2,0:nnx+1,0:nny+1), Phi(0:nnx+1,0:nny+1), Rho(0:nnx+1,0:nny+1) !, Rho_s(0:nnx+1,0:nny+1,1:ispe_max)

INTEGER							Output_Option

INTEGER			i, j, nindex, isp
CHARACTER*20	fname, sname, filename

PRINT*, "Saving Data in Tecplot IJ Format ...."

WRITE(fname,999) 0
999 FORMAT(I6.6)

Check_Option: SELECT CASE (Output_Option)
CASE (1)

filename = 'field_IJ_'//TRIM(fname)//'.dat'
PRINT*, 'Writing to file: ', TRIM(filename)
OPEN(1, ACTION = 'WRITE', FILE = TRIM(filename))
WRITE(1,*) 'TITLE = "Field Plot"'
WRITE(1,*) 'VARIABLES = "x" "y" "Phi" "Rho"'
WRITE(1,40) nnx, nny
	DO j=1,nny
		DO i=1,nnx
			WRITE(1,50) VertX(1:2,i,j), Phi(i,j), Rho(i,j)
		END DO
	END DO
40 FORMAT (' ZONE I = ',I6,', J= ',I6)
50 FORMAT (E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
CLOSE(1)





END SELECT Check_Option

PRINT*, "Saving Data in Tecplot IJK Format Done."

END SUBROUTINE
