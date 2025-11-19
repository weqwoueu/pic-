SUBROUTINE Linear_IFE_Basis_Eval_2D(coef, x, y, region_ind,				&
								    basis_ind, d_ind_x, d_ind_y, phi)

USE IFE_MAIN_PARAM

IMPLICIT NONE

REAL(8), INTENT(IN)					::	coef(4,8)
REAL(8), DIMENSION(:), INTENT(IN)	::	x, y
INTEGER, INTENT(IN)					::	basis_ind, d_ind_x, d_ind_y, &
										region_ind
REAL(8), DIMENSION(:), INTENT(OUT)	::	phi

IF (region_ind==1) THEN
	IF (d_ind_x == 1) THEN			! x derivative
	!	phi = coef(basis_ind, 1)
		phi = coef(basis_ind, 2) + coef(basis_ind, 4)* y
	ELSEIF (d_ind_y == 1) THEN		! y derivative
	!	phi = coef(basis_ind, 2)
        phi = coef(basis_ind, 3) + coef(basis_ind, 4)* x
	ELSE							! NO derivative
	!	phi = coef(basis_ind, 1)*x + coef(basis_ind, 2)*y + coef(basis_ind, 3)
        phi = coef(basis_ind, 1) + coef(basis_ind, 2)*x + coef(basis_ind, 3)*y +  coef(basis_ind, 4) * x* y
	END IF
ELSEIF (region_ind==2) THEN
	IF (d_ind_x == 1) THEN			! x derivative
!		phi = coef(basis_ind, 4)
		phi = coef(basis_ind, 6) + coef(basis_ind, 8)* y

	ELSEIF (d_ind_y == 1) THEN		! y derivative
!			phi = coef(basis_ind, 5)
        phi = coef(basis_ind, 7) + coef(basis_ind, 8)* x	
	
	ELSE							! NO derivative
!		phi = coef(basis_ind, 4)*x + coef(basis_ind, 5)*y + coef(basis_ind, 6)
        phi = coef(basis_ind, 5) + coef(basis_ind, 6)*x + coef(basis_ind, 7)*y +  coef(basis_ind, 8) * x* y
	END IF
END IF

END
