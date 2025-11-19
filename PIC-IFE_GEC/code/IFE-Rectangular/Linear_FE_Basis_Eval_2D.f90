SUBROUTINE Linear_FE_Basis_Eval_2D(coef, x, y, basis_ind, d_ind_x, d_ind_y, phi)

USE IFE_MAIN_PARAM

IMPLICIT NONE
																					   &
REAL(8), INTENT(IN)					::	coef(4,4)
REAL(8), DIMENSION(:), INTENT(IN)	::	x, y  !, z
INTEGER, INTENT(IN)					::	basis_ind, d_ind_x, d_ind_y  !, d_ind_z
REAL(8), DIMENSION(:), INTENT(OUT)	::	phi




IF (d_ind_x == 1) THEN			! x derivative
	phi = coef(basis_ind, 2) + coef(basis_ind, 4)* y
ELSEIF (d_ind_y == 1) THEN		! y derivative
    phi = coef(basis_ind, 3) + coef(basis_ind, 4)* x

ELSE							! NO derivative

    phi = coef(basis_ind, 1) + coef(basis_ind, 2)*x + coef(basis_ind, 3)*y +  coef(basis_ind, 4) * x* y

END IF


END
