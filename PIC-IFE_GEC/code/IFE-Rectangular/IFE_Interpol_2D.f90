SUBROUTINE IFE_Interpol_2D(vert, u_val, x_int, y_int, &
						Intrs_pts, el_type, el_beta, sub_region_ind, dind, f)

USE IFE_MAIN_PARAM
USE IFE_INTERFACE, ONLY: Linear_IFE_Basis_Eval_2D

IMPLICIT NONE

REAL(8), INTENT(IN)					::	vert(2,4), el_beta(2)
REAL(8), DIMENSION(:,:), INTENT(IN)	::	Intrs_pts
REAL(8), DIMENSION(:), INTENT(IN)	::	u_val, x_int, y_int
INTEGER, INTENT(IN)					::	dind(2), el_type, sub_region_ind
REAL(8), DIMENSION(:), INTENT(OUT)	::	f

INTEGER	i
REAL(8), DIMENSION(:), ALLOCATABLE		::	basis
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	basis_coef, intrs_cpy

f = Zero

ALLOCATE(basis(SIZE(x_int)))

basis = Zero

ALLOCATE(basis_coef(4,8))
ALLOCATE(intrs_cpy(2,2))

CALL Linear_IFE_Basis_Coeff_2D(vert, Intrs_pts, el_type, el_beta, basis_coef)

DO i=1,4
	
	CALL Linear_IFE_Basis_Eval_2D(	basis_coef, x_int, y_int, &
					                sub_region_ind, i, dind(1), dind(2), basis)

	f = f + u_val(i)*basis

END DO

DEALLOCATE(basis, basis_coef)
DEALLOCATE(intrs_cpy)

END