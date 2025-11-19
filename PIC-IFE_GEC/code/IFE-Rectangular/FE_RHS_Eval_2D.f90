SUBROUTINE FE_RHS_Eval_2D( RHS_FUN, U_val, R_val, el_region, n_nodes_in_elem_1, &
						   rhs_hw, interpol_hw, vector, gnodes, delta )

USE IFE_MAIN_PARAM

IMPLICIT NONE

EXTERNAL						RHS_FUN
INTEGER, INTENT(IN)			::	el_region, n_nodes_in_elem_1
REAL(8), INTENT(IN)			::	rhs_hw(4,4), interpol_hw(4,4)

REAL(8), INTENT(IN)			::	gnodes(2,4)

REAL(8), INTENT(OUT)		::	vector(4)
REAL(8), INTENT(IN)			::	U_val(n_nodes_in_elem_1), R_val(n_nodes_in_elem_1)

REAL(8)									::	eint
INTEGER									::	i,delta
REAL(8), DIMENSION(:), ALLOCATABLE		::	rhs_val, u_gnodes, r_gnodes
											

ALLOCATE(u_gnodes(4), r_gnodes(4))
u_gnodes = Zero
r_gnodes = Zero
DO i=1,4
	u_gnodes = u_gnodes + u_val(i)*interpol_hw(i,:)
	r_gnodes = r_gnodes + r_val(i)*interpol_hw(i,:)
END DO
ALLOCATE(rhs_val(4))

CALL RHS_FUN(r_gnodes, u_gnodes, el_region, rhs_val, gnodes(1, :), gnodes(2, :), SIZE(u_gnodes))

vector = Zero

DO i=1,n_nodes_in_elem_1


	eint = SUM(rhs_hw(i,:)*rhs_val)   

	vector(i) = vector(i) + eint

END DO

DEALLOCATE(rhs_val)
DEALLOCATE(u_gnodes, r_gnodes)

END