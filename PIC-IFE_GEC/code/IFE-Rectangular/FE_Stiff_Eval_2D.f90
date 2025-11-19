SUBROUTINE FE_Stiff_Eval_2D( el_beta, n_nodes_in_elem_1, n_nodes_in_elem_2, stiff_hw, matrix )

USE IFE_MAIN_PARAM

IMPLICIT NONE

INTEGER, INTENT(IN)						::	n_nodes_in_elem_1, n_nodes_in_elem_2
REAL(8), INTENT(IN)						::	el_beta, stiff_hw(4,4)
REAL(8), INTENT(OUT)					::	matrix(4,4)

INTEGER									::	i, j
											
matrix = Zero
DO i=1,n_nodes_in_elem_1
   		            
   	DO j=1,n_nodes_in_elem_2

		matrix(i,j) = el_beta*stiff_hw(i,j)
	END DO

END DO


END