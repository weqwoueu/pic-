SUBROUTINE FE_RHS_HW_2D( vert, gwght, gnodes, n_nodes_in_elem_1, hardwire, delta)

USE IFE_MAIN_PARAM
USE IFE_INTERFACE, ONLY: Linear_FE_Basis_Eval_2D, Linear_FE_Basis_Coeff_2D

IMPLICIT NONE

REAL(8), INTENT(IN)						::	vert(2,4), gnodes(2,4),gwght(4)
INTEGER, INTENT(IN)						::	n_nodes_in_elem_1, delta

REAL(8), INTENT(OUT)					::	hardwire(4,4)

REAL(8), DIMENSION(:), ALLOCATABLE		::	eint, ibas_val
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	basis_coef
INTEGER									::	i
											
INTEGER, DIMENSION(:), POINTER			::	dind



ALLOCATE(ibas_val(4))

ALLOCATE(dind(2))
dind = (/0, 0/);

ALLOCATE(eint(4))
!
ALLOCATE(basis_coef(4,4))
CALL Linear_FE_Basis_Coeff_2D(vert, basis_coef)

hardwire = Zero
DO i=1,n_nodes_in_elem_1

	CALL Linear_FE_Basis_Eval_2D(basis_coef, gnodes(1,:), gnodes(2,:), i,		&
							     dind(1), dind(2), ibas_val)

 	IF (delta == 0) THEN
   	  eint = gwght*ibas_val;      
    ELSEIF(delta == 1) THEN
   	  eint = gwght*ibas_val*gnodes(2,:);      
	ELSE 
	  PRINT*, ' delta value is wrong, check again, stop'
	  STOP
	ENDIF
                 
	hardwire(i,:) = hardwire(i,:) + eint;

END DO

DEALLOCATE(ibas_val)
DEALLOCATE(basis_coef)
DEALLOCATE(dind)
DEALLOCATE(eint)

END