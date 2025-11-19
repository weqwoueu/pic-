SUBROUTINE FE_Stiff_HW_2D( Coeff_FUN, vert, n_nodes_in_elem_1, n_nodes_in_elem_2, hardwire, delta )

USE IFE_MAIN_PARAM
USE IFE_Data
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Linear_FE_Basis_Eval_2D, Generate_Gauss_local

IMPLICIT NONE

EXTERNAL									Coeff_FUN
REAL(8), INTENT(IN)						::	vert(2,4)
INTEGER, INTENT(IN)						::	n_nodes_in_elem_1, n_nodes_in_elem_2
REAL(8), INTENT(OUT)					::	hardwire(4,4)

REAL(8)									::	V, eint_x, eint_y, eint
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	bas_val_x, bas_val_y
REAL(8), DIMENSION(:), ALLOCATABLE		::	gwght, coeff_val
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	gnodes, basis_coef
INTEGER									::	i, j
											
INTEGER, DIMENSION(:), ALLOCATABLE		::	dind
INTEGER                                 ::  delta


!gwght = (/1.0/Three, 1.0/Three, 1.0/Three/)
!=======================================================================================
INTEGER                                           Gauss_point_number
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference_triangle
!REAL(8)	                                          h_partition(2)
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_local
!=========================================================================================
ALLOCATE(gwght(4))  


CALL Generate_Gauss_local(Gauss_coefficient_reference_Four,Gauss_point_reference_Four, &
                            vert(1:2,1),h_partition, Gauss_coefficient_local_Four, Gauss_point_local_Four)

gwght = Gauss_coefficient_local_Four

ALLOCATE(bas_val_x(4,4), bas_val_y(4,4))

ALLOCATE(gnodes(2,4))

gnodes=TRANSPOSE(Gauss_point_local_Four)

ALLOCATE(coeff_val(4))
CALL Coeff_FUN(gnodes(1,:), gnodes(2,:), coeff_val, SIZE(gnodes, 2))

ALLOCATE(dind(2))

ALLOCATE(basis_coef(4,4))
CALL Linear_FE_Basis_Coeff_2D(vert, basis_coef)

hardwire = Zero
DO i=1,n_nodes_in_elem_1
    
	dind = (/1, 0/)
	CALL Linear_FE_Basis_Eval_2D(	basis_coef, gnodes(1,:), gnodes(2,:), &
								    i, dind(1), dind(2), bas_val_x(i,:))

    dind = (/0, 1/)
	CALL Linear_FE_Basis_Eval_2D(	basis_coef, gnodes(1,:), gnodes(2,:), &
							        i, dind(1), dind(2), bas_val_y(i,:))

END DO

DO i=1,n_nodes_in_elem_1
   	DO j=1,n_nodes_in_elem_2
		IF (delta == 0) THEN
		    eint_x = SUM(gwght*coeff_val*bas_val_x(i,:)*bas_val_x(j,:))
		    eint_y = SUM(gwght*coeff_val*bas_val_y(i,:)*bas_val_y(j,:))
		ELSEIF (delta ==1) THEN
		    eint_x = SUM(gwght*coeff_val*gnodes(2,:)*bas_val_x(i,:)*bas_val_x(j,:))
		    eint_y = SUM(gwght*coeff_val*gnodes(2,:)*bas_val_y(i,:)*bas_val_y(j,:))
		ELSE 
		  PRINT*, ' delta value is wrong, check again, stop'
		  STOP
		ENDIF
		eint = eint_x + eint_y 

		hardwire(i,j) = hardwire(i,j) + eint

	END DO
END DO

DEALLOCATE(gwght)
DEALLOCATE(bas_val_x, bas_val_y)
DEALLOCATE(gnodes, basis_coef)
DEALLOCATE(coeff_val)
DEALLOCATE(dind)


END