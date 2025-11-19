SUBROUTINE IFE_Stiff_2D(	tri_sign, pointer_reference_to_local, Coeff_FUN, vert, Intrs_pts, el_type, el_beta, &
					  	    n_nodes_in_elem_1, n_nodes_in_elem_2, matrix, delta )

USE IFE_MAIN_PARAM
USE IFE_Data
USE IFE_INTERFACE, ONLY: Inter_Tetra_Partition_2D, Gauss_Nodes_2D, Linear_IFE_Basis_Eval_2D, &
                        Generate_Gauss_reference, Generate_Gauss_local, Linear_IFE_Basis_Coeff_2D

IMPLICIT NONE

EXTERNAL									Coeff_FUN
REAL(8), INTENT(IN)						::	vert(2,4), el_beta(2)
INTEGER, INTENT(IN)						::	n_nodes_in_elem_1, n_nodes_in_elem_2,	&
											el_type, delta
INTEGER									    tri_sign(2)
REAL(8), INTENT(IN)						::	Intrs_pts(2, 2)
REAL(8), INTENT(OUT)					::	matrix(4,4)

REAL(8)									::	V, eint_x, eint_y, eint,intrs_cpy_1(2,2)
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	bas_val_x, bas_val_y
REAL(8), DIMENSION(:), ALLOCATABLE		::	gwght, coeff_val
REAL(8), DIMENSION(:,:), POINTER		::	p_sub_pt, g_sub_x, g_sub_y
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	gnodes, vert_st, basis_coef, intrs_cpy
INTEGER									::	i, j, n_sub_pt, st, sub_region_ind									
INTEGER, DIMENSION(:), ALLOCATABLE		::	dind
INTEGER, DIMENSION(:,:), POINTER		::	t_sub_pt
!=======================================================================================
INTEGER                                           Gauss_point_number
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference_triangle

REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_local
INTEGER                                     ::	  pointer_reference_to_local(4)
REAL(8)                                     ::    matrix1(4,4)
!=========================================================================================
ALLOCATE(gwght(3))
gwght = (/1.0/3.0, 1.0/3.0, 1.0/3.0/)


CALL Inter_Tetra_Partition_2D(tri_sign, pointer_reference_to_local, el_type, vert, Intrs_pts, t_sub_pt, p_sub_pt)

n_sub_pt = SIZE(t_sub_pt,2)

CALL Gauss_Nodes_2D(p_sub_pt, t_sub_pt, g_sub_x, g_sub_y)


ALLOCATE(vert_st(2,3))
ALLOCATE(gnodes(2,3))
ALLOCATE(coeff_val(3))
ALLOCATE(dind(2))
ALLOCATE(bas_val_x(4,3), bas_val_y(4,3))

ALLOCATE(intrs_cpy(2,2))

intrs_cpy = Intrs_pts(:,1:2)
ALLOCATE(basis_coef(4,8))

intrs_cpy_1=TRANSPOSE(intrs_cpy)

CALL Linear_IFE_Basis_Coeff_2D( vert(:,pointer_reference_to_local), intrs_cpy, el_type, el_beta, basis_coef)

matrix = Zero

DO st=1,n_sub_pt

	vert_st = p_sub_pt(:,t_sub_pt(1:3,st))
	CALL Tetra_Volume_2D(vert_st, V)

	gnodes = RESHAPE((/g_sub_x(:,st),g_sub_y(:,st)/), (/2,3/), ORDER=(/2,1/))

	CALL Coeff_FUN(gnodes(1,:), gnodes(2,:), coeff_val, SIZE(gnodes, 2))

	sub_region_ind = t_sub_pt(4, st)

	DO i=1,n_nodes_in_elem_1    
		dind = (/1, 0/)
		CALL Linear_IFE_Basis_Eval_2D( basis_coef, gnodes(1,:), gnodes(2,:),		&
							           sub_region_ind, i, dind(1), dind(2), bas_val_x(i,:))

		dind = (/0, 1/)
		CALL Linear_IFE_Basis_Eval_2D( basis_coef, gnodes(1,:), gnodes(2,:),		&
							           sub_region_ind, i, dind(1), dind(2), bas_val_y(i,:))
	END DO

	DO i=1, n_nodes_in_elem_1						
   		DO j=1,i 

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

			matrix(i,j) = matrix(i,j) + el_beta(sub_region_ind)*V*eint

		END DO
	END DO

END DO


DO i=1,n_nodes_in_elem_1						
	DO j=i+1,n_nodes_in_elem_2    					
		matrix(i,j) = matrix(j,i)
	END DO
END DO

matrix1=0
DO  i=1, n_nodes_in_elem_1	
    DO j=1,n_nodes_in_elem_1	
    matrix1(pointer_reference_to_local(j),pointer_reference_to_local(i)) =  matrix(i,j)

    ENDDO
ENDDO

!==================delete teporate===================
matrix=0
matrix=matrix1
!======================================================


DEALLOCATE(gwght)
DEALLOCATE(t_sub_pt, p_sub_pt)
DEALLOCATE(g_sub_x, g_sub_y)
DEALLOCATE(vert_st)
DEALLOCATE(gnodes, basis_coef)
DEALLOCATE(coeff_val)
DEALLOCATE(dind)
DEALLOCATE(bas_val_x, bas_val_y)
DEALLOCATE(intrs_cpy)

END