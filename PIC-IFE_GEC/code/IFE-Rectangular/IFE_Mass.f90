SUBROUTINE	IFE_Mass(	RHS_FUN, U_val, R_val, vert, Intrs_pts,					&
						el_region, n_nodes_in_elem_1, n_nodes_in_elem_2, matrix, information_vector_1, information_vector_2)


USE	IFE_MAIN_PARAM
USE IFE_INTERFACE, ONLY: Inter_Tetra_Partition_2D, Gauss_Nodes_2D, IFE_Interpol_2D, &
                            Linear_IFE_Basis_Eval_2D, Linear_IFE_Basis_Coeff_2D

IMPLICIT NONE

EXTERNAL				RHS_FUN
REAL(8), INTENT(IN)		::	vert(2,4) !, el_beta(2)
INTEGER, INTENT(IN)		::	n_nodes_in_elem_1, n_nodes_in_elem_2,	&
                            el_region(2)
REAL(8), INTENT(IN)		::	Intrs_pts(2, 2)
REAL(8), INTENT(IN)		::	U_val(n_nodes_in_elem_1), R_val(n_nodes_in_elem_1)
REAL(8), INTENT(OUT)	::	matrix(4,4)


REAL(8)									::	V, eint
REAL(8), DIMENSION(:), ALLOCATABLE		::	gwght, rhs_val, u_gnodes, r_gnodes
REAL(8), DIMENSION(:,:), POINTER		::	p_sub_pt, g_sub_x, g_sub_y, bas_val
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	gnodes, vert_st, basis_coef, intrs_cpy
INTEGER									::	i, j, k, n_sub_pt, st, sub_region_ind									
INTEGER, DIMENSION(:), ALLOCATABLE		::	dind
INTEGER, DIMENSION(:,:), POINTER		::	t_sub_pt
!=======================================================================================
INTEGER                                           information_vector_1(18) , tri_sign(2)
REAL(8)                                           information_vector_2(6) !, Intrs_pts_1(2,2),Intrs_pts(2,2)
INTEGER                                           pointer_reference_to_local(4)
REAL(8)                                           el_beta(2)
INTEGER                                           el_type


el_beta(1)=information_vector_2(1)
el_beta(2)=information_vector_2(2)
el_type	      =information_vector_1(6)
tri_sign(1)=information_vector_1(15)
tri_sign(2)=information_vector_1(16)
pointer_reference_to_local=information_vector_1(11:14)

ALLOCATE(gwght(3))
gwght = (/1.0/3.0, 1.0/3.0, 1.0/3.0/)


CALL Inter_Tetra_Partition_2D(tri_sign, pointer_reference_to_local, el_type, vert, Intrs_pts, t_sub_pt, p_sub_pt)
n_sub_pt = SIZE(t_sub_pt,2)

CALL Gauss_Nodes_2D(p_sub_pt, t_sub_pt, g_sub_x, g_sub_y)


ALLOCATE(bas_val(4,3))

ALLOCATE(vert_st(2,3))
ALLOCATE(gnodes(2,3))
ALLOCATE(rhs_val(3), u_gnodes(3), r_gnodes(3))
ALLOCATE(dind(2))


ALLOCATE(basis_coef(4,8))
ALLOCATE(intrs_cpy(2,2))

CALL Linear_IFE_Basis_Coeff_2D( vert(:,pointer_reference_to_local), Intrs_pts, el_type, el_beta, basis_coef)


matrix = 0
DO st=1,n_sub_pt

	vert_st = p_sub_pt(:,t_sub_pt(1:3,st))

	CALL Tetra_Volume_2D(vert_st, V)
	sub_region_ind = t_sub_pt(4,st)

	gnodes = RESHAPE(				&
			 (/g_sub_x(:,st),g_sub_y(:,st)/), (/2,3/), ORDER=(/2,1/) )

	dind = (/0, 0/)

	CALL IFE_Interpol_2D(	vert(:,pointer_reference_to_local), U_val, gnodes(1,:), gnodes(2,:), &
						Intrs_pts, el_type, el_beta, sub_region_ind, dind, u_gnodes)

	CALL IFE_Interpol_2D(	vert(:,pointer_reference_to_local), R_val, gnodes(1,:), gnodes(2,:), &
						Intrs_pts, el_type, el_beta, sub_region_ind, dind, r_gnodes)

    CALL RHS_FUN(r_gnodes, u_gnodes, el_region, rhs_val, SIZE(u_gnodes, 1))

	DO i=1,n_nodes_in_elem_1
		dind = (/0, 0/)    

!============================================
 !       j=pointer_reference_to_local(i)
!==============================================
		CALL Linear_IFE_Basis_Eval_2D(	basis_coef, gnodes(1,:), gnodes(2,:),		&
							            sub_region_ind, i, dind(1), dind(2), bas_val(i,:))
	END DO

	DO i=1,n_nodes_in_elem_1
!==============================================
		      k =pointer_reference_to_local(i)
!==============================================	
		
   		DO j=1,k !n_nodes_in_elem_2
!==============================================
!		      k =pointer_reference_to_local(i)
!==============================================

			eint = V*SUM(gwght*rhs_val*bas_val(k,:)*bas_val(j,:))

			matrix(k,j) = matrix(k,j) + eint		

		END DO
	END DO

END	DO

DO i=1,n_nodes_in_elem_1						
	DO j=i+1,n_nodes_in_elem_2    					
		matrix(i,j) = matrix(j,i)
	END DO
END DO


DEALLOCATE(gwght)
DEALLOCATE(t_sub_pt, p_sub_pt)
DEALLOCATE(g_sub_x, g_sub_y)
DEALLOCATE(vert_st)
DEALLOCATE(gnodes, basis_coef)
DEALLOCATE(rhs_val, u_gnodes, r_gnodes)
DEALLOCATE(dind)
DEALLOCATE(bas_val)
DEALLOCATE(intrs_cpy)

END