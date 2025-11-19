SUBROUTINE	Global_Mass( RHS_FUN, U_full, R_full, p_basic, t_basic_int, &
    					 h_partition, node_type, num_of_unknowns, Mass_HW, Interpol_HW, matrix, element_index,information_1,information_2)
USE IFE_MAIN_PARAM
USE IFE_INTERFACE, ONLY: Gauss_Nodes_Elem_2D, Global_HardWire_2D, FE_Mass, IFE_Mass

IMPLICIT NONE

EXTERNAL								RHS_FUN
REAL(8), DIMENSION(:,:), POINTER	::	p_basic, p_int_x, p_int_y
REAL(8), DIMENSION(:),	INTENT(IN)	::	U_full, R_full
REAL(8), DIMENSION(:),	POINTER		::	beta
INTEGER, DIMENSION(:,:), POINTER	::	t_basic_int, t_iel, node_type
INTEGER, INTENT(IN)					::	num_of_unknowns
TYPE(SPARSE), DIMENSION(:), POINTER	::	matrix
REAL(8)								::	Mass_HW(4,4,4), Interpol_HW(4,4)

REAL(8)								::	Vert(2,4), V, el_beta(2),					&
										El_Mass(4,4), U_full_El(4), R_full_El(4),	&
										El_Mass_HW(4,4,4), El_Interpol_HW(4,4)
REAL(8), DIMENSION(:,:), POINTER	::	Intrs_pts
INTEGER								::	n_t, m_t, k, i, j, m,			&
										n_nodes_in_elem, m_unknowns, 	&
										e, ie, el_type, el_region(2) !, eindex
!REAL(8)                          	::	gnodes(2,3)
!=======================================================================================
INTEGER                                           Gauss_point_number
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference_triangle
REAL(8)	                                          h_partition(2)
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_local
INTEGER, DIMENSION(:), POINTER				::	  element_index
INTEGER,DIMENSION(:,:),POINTER              ::    information_1
REAL(8),DIMENSION(:,:),POINTER              ::    information_2
EXTERNAL									      FUN_2D_ONE
REAL(8)		                                      G_Stiff_HW1(2,4,4), G_RHS_HW1(2,4,4), G_Mass_HW1(2,4,4,4)!, Interpol_HW(4,4)!, G_Interpol_HW(4,4)
!======================================================================================

n_t = SIZE(t_basic_int,2)				!the number of elements


n_nodes_in_elem = 4					!retangular element

m_unknowns = num_of_unknowns

DO i =1, SIZE(matrix)
	matrix(i)%K	= Zero
END DO

!===========================================El_Interpol_HW=========================================================
CALL Global_HardWire_2D( FUN_2D_ONE, p_basic, t_basic_int,	element_index,			&
   					     G_Stiff_HW1, G_RHS_HW1, G_Mass_HW1, Interpol_HW, 0)
!==================================================================================================================

DO e=1,n_t
	El_Mass = Zero

	IF (element_index(e) < -1 .OR. element_index(e) == -1) THEN		! Non-interface element

		el_region = element_index(e)

		U_full_El = U_full(t_basic_int(1:4,e))
		R_full_El = R_full(t_basic_int(1:4,e))

        El_Mass_HW = Mass_HW(:,:,:)

		El_Interpol_HW = Interpol_HW(:,:)
		CALL	FE_Mass(RHS_FUN, U_full_El, R_full_El, el_region(1), &
					    n_nodes_in_elem, n_nodes_in_elem, El_Mass_HW, El_Interpol_HW, El_Mass)

	ELSEIF	( t_basic_int(4,e)>0 )	THEN	!interface element

		ie		= element_index(e)

		vert	= p_basic(:,t_basic_int(1:4,e));

		el_region(information_1(15,ie)) = -1
        el_region(information_1(16,ie)) = element_index(e)

        ALLOCATE(Intrs_pts(2,2))
        Intrs_pts(1,1) =information_2(3,ie)
        Intrs_pts(1,2) =information_2(4,ie)
        Intrs_pts(2,1) =information_2(5,ie)
		Intrs_pts(2,2) =information_2(6,ie)
    	U_full_El =  U_full(t_basic_int(information_1(11:14,ie),e))
		R_full_El =  R_full(t_basic_int(information_1(11:14,ie),e))

		CALL IFE_Mass(	RHS_FUN, U_full_El, R_full_El, vert, Intrs_pts,					&
                        el_region, n_nodes_in_elem, n_nodes_in_elem, El_Mass, &
						information_1(:,element_index(e)), information_2(:,element_index(e)))

		DEALLOCATE(Intrs_pts)
	ENDIF

	DO i=1,n_nodes_in_elem
		IF ( node_type(2,t_basic_int(i,e)) > 0 ) THEN	
		! This is an Unknown node
    		DO j=1,n_nodes_in_elem
    			IF ( node_type(2,t_basic_int(j,e)) > 0 ) THEN	
				! This is an Unknown node
					DO m=matrix( node_type(2,t_basic_int(i,e)) )%SROW,					&
							matrix( node_type(2,t_basic_int(i,e))+1 )%SROW-1
						IF( matrix(m)%JCOL==node_type(2,t_basic_int(j,e)) )THEN
							matrix(m)%K = matrix(m)%K+El_Mass(i,j)
							EXIT
						END IF
					END DO

				END IF
			END DO
		END IF
	END DO

END DO

END