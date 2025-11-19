SUBROUTINE Global_RHS_PDE_2D(	RHS_FUN, RHS_FUN1,U_full, R_full, p_basic, t_c1,		&
                                h_partition, node_type1, num_of_unknowns, RHS_HW, Interpol_HW, &
                                vector,delta,element_index,information_1,information_2)


! Partial differential equation contribution to IFE RHS vector
! F_PDE = Int[f(x,y,z,u) dOmega]
USE IFE_MAIN_PARAM
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Gauss_Nodes_Elem_2D, Generate_Gauss_reference, Generate_Gauss_local, &
                          Global_HardWire_2D, IFE_RHS_2D, FE_RHS_HW_2D

IMPLICIT NONE

EXTERNAL								RHS_FUN, RHS_FUN1
REAL(8), DIMENSION(:,:), POINTER	::	p_basic, p_int_x, p_int_y
REAL(8), DIMENSION(:),	INTENT(IN)	::	U_full, R_full
REAL(8), DIMENSION(:),	POINTER		::	beta
INTEGER, DIMENSION(:,:), POINTER	::	t_c1, t_iel, node_type1
INTEGER, INTENT(IN)					::	num_of_unknowns,delta
REAL(8), DIMENSION(:), POINTER		::	vector !,vector1
REAL(8)									RHS_HW(2,4,4), Interpol_HW(4,4)
REAL(8)                   			::	gnodes(2,4)

REAL(8)								::	vert(2,4), V, el_beta(2),		&
										El_RHS(4), U_full_El(4), R_full_El(4), &
										El_RHS_HW(4,4), El_Interpol_HW(4,4)
REAL(8), DIMENSION(:,:), POINTER	::	Intrs_pts
INTEGER								::	n_t, m_t, k, i, j, m,			&
										n_nodes_in_elem, m_unknowns ,	&
										e, ie, el_type, el_region(2), eindex

!================================dukun New add START=========================================
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
! Hardwire Global Arrays
REAL(8)		G_Stiff_HW(2,4,4), G_RHS_HW(2,4,4), G_Mass_HW(2,4,4,4), G_Interpol_HW(4,4)
!================================dukun New add END=========================================

!$ ========================== ab.ZWZ 2021/4/22 after changing the EBC ========================================\\
Gauss_point_number = 4
CALL Generate_Gauss_reference( Gauss_point_number, Gauss_coefficient_reference_Four, Gauss_point_reference_Four)
!$ ========================== ab.ZWZ 2021/4/22 after changing the EBC ========================================//
				
n_t	= SIZE(t_c1,2) ! n_t is the number of element

n_nodes_in_elem = 4 ! tangular element
   
m_unknowns = num_of_unknowns

ALLOCATE(vector(m_unknowns))! ,vector1(m_unknowns))
DO i =1, SIZE(vector)
	vector(i)	= Zero
END DO

!===========================================El_Interpol_HW======================================================================
CALL Global_HardWire_2D( FUN_2D_ONE, p_basic, t_c1,	element_index,			&
   					     G_Stiff_HW, G_RHS_HW, G_Mass_HW, Interpol_HW, delta)
!==================================================================================================================

DO e=1,n_t 

	El_RHS = Zero
	! here
    vert	= p_basic(:,t_c1(1:4,e));
    ! wsy add for sidg
    h_partition(1) = p_basic(1,t_c1(2,e)) - p_basic(1,t_c1(1,e))
    h_partition(2) = p_basic(2,t_c1(4,e)) - p_basic(2,t_c1(1,e))
    
	IF (element_index(e) < -1 .OR. element_index(e) == -1) THEN		! Non-interface element

		U_full_El = U_full(t_c1(1:4,e))
		R_full_El = R_full(t_c1(1:4,e))
		el_region = element_index(e)
		El_Interpol_HW	= Interpol_HW(:,:)


        CALL Generate_Gauss_local(Gauss_coefficient_reference_Four,Gauss_point_reference_Four, &
                             p_basic(1:2,t_c1(1,e)),h_partition, Gauss_coefficient_local_Four, Gauss_point_local_Four)

        gnodes=TRANSPOSE(Gauss_point_local_Four)

		IF   (delta == 0) THEN

		   CALL FE_RHS_HW_2D( vert, Gauss_coefficient_local_Four, gnodes, n_nodes_in_elem, El_RHS_HW, delta)

		ELSEIF(delta == 1) THEN
           vert = p_basic(:, t_c1(1:4,e));
		   CALL FE_RHS_HW_2D( vert, Gauss_coefficient_local_Four, gnodes, n_nodes_in_elem,El_RHS_HW, delta)
		ENDIF


		CALL FE_RHS_Eval_2D(RHS_FUN, U_full_El, R_full_El, el_region(1), &
					        n_nodes_in_elem, El_RHS_HW, El_Interpol_HW, El_RHS, gnodes, delta )   

   	ELSEIF (element_index(e) > 0) THEN	! Interface element
  
		ie		= element_index(e);
		

		el_region(information_1(15,ie)) = -1
        el_region(information_1(16,ie)) = element_index(e)

!================================new add=============================
			U_full_El =  U_full(t_c1(information_1(11:14,ie),e))
		    R_full_El =  R_full(t_c1(information_1(11:14,ie),e))
!====================================================================

        ALLOCATE(Intrs_pts(2,2))
        Intrs_pts(1,1) =information_2(3,ie)
        Intrs_pts(1,2) =information_2(4,ie)
        Intrs_pts(2,1) =information_2(5,ie)
		Intrs_pts(2,2) =information_2(6,ie)

		CALL IFE_RHS_2D(	RHS_FUN1, U_full_El, R_full_El, vert, Intrs_pts,			&
						    el_region, n_nodes_in_elem, El_RHS, delta, &
							information_1(:,element_index(e)), information_2(:,element_index(e))  )

        DEALLOCATE(Intrs_pts) !$ ab.ZWZ
	END IF
	DO j=1,n_nodes_in_elem
        IF ( node_type1(2,t_c1(j,e)) > 0 ) THEN	

			vector(node_type1(2,t_c1(j,e))) = vector(node_type1(2,t_c1(j,e))) + El_RHS(j)
            !print*,EL_RHS(j)
            !print*,j
            !pause

		END IF
	END DO

END DO


END