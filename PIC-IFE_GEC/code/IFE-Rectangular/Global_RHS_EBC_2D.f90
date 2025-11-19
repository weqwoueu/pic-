SUBROUTINE Global_RHS_EBC_2D (	information_1,information_2, element_index, Coeff_FUN, U_full, p_basic, t_basic_int, e_basic,	&
							    el_type1, p_int_x, p_int_y, beta, node_type,	&
							    bnd_elem_index, Stiff_HW, vector, delta)
! Essential boundary conditions contribution to IFE RHS vector
! F_EBC = Int[Ud_j*Grad(Epsi_i).Grad(Epsi_j) dOmega]

USE IFE_MAIN_PARAM
USE IFE_Boundary
USE Gauss_Data
USE IFE_INTERFACE, ONLY: FE_Stiff_Eval_2D, FE_Stiff_HW_2D, IFE_Stiff_2D, Userdef_EBC_Value, Generate_Gauss_reference

IMPLICIT NONE

EXTERNAL								Coeff_FUN
INTEGER, DIMENSION(:), POINTER		::  element_index
INTEGER, DIMENSION(:), POINTER		::	bnd_elem_index, el_type1
INTEGER, DIMENSION(:,:), POINTER	::	t_basic_int, t_iel, node_type, e_basic, pointer_reference_to_local !, BETA_SIGN
REAL(8), DIMENSION(:), POINTER		::	U_full, vector, beta
REAL(8), DIMENSION(:,:), POINTER	::	p_basic , p_int_x, p_int_y
REAL(8)									Stiff_HW(2,4,4)
INTEGER                                 delta, pointer_reference_to_local_1(4), tri_sign(2)

REAL(8)		vert(2,4), V, el_beta(2), El_Stiff(4,4),	&							
			U_full_el(4), rhs_bc(4), El_Stiff_HW(4,4)
INTEGER		num_bound_elem, n_nodes_in_elem, m_unknowns,		&
			i, j, k, be, e, ie, el_type, el_region(2), eindex !, n_t, m_t, 
REAL(8), DIMENSION(:,:), POINTER			::	Intrs_pts
INTEGER, DIMENSION(:,:), POINTER            ::  information_1
REAL(8), DIMENSION(:,:), POINTER            ::  information_2

REAL(8), DIMENSION(:,:), POINTER			::	EBC_Value
INTEGER                                         Gauss_point_number



num_bound_elem	= SIZE(bnd_elem_index)

n_nodes_in_elem = 4 ! tangular element
   
m_unknowns = MAXVAL(node_type)

ALLOCATE(vector(m_unknowns))
DO i =1, SIZE(vector)
	vector(i)	= Zero
END DO

ALLOCATE(EBC_Value(n_nodes_in_elem,SIZE(t_basic_int,2)))
!CALL Userdef_EBC_Value(p_basic, t_basic_int, node_type, EBC_Value)

Gauss_point_number = 4

CALL Generate_Gauss_reference( Gauss_point_number, Gauss_coefficient_reference_Four, Gauss_point_reference_Four)

DO be=1,num_bound_elem

	e = bnd_elem_index(be)

	U_full_el = 0
	DO i = 1, SIZE(e_basic,2)

	     DO j=1,n_nodes_in_elem
!	set the initial value of U_full_el, change at here when debug different situation
		    IF((e_basic(1,i)==t_basic_int(j,e)).AND.(bc_index(e_basic(4,i))==1)) THEN
!2007.7.29, Boundary Condition Setup from object.inp
				!U_full_el(j) = EBC_Value(j,e)
				U_full_el(j) = bc_value(e_basic(4,i))

				U_full(t_basic_int(j,e)) = U_full_el(j)  
!!Yong Modified, Boundary Condition Setup from pic.inp

               IF (node_type(1,t_basic_int(j,e)) /= 1) THEN
                  PRINT*,'node_type(1,t_basic_int(j,e)) /= 1', e_basic(1,i), t_basic_int(j,e)
                  STOP
               ENDIF
			ENDIF
!	set the initial value of U_full_el, change at here when debug different situation
		    IF((e_basic(2,i)==t_basic_int(j,e)).AND.(bc_index(e_basic(5,i))==1)) THEN
!!2007.7.29, Boundary Condition Setup from object.inp
					!U_full_el(j) = EBC_Value(j,e)			   
					U_full_el(j) = bc_value(e_basic(5,i))

					U_full(t_basic_int(j,e)) = U_full_el(j)  

               IF (node_type(1,t_basic_int(j,e)) /= 1) THEN
                   PRINT*,'node_type(1,t_basic_int(j,e)) /= 1', e_basic(1,i), t_basic_int(j,e)
                   STOP
               ENDIF
			ENDIF  
		 ENDDO
	ENDDO


	El_Stiff = Zero

    IF (element_index(e) <= -1) THEN		! Non-interface element

        IF (element_index(e) == -1) THEN
            el_beta(1) = BETA(1)
		ELSEIF(element_index(e) < -1) THEN
            el_beta(1) = BETA(2)
		ENDIF
		el_region	= element_index(e)

		IF   (delta == 0) THEN
! delete temp by YC  		   El_Stiff_HW	= Stiff_HW(eindex,:,:)

           vert = p_basic(:, t_basic_int(1:4,e)); 

 		   CALL FE_Stiff_HW_2D( Coeff_FUN, vert, n_nodes_in_elem, n_nodes_in_elem, El_Stiff_HW, delta)
		ELSEIF(delta == 1) THEN
           vert = p_basic(:, t_basic_int(1:4,e));
		   CALL FE_Stiff_HW_2D( Coeff_FUN, vert, n_nodes_in_elem, n_nodes_in_elem, El_Stiff_HW, delta)

		ENDIF

	    CALL FE_Stiff_Eval_2D( el_beta(1), n_nodes_in_elem, n_nodes_in_elem, El_Stiff_HW, El_Stiff )

   	ELSEIF (element_index(e) > 0) THEN	! Interface element

		ie		= element_index(e);     
		vert	= p_basic(:,t_basic_int(1:4,e));

	    el_type = 2

		el_region	= element_index(e) !t_iel(6:7,ie);	! (T1 region, T2 region) = -1, -2, ...

!=====================WARNING THERE IS ONLY ONE BETA===========================	

        el_beta(1)= beta(information_1(15,ie))
        el_beta(2)= beta(information_1(16,ie))

		tri_sign(1)=information_1(15,ie)
		tri_sign(2)=information_1(16,ie)    !object in 2, vacumn in 1

        ALLOCATE(Intrs_pts(2,el_type))
        Intrs_pts(1,1) =information_2(3,ie)
        Intrs_pts(1,2) =information_2(4,ie)
        Intrs_pts(2,1) =information_2(5,ie)
		Intrs_pts(2,2) =information_2(6,ie)

		pointer_reference_to_local_1(1:4) = information_1(11:14,ie)

!==============================================================================

		IF (el_type == 2 ) THEN

		     CALL IFE_Stiff_2D(	tri_sign, pointer_reference_to_local_1, Coeff_FUN, vert, Intrs_pts, el_type1(ie), el_beta,	&
 			    	    	    n_nodes_in_elem, n_nodes_in_elem, El_Stiff, delta )
		     DEALLOCATE(Intrs_pts)
		ELSEIF(el_type == 0) THEN
		     PRINT*, 'el_type has some problem 1, check it, stop'
		     STOP
		ENDIF

	END IF

!========================NEW ADD========================================

    U_full_el = U_full(t_basic_int(1:4,e))

!======================================================================

	DO i=1,n_nodes_in_elem
		IF (node_type(2,t_basic_int(i,e)) > 0) THEN
		! This is an Unknown node
	         DO j=1,n_nodes_in_elem
	            IF (node_type(1,t_basic_int(j,e)) == 1) THEN	

					vector(node_type(2,t_basic_int(i,e))) = vector(node_type(2,t_basic_int(i,e))) &
													+ U_full_el(j)*El_Stiff(i,j)

				END IF
			END DO
		END IF
	END DO
   
END DO

END


