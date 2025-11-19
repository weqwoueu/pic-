SUBROUTINE Error_Analysis(beta, U_FUN_ErrAna, U, p_basic, t_basic_int,	&
                          h_partition, node_type, num_of_unknowns, delta, &
					      element_index,information_1,information_2, &
					      Error_2, Error_0, Error_1)

USE IFE_INTERFACE, ONLY: Generate_Gauss_reference, Generate_Gauss_local, Linear_FE_Basis_Eval_2D, Linear_FE_Basis_Coeff_2D, &
						 Inter_Tetra_Partition_2D, Gauss_Nodes_2D, Linear_IFE_Basis_Eval_2D, Linear_IFE_Basis_Coeff_2D

IMPLICIT NONE

REAL(8), DIMENSION(:,:), POINTER	::	p_basic, p_int_x, p_int_y !, p_int_z
REAL(8), DIMENSION(:),	POINTER		::	U, Error_2, Error_0
REAL(8), DIMENSION(:,:),	POINTER	::	Error_1
INTEGER, DIMENSION(:,:), POINTER	::	t_basic_int, t_iel

INTEGER                                           Gauss_point_number, n_nodes_in_elem, &
                                                  tri_sign(2), pointer_reference_to_local(4)
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_reference_triangle
REAL(8)	                                          h_partition(2), beta(2)
REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_local
INTEGER, DIMENSION(:), POINTER				::	  element_index
INTEGER,DIMENSION(:,:),POINTER              ::    information_1
REAL(8),DIMENSION(:,:),POINTER              ::    information_2

INTEGER, DIMENSION(:,:), POINTER	::	node_type
INTEGER, INTENT(IN)					::	num_of_unknowns,delta

REAL(8)									RHS_HW(2,4,4), Interpol_HW(4,4)

!=========================================================================================

INTERFACE
	SUBROUTINE U_FUN_ErrAna(x, y, region, dind_x, dind_y, f)
	REAL(8), DIMENSION(:), INTENT(IN)	::	x, y !, z
	INTEGER, INTENT(IN)					::	region, dind_x, dind_y !, dind_z
	REAL(8), DIMENSION(:), INTENT(OUT)	::	f
	END SUBROUTINE
END INTERFACE

REAL(8)							::	FE_coef(4,4), IFE_coef(4,8), vert_st(2,3)

REAL(8)							::	vert(2,4), V, el_beta(2),			&
									Intrs_pts(2,2), & !!vert_st(3,4),		&
									u_ex_val(4), u_ex_val1(3) ,u_h_val(4), u_h_val1(3), &
									u_h_el(4), u_h_el1(3), gnodes(2,4), gnodes1(2,3), &
									basis(4), basis1(3), gwght(4), gwght1(3), err_val(4), err_val1(3), eint,	&
									eintx, einty, eintz

INTEGER							::	n_t, i, j, k, dind(2),							&
									e, ie, el_type, el_region(2), st, n_sub_pt,		&
									sub_region_ind
INTEGER, DIMENSION(:,:),	POINTER		::	t_sub_pt
REAL(8), DIMENSION(:,:),	POINTER		::	p_sub_pt, g_sub_x, g_sub_y, g_sub_z

	

n_t	= SIZE(t_basic_int,2); ! n_t is the number of element
n_nodes_in_elem = size(t_basic_int,1) !the number of node in a element===dunkun

ALLOCATE(Error_2(n_t), Error_0(n_t), Error_1(n_t,3))
Error_2 = 0.D0
Error_0 = 0.D0
Error_1 = 0.D0

gnodes	= 0.D0
basis	= 0.D0
eint =0
eintx=0
einty=0

DO e=1,n_t 

	IF (element_index(e) < -1 .OR. element_index(e) == -1) THEN		! Non-interface element

		vert = p_basic(:, t_basic_int(:,e));

		CALL Linear_FE_Basis_Coeff_2D(vert, FE_coef)

!============================dukun add=======================================================================
   Gauss_point_number = 4

   CALL Generate_Gauss_reference( Gauss_point_number, Gauss_coefficient_reference, Gauss_point_reference)

   CALL Generate_Gauss_local(Gauss_coefficient_reference,Gauss_point_reference, &
                             p_basic(1:2,t_basic_int(1,e)),h_partition, Gauss_coefficient_local, Gauss_point_local)

   gwght=Gauss_coefficient_local
   gnodes=TRANSPOSE(Gauss_point_local)

!=============================================================================================================

		el_region   = element_index(e)

		u_h_el		= U(t_basic_int(1:4,e))

!!!!========================(0,0)==========================================================		
		dind = (/0, 0/)
		u_h_val		= 0.D0

        DO i=1,n_nodes_in_elem

	    CALL Linear_FE_Basis_Eval_2D(FE_coef, gnodes(1,:), gnodes(2,:), i,		&
							         dind(1), dind(2), basis)

        u_h_val = u_h_val + u_h_el(i)*basis

		ENDDO

		CALL U_FUN_ErrAna(	gnodes(1,:), gnodes(2,:), el_region(1),	&
							dind(1), dind(2), u_ex_val)

		err_val		= (DABS(u_h_val-u_ex_val))**2.

		eint = SUM(gwght*err_val);

!=============================(1,0)========================================================
		dind = (/1, 0/)
		u_h_val		= 0.D0
		DO i=1,n_nodes_in_elem
			CALL Linear_FE_Basis_Eval_2D(FE_coef, gnodes(1,:), gnodes(2,:), i,		&
							             dind(1), dind(2), basis)		
			u_h_val = u_h_val + u_h_el(i)*basis
		END DO

		CALL U_FUN_ErrAna(	gnodes(1,:), gnodes(2,:), el_region(1),	&
							dind(1), dind(2), u_ex_val)

		err_val		= (DABS(u_h_val-u_ex_val))**2.
		eintx = SUM(gwght*err_val);

!======================(0,1)=============================================================
		dind = (/0, 1/)
		u_h_val		= 0.D0
		DO i=1,n_nodes_in_elem
			CALL Linear_FE_Basis_Eval_2D(FE_coef, gnodes(1,:), gnodes(2,:), i,		&
							         dind(1), dind(2), basis)			
			u_h_val = u_h_val + u_h_el(i)*basis
		END DO
		CALL U_FUN_ErrAna(	gnodes(1,:), gnodes(2,:), el_region(1),	&
							dind(1), dind(2), u_ex_val)
		err_val		= (DABS(u_h_val-u_ex_val))**2.
		einty = SUM(gwght*err_val);
!=========================================================================================
				
   	ELSEIF (element_index(e) > 0) THEN	! Interface element 

        gwght1  = (/1.0/3.0, 1.0/3.0, 1.0/3.0/)
        ie		= (element_index(e)); 
        el_type = information_1(6,ie)

        el_region(1)    = -1;
		el_region(2)	= element_index(e);

		vert    = p_basic(:, t_basic_int(:,e));

		el_beta(1)= beta(information_1(15,ie))
   
        el_beta(2)= beta(information_1(16,ie))

		tri_sign(1)=information_1(15,ie)
		tri_sign(2)=information_1(16,ie)    !object in 2, vacumn in 1

		pointer_reference_to_local(1:4) = information_1(11:14,ie)

        Intrs_pts(1,1) =information_2(3,ie)
        Intrs_pts(1,2) =information_2(4,ie)
        Intrs_pts(2,1) =information_2(5,ie)
		Intrs_pts(2,2) =information_2(6,ie)

        CALL Linear_IFE_Basis_Coeff_2D( vert(:,pointer_reference_to_local), intrs_pts, el_type, el_beta, IFE_coef)	

        CALL Inter_Tetra_Partition_2D(tri_sign, pointer_reference_to_local, el_type, vert, Intrs_pts, t_sub_pt, p_sub_pt)

		n_sub_pt = SIZE(t_sub_pt,2)

		CALL Gauss_Nodes_2D(p_sub_pt, t_sub_pt, g_sub_x, g_sub_y)
		u_h_el		= U(t_basic_int(1:4,e))

		eint = 0.D0
		eintx = 0.D0
		einty = 0.D0

		DO st=1,n_sub_pt

			vert_st = p_sub_pt(:,t_sub_pt(1:3,st));

     	    CALL Tetra_Volume_2D(vert_st, V)

            gnodes1 = RESHAPE((/g_sub_x(:,st),g_sub_y(:,st)/), (/2,3/), ORDER=(/2,1/))

	        sub_region_ind = t_sub_pt(4, st)

            dind = (/0, 0/)
			u_h_val1		= 0.D0
			DO i=1,n_nodes_in_elem

                j=pointer_reference_to_local(i)

		        CALL Linear_IFE_Basis_Eval_2D( IFE_coef, gnodes1(1,:), gnodes1(2,:),		&
							                   sub_region_ind, i, dind(1), dind(2), basis1)
				u_h_val1 = u_h_val1 + u_h_el(j)*basis1

			END DO	

		    CALL U_FUN_ErrAna(	gnodes1(1,:), gnodes1(2,:), el_region(tri_sign(sub_region_ind)),	&
							    dind(1), dind(2), u_ex_val1)

			err_val1		= (DABS(u_h_val1-u_ex_val1))**2.
			eint = eint+V*SUM(gwght1*err_val1);	

            dind = (/1, 0/)
			u_h_val1		= 0.D0
			DO i=1,n_nodes_in_elem

                j=pointer_reference_to_local(i)
		        CALL Linear_IFE_Basis_Eval_2D( IFE_coef, gnodes1(1,:), gnodes1(2,:),		&
							                   sub_region_ind, i, dind(1), dind(2), basis1)
				u_h_val1 = u_h_val1 + u_h_el(j)*basis1
			END DO

		    CALL U_FUN_ErrAna(	gnodes1(1,:), gnodes1(2,:), el_region(tri_sign(sub_region_ind)),	&
							    dind(1), dind(2), u_ex_val1)

			err_val1		= (DABS(u_h_val1-u_ex_val1))**2.

			eintx = eintx+V*SUM(gwght1*err_val1);

            dind = (/0, 1/)
			u_h_val1		= 0.D0
			DO i=1,n_nodes_in_elem

                j=pointer_reference_to_local(i)
		        CALL Linear_IFE_Basis_Eval_2D( IFE_coef, gnodes1(1,:), gnodes1(2,:),		&
							                   sub_region_ind, i, dind(1), dind(2), basis1)
				u_h_val1 = u_h_val1 + u_h_el(j)*basis1
			END DO

		    CALL U_FUN_ErrAna(	gnodes1(1,:), gnodes1(2,:), el_region(tri_sign(sub_region_ind)),	&
							    dind(1), dind(2), u_ex_val1)
			err_val1		= (DABS(u_h_val1-u_ex_val1))**2.

			einty = einty+V*SUM(gwght1*err_val1);	
				
		END DO

		DEALLOCATE(t_sub_pt, p_sub_pt, g_sub_x, g_sub_y)
	
	END IF

		dind = (/0, 0/)
		DO i=1,4			

            IF (DSQRT(vert(1,i)**2+vert(2,i)**2)<=0.5) THEN
				CALL U_FUN_ErrAna(vert(1,i:i), vert(2,i:i), -2, &
								dind(1), dind(2), u_ex_val(i:i))

			ELSE
				CALL U_FUN_ErrAna(vert(1,i:i), vert(2,i:i),  -1, &
								dind(1), dind(2), u_ex_val(i:i))
			END IF
		END DO

		Error_0(e) = MAXVAL(DABS(u_h_el-u_ex_val))

		Error_2(e) = eint;
		Error_1(e,:) = (/eintx, einty/)

END DO

END