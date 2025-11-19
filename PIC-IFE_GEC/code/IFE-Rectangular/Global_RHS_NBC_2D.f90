SUBROUTINE Global_RHS_NBC_2D_dj (element_index	,el_type1, information_1, information_2, p_basic, t_basic_int, e_basic,	&
					  	        p_int_x, p_int_y,  beta,	&
						        node_type, bnd_elem_index, vector, delta)	

! Neumann boundary conditions contribution to IFE RHS vector
! F_EBC = Int[Ud_j*Grad(Epsi_i).Grad(Epsi_j) dOmega]

USE IFE_MAIN_PARAM
USE IFE_Boundary
USE IFE_INTERFACE, ONLY: Inter_Tetra_Partition_2D, Gauss_Nodes_Edge_2D, &
                            Linear_FE_Basis_Eval_2D, Linear_IFE_Basis_Eval_2D, &
                            Linear_IFE_Basis_Coeff_2D, Check_Sub_Nodes_2D

IMPLICIT NONE

REAL(8), DIMENSION(:,:), POINTER	::	p_basic, p_int_x, p_int_y
INTEGER, DIMENSION(:,:), POINTER	::	t_basic_int, t_iel, node_type, e_basic
INTEGER, DIMENSION(:), POINTER		::	bnd_elem_index
REAL(8), DIMENSION(:), POINTER		::	vector, beta

INTEGER                                 delta

REAL(8)		vert(2,4), el_beta(2)

INTEGER		num_bound_elem, n_nodes_in_elem, m_unknowns,		&
			i, j, k, be, e, ie, el_type, el_region(2)      !, eindex     !, n_t, m_t, 

INTEGER		Nedge, n_nodes_in_edge
 
REAL(8), DIMENSION(:,:), POINTER	   ::	Intrs_pts
REAL(8), DIMENSION(:), ALLOCATABLE	   ::	gwght, BcEval       !, ibas_val !, eint
INTEGER, DIMENSION(:), POINTER		   ::	dind
REAL(8), DIMENSION(:,:), ALLOCATABLE   ::	gnodes, basis_coef, bas_val
                                            

REAL(8), DIMENSION(:,:), POINTER		::	p_sub_pt
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	vert_st, IFE_basis_coef, IFE_bas_val, intrs_cpy
INTEGER									::	n_sub_pt, st, IN1, IN2, sub_region_ind									
INTEGER, DIMENSION(:,:), POINTER		::	t_sub_pt
INTEGER, DIMENSION(:,:), ALLOCATABLE 	::	edges

!================================NEW ADD=========================================================
INTEGER, DIMENSION(:,:), POINTER             ::  information_1
REAL(8), DIMENSION(:,:), POINTER             ::  information_2
INTEGER, DIMENSION(:), POINTER		         ::  el_type1
INTEGER                                          pointer_reference_to_local_1(4), tri_sign(2)
INTEGER, DIMENSION(:), POINTER		         ::  element_index

!================================================================================================
REAL(8)   Norm_2, Length


ALLOCATE(gwght(3))
gwght = (/5.0/9.0, 8.0/9.0, 5.0/9.0/);

ALLOCATE(dind(2))
dind = (/0, 0/)

ALLOCATE(BcEval(3))

ALLOCATE(gnodes(2,3))
ALLOCATE(basis_coef(4,4))
ALLOCATE(bas_val(4, 3))
ALLOCATE(IFE_basis_coef(4,8))
ALLOCATE(IFE_bas_val(4, 3))
ALLOCATE(edges(3,2))
ALLOCATE(vert_st(2,3))
!=======================NEW=================================
edges = RESHAPE((/	1,	2,	&
	                2,	3,	&
	                3,	1	/), (/3,2/), ORDER=(/2,1/))
!edges = RESHAPE((/	1,	2,	&
!	                2,	3,	&
 !                   3,	4,  &
!	                4,  1	/), (/4,2/), ORDER=(/2,1/))
!=================================================================

num_bound_elem	= SIZE(bnd_elem_index)

n_nodes_in_elem = 4  ! Retanglar element
   
m_unknowns = MAXVAL(node_type)

n_nodes_in_edge = 2

ALLOCATE(vector(m_unknowns))
DO i =1, SIZE(vector)
	vector(i)	= Zero
END DO

DO be=1,num_bound_elem
	e = bnd_elem_index(be)

	DO i = 1, SIZE(e_basic,2)

!bc_index = 1: Direchlet, = 0: Neuuman

	   IF ( e == e_basic(3, i).AND.(bc_index(e_basic(4, i))==0 .AND. bc_index(e_basic(5, i))==0 )) THEN	

            IF (element_index(e) <= -1) THEN    ! Non-interface element

			   IF(bc_value(e_basic(4, i)) /= bc_value(e_basic(5, i))) THEN
			   PRINT*, 'bc_value of Node1=', bc_value(e_basic(4, i))
			   PRINT*, 'bc_value of Node2=', bc_value(e_basic(5, i))
			   PRINT*, 'bc_value of Node1 /= bc_value of Node2, STOP '
			   STOP
			   ENDIF

			   CALL Gauss_Nodes_Edge_2D(p_basic(:,e_basic(1, i)), p_basic(:,e_basic(2, i)), gnodes)

                vert = p_basic(:, t_basic_int(1:4,e)); 

                CALL Linear_FE_Basis_Coeff_2D(vert, basis_coef)

				DO j=1,n_nodes_in_elem
	               CALL Linear_FE_Basis_Eval_2D(basis_coef, gnodes(1,:), gnodes(2,:), j,		&
							                    dind(1), dind(2), bas_val(j,:))
  
				ENDDO
			    Length = Norm_2(p_basic(:,e_basic(1, i))-p_basic(:,e_basic(2, i)))

                 CALL Node_NBC_Eval_2D(BcEval, gnodes(1,:), gnodes(2,:), bc_value(e_basic(4, i)), SIZE(BcEval))

	            DO j=1,n_nodes_in_elem
		           IF (node_type(2,t_basic_int(j,e)) > 0 ) THEN
!	               ! This is an Unknown node
                      DO k = 1, n_nodes_in_edge
					     IF (e_basic(k, i) == t_basic_int(j,e)) THEN
						 ! This is an Neumann boundary conditions node
 							  IF (delta == 0) THEN
           		  					vector(node_type(2,t_basic_int(j,e))) = vector(node_type(2,t_basic_int(j,e))) &
													+0.5*Length*SUM(gwght*BcEval*bas_val(j,:))
							  ELSEIF(delta == 1) THEN
           		  					vector(node_type(2,t_basic_int(j,e))) = vector(node_type(2,t_basic_int(j,e))) &
													+0.5*Length*SUM(gwght*BcEval*gnodes(2,:)*bas_val(j,:))
							  ELSE 
								  PRINT*, ' delta value is wrong, check again, stop'
								  STOP
							  ENDIF
					      ENDIF
					  ENDDO  
                  ENDIF
	            END DO

   	        ELSEIF (element_index(e) > 0) THEN	! Interface element

!======================dukun add===========================================


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

!=============================dukun add======================================================

		      IF (el_type == 2 ) THEN

                 ALLOCATE(intrs_cpy(2,2))
                 intrs_cpy = Intrs_pts(:,1:2)

				 CALL Inter_Tetra_Partition_2D(tri_sign, pointer_reference_to_local_1, el_type1(ie), vert, Intrs_pts, t_sub_pt, p_sub_pt)

                  CALL Linear_IFE_Basis_Coeff_2D(vert, intrs_cpy, el_type, el_beta, IFE_basis_coef)

				 n_sub_pt = SIZE(t_sub_pt,2)

		         DO st=1,n_sub_pt
			    	vert_st = p_sub_pt(:,t_sub_pt(1:3,st))
	                sub_region_ind = t_sub_pt(4, st)

					DO Nedge = 1, 3

					   CALL Check_Sub_Nodes_2D(vert_st(:, edges(Nedge,1)),vert_st(:, edges(Nedge,2)),	  &
                                               p_basic(:,e_basic(1, i)), p_basic(:,e_basic(2, i)), IN1, IN2)

					   IF(IN1 /= -1 .AND. IN2 /= -1) THEN
					      IF(IN1 == 0 .AND. IN2 == 0) THEN
						  	 PRINT*, 'IN1 == 0 .AND. IN2 == 0, Global_RHS_NBC_2D. f90, STOP'  
						  	 STOP
						  ENDIF
			              CALL Gauss_Nodes_Edge_2D(vert_st(:, edges(Nedge,1)), vert_st(:, edges(Nedge,2)), gnodes)

				          DO j=1,n_nodes_in_elem
		                      CALL Linear_IFE_Basis_Eval_2D( IFE_basis_coef, gnodes(1,:), gnodes(2,:),		&
							                                 sub_region_ind, j, dind(1), dind(2), IFE_bas_val(j,:))

				          ENDDO

			              Length = Norm_2(vert_st(:, edges(Nedge,2))-vert_st(:, edges(Nedge,1)))

                          CALL Node_NBC_Eval_2D(BcEval, gnodes(1,:), gnodes(2,:), bc_value(e_basic(4, i)), SIZE(BcEval))

						  DO j=1,n_nodes_in_elem
		                    IF (node_type(2,t_basic_int(j,e)) > 0 ) THEN
          	                  ! This is an Unknown node
                               DO k = 1, n_nodes_in_edge
					             IF (e_basic(k, i) == t_basic_int(j,e)) THEN

 							       IF (delta == 0) THEN
           		  					  vector(node_type(2,e_basic(k, i))) = vector(node_type(2,e_basic(k, i))) &
													+0.5*Length*SUM(gwght*BcEval*IFE_bas_val(j,:))
							       ELSEIF(delta == 1) THEN
           		  					  vector(node_type(2,e_basic(k, i))) = vector(node_type(2,e_basic(k, i))) &
													+0.5*Length*SUM(gwght*BcEval*gnodes(2,:)*IFE_bas_val(j,:))
							       ELSE 
								      PRINT*, ' delta value is wrong, check again, stop'
								      STOP
							       ENDIF
							     ENDIF
						       ENDDO
							ENDIF
						  ENDDO

					   ENDIF  

					ENDDO

				 ENDDO

				 DEALLOCATE(Intrs_pts)
                 DEALLOCATE(intrs_cpy)
				 DEALLOCATE(p_sub_pt)
				 DEALLOCATE(t_sub_pt)


		      ELSEIF(el_type == 0) THEN

		      ENDIF
			ENDIF
	   ENDIF

	ENDDO

END DO

DEALLOCATE(gwght)
DEALLOCATE(dind)
DEALLOCATE(gnodes)
DEALLOCATE(basis_coef)
DEALLOCATE(bas_val)
DEALLOCATE(IFE_basis_coef)
DEALLOCATE(IFE_bas_val)
DEALLOCATE(edges)
DEALLOCATE(vert_st)
  
END


