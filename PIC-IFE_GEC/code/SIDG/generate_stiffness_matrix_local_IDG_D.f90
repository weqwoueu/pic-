SUBROUTINE generate_stiffness_matrix_local_IDG_D( Global_Beta, node_type,num_of_unknowns, &
                                                information_1, information_2, information_3_D, HP, element_index, &
                                                edge_index_D, HT, E_Dirichlet, HE, dimensions, &
                                                Gauss_coefficient_reference_1D, Gauss_point_reference_1D,&
                                                trial_basis_type, &
                                                test_basis_type, matrix, matrix_xt, delta, con_penalty)

USE IFE_MAIN_PARAM
USE Gauss_Data
USE IFE_INTERFACE, ONLY: generate_Gauss_local_1D_Linear, Sparse_Structure, Sparse_Structure_xt, &
                         gauss_integration_local_stiffness_DG_FE_D, gauss_integration_local_stiffness_IDG_D


IMPLICIT NONE
INTEGER,DIMENSION(:,:),POINTER          ::    information_1
REAL(8),DIMENSION(:,:),POINTER          ::    information_2, information_3_D
INTEGER,DIMENSION(:), POINTER           ::	  element_index, edge_index_D
REAL(8),DIMENSION(:)                    ::    Gauss_coefficient_reference_1D
REAL(8),DIMENSION(:)                    ::    Gauss_point_reference_1D
INTEGER, DIMENSION(:,:), POINTER        ::    HT, E_Dirichlet, HE
REAL(8)                                 ::    dimensions(2,2)
INTEGER                                 ::    out_index, in_index
INTEGER                                 ::    trial_basis_type, test_basis_type
REAL(8)                                 ::    con_penalty
    
INTEGER                                 ::    i, j ,k, m, n !, NZ
REAL(8),DIMENSION(:),POINTER            ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:),POINTER          ::    Gauss_point_local
REAL(8), DIMENSION(:,:), POINTER        ::	  HP 
INTEGER                                 ::    information_vector_1(18)
REAL(8)                                 ::    information_vector_2(8)
REAL(8)                                 ::    vertices_linear(2,2)
REAL(8)                                 ::    vertices_Rectangular(2,4)
REAL(8), DIMENSION(:), POINTER          ::	  Global_Beta
INTEGER,DIMENSION(4,2)                  ::    this_flag1, this_flag2

!-----------------------use the sparse matrix to store the Global stiff-----------------------
INTEGER, DIMENSION(:,:), POINTER        ::	VROW
TYPE(SPARSE), DIMENSION(:), POINTER     ::	matrix
TYPE(SPARSE), DIMENSION(:), ALLOCATABLE	::	matrix_xt
INTEGER, DIMENSION(:,:), POINTER	      ::	node_type
INTEGER                                 ::  NZ
INTEGER, INTENT(IN)                     ::	num_of_unknowns, delta
INTEGER                                 ::  n_nodes_in_elem, loc1, loc2
REAL(8)	                                ::  El_Stiff(4,4)
REAL(8)                                 ::  int_value

out_index = -1
in_index = -2
El_Stiff = 0    ! big problem 
n_nodes_in_elem = 4

CALL Sparse_Structure( HT, HE, node_type, VROW, matrix, NZ )
CALL Sparse_Structure_xt( HT, HE, node_type, VROW, matrix_xt, NZ )


do n=1,SIZE(E_Dirichlet,2)
    
  !LY Add for Boundary Condition, 2021-11-18
  IF (E_Dirichlet(6,n) == 1) THEN
        vertices_linear(:,1) = HP(:, E_Dirichlet(1,n))
        vertices_linear(:,2) = HP(:, E_Dirichlet(2,n))
        vertices_Rectangular(:, :) = HP(:,HT(1:4,E_Dirichlet(5,n)))
        this_flag1(1,:) = 0
        this_flag1(2:3,1) = E_Dirichlet(3:4,n)
        this_flag1(2:3,2) = E_Dirichlet(3:4,n)
        this_flag1(4,:) = 1
        
        El_Stiff = 0
        
    CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,vertices_linear(:,1),&
                                    vertices_linear(:,2), Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
        
        
        if (element_index(E_Dirichlet(5, n)) == out_index) then
            
            
            do j = 1, 4
              do k = 1, 4
                CALL gauss_integration_local_stiffness_DG_FE_D(delta, Global_Beta(1),Gauss_Coefficient_Local_1D_Linear_Eight,&
                                                        Gauss_Point_Local_1D_Linear_Eight, &
	                                                    this_flag1,vertices_Rectangular,trial_basis_type, j, &
											            test_basis_type, k, con_penalty, int_value, E_Dirichlet(5,n))
            
                    El_Stiff(j,k) = El_Stiff(j,k) + int_value
                enddo
            enddo
            
        elseif (element_index(E_Dirichlet(5, n)) == in_index) then
            
                        
            do j = 1, 4
              do k = 1, 4
                CALL gauss_integration_local_stiffness_DG_FE_D(delta, Global_Beta(2),Gauss_Coefficient_Local_1D_Linear_Eight,&
                                                        Gauss_Point_Local_1D_Linear_Eight, &
	                                                    this_flag1,vertices_Rectangular,trial_basis_type, j, &
											            test_basis_type, k, con_penalty, int_value, E_Dirichlet(5,n))
            
                    El_Stiff(j,k) = El_Stiff(j,k) + int_value
                enddo
            enddo
            
            
        elseif (element_index(E_Dirichlet(5, n)) > 0) then
            
            information_vector_1 = information_1(:, element_index(E_Dirichlet(5, n)))
            information_vector_2 = information_2(:, element_index(E_Dirichlet(5, n)))
            
            if (edge_index_D(n) < 0) then
                
                do j = 1, 4
                    do k = 1, 4
                        CALL gauss_integration_local_stiffness_IDG_D(delta, Global_Beta, information_vector_1, &
                                                information_vector_2, Gauss_Coefficient_Local_1D_Linear_Eight, &
                                                Gauss_Point_Local_1D_Linear_Eight, &
	                                            this_flag1, vertices_Rectangular, trial_basis_type, j, &
											    test_basis_type, k, edge_index_D(n), con_penalty, int_value, E_Dirichlet(5,n))
                
                        El_Stiff(j,k) = El_Stiff(j,k) + int_value
                    enddo
                enddo
                
            elseif (edge_index_D(n) > 0)  then
                
                CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                    vertices_linear(:,1), information_3_D(4:5, edge_index_D(n)), &
                                    Gauss_Coefficient_Local_1D_Linear_Eight, Gauss_Point_Local_1D_Linear_Eight)
                
                loc1 = information_3_D(6, edge_index_D(n))
                
                do j = 1, 4
                    do k = 1, 4
                        CALL gauss_integration_local_stiffness_IDG_D(delta, Global_Beta, information_vector_1, &
                                                information_vector_2, Gauss_Coefficient_Local_1D_Linear_Eight, &
                                                Gauss_Point_Local_1D_Linear_Eight, &
	                                            this_flag1, vertices_Rectangular, trial_basis_type, j, &
											    test_basis_type, k, loc1, con_penalty, int_value, E_Dirichlet(5,n))
                
                        El_Stiff(j,k) = El_Stiff(j,k) + int_value   ! unsure
                    enddo
                enddo
                
                CALL generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
                                    information_3_D(4:5, edge_index_D(n)), &
                                    vertices_linear(:,2),Gauss_Coefficient_Local_1D_Linear_Eight, &
                                    Gauss_Point_Local_1D_Linear_Eight)
                
                loc2 = information_3_D(7, edge_index_D(n))
                
                do j = 1, 4
                    do k = 1, 4
                        CALL gauss_integration_local_stiffness_IDG_D(delta, Global_Beta, information_vector_1, &
                                                information_vector_2, Gauss_Coefficient_Local_1D_Linear_Eight, &
                                                Gauss_Point_Local_1D_Linear_Eight, &
	                                            this_flag1, vertices_Rectangular, trial_basis_type, j, &
											    test_basis_type, k, loc2, con_penalty, int_value, E_Dirichlet(5,n))
                
                        El_Stiff(j,k) = El_Stiff(j,k) + int_value
                    enddo
                enddo
                
            else
                
                WRITE(6,*) 'error edge_index in generate_stiffness_matrix_local_IDG_D'
                STOP
                
            end if
            
        else
            
            WRITE(6,*) 'error element_index in generate_stiffness_matrix_local_IDG_D'
            STOP
            
        end if
        
            !=======================STORE MATRIX================================
            i = E_Dirichlet(5,n)

    	    DO k=1,n_nodes_in_elem  
		    IF ( node_type(2,HT(k,i)) > 0 ) THEN	! ˵���õ����ڲ��ĵ� iΪ��ǰ������Ԫ���
		    ! This is an Unknown node           
    		    DO j=1,n_nodes_in_elem
    			    IF ( node_type(2,HT(j,i)) > 0 ) THEN	
				    ! This is an Unknown node
                        IF (size(matrix%K)==1)  Then
                            m=matrix( node_type(2,HT(k,i)) )%SROW
					        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff(k,j)
							        EXIT
					        END IF
					    ElSE
					        DO m=matrix( node_type(2,HT(k,i)) )%SROW,  &
							        matrix( node_type(2,HT(k,i))+1 )%SROW-1
						        IF( matrix(m)%JCOL==node_type(2,HT(j,i)) )THEN
							        matrix(m)%K = matrix(m)%K+El_Stiff(k,j)
							        EXIT
						        END IF
					        END DO
                        ENDIF
                    ELSEIF( node_type(2,HT(j,i)) < 0 ) THEN
                        ! This is an known node��1
                        IF (size(matrix_xt%K)==1)  Then
                            m=matrix_xt( node_type(2,HT(k,i)) )%SROW
                            IF( matrix_xt(m)%JCOL==node_type(2,HT(j,i)) )THEN  ! ?????
                                matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff(k,j)
							    EXIT
                            END IF
                        ElSE
                            DO m=matrix_xt( node_type(2,HT(k,i)) )%SROW, &
                                matrix_xt( node_type(2,HT(k,i))+1 )%SROW-1
                                IF( matrix_xt(m)%JCOL==HT(j,i) )THEN
                                    matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff(k,j)
                                    EXIT
                                END IF
                            END DO
                        ENDIF
				    END IF
			    END DO
		    END IF
	    END DO

  END IF
    
enddo

end