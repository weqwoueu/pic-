SUBROUTINE Generate_stiffness_matrix_local_IFE( minus_coefficient_function_name, Global_Beta, node_type,num_of_unknowns, &
                                                information_1, information_2, p_basic, element_index, &
									            t_c, nnx, nny, dimensions, &
												Gauss_coefficient_reference, &
												Gauss_point_reference, Gauss_coefficient_reference_triangle, &
												Gauss_point_reference_triangle, &
                                                trial_basis_type, &
												trial_derivative_degree_x,trial_derivative_degree_y, &
												test_basis_type,test_derivative_degree_x,test_derivative_degree_y, matrix, matrix_xt, delta)

USE IFE_MAIN_PARAM
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Generate_Gauss_local, Gauss_integration_local_stiffness, &
                            Gauss_integration_local_stiffness_IFE, Sparse_Structure, Sparse_Structure_xt


IMPLICIT NONE

INTEGER,DIMENSION(:,:),POINTER              ::    information_1
REAL(8),DIMENSION(:,:),POINTER              ::    information_2
INTEGER,DIMENSION(:), POINTER	        	::	  element_index
REAL(8),DIMENSION(:)						::    Gauss_coefficient_reference
REAL(8),DIMENSION(:,:)						::    Gauss_point_reference
REAL(8),DIMENSION(:)					    ::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:)			            ::    Gauss_point_reference_triangle
REAL(8),DIMENSION(:,:),POINTER              ::    r
INTEGER, DIMENSION(:,:), POINTER			::    t_c
REAL(8)                                           minus_coefficient_function_name,plus_coefficient_function_name
REAL(8)	                                          dimensions(2,2)
INTEGER                                           nnx,nny,num_of_nodes
INTEGER                                           trial_basis_type, trial_derivative_degree_x,trial_derivative_degree_y, &
												  test_basis_type,test_derivative_degree_x,test_derivative_degree_y

INTEGER                                           i, j ,k,m !, NZ
REAL(8)								        	  x_s, x_l, y_s, y_l, h_x, h_y,temp
REAL(8)	                                          h_partition(2)
INTEGER								        	  n_x, n_y

REAL(8),DIMENSION(:),POINTER                ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:),POINTER              ::    Gauss_point_local
REAL(8), DIMENSION(:,:), POINTER	        ::	  p_basic
INTEGER                                           information_vector_1(18)
REAL(8)                                           information_vector_2(8)
REAL(8)                                           vertices(2,4)
REAL(8), DIMENSION(:), POINTER	        	::	  Global_Beta
INTEGER                                           a !, b


!-----------------------use the sparse matrix to store the Global stiff-----------------------
INTEGER, DIMENSION(:,:), POINTER	::	VROW			
TYPE(SPARSE), DIMENSION(:), POINTER	::	matrix
TYPE(SPARSE), DIMENSION(:), ALLOCATABLE	::	matrix_xt   !$ ab.ZWZ 2021/7/9
INTEGER, DIMENSION(:,:), POINTER	::	node_type
INTEGER                                 NZ
INTEGER, INTENT(IN)					::	num_of_unknowns, delta
INTEGER                                 n_nodes_in_elem
REAL(8)	                            ::  El_Stiff(4,4)

n_nodes_in_elem = 4
!CALL Sparse_Structure( t_c, node_type, VROW, matrix, NZ )
!CALL Sparse_Structure_xt( t_c, node_type, VROW, matrix_xt, NZ ) !$ ab.ZWZ 2021/7/9

!--------------------------------------------------------------------------------------------
!k = 0
!matrix(1)%SROW = 1
!DO i=1,num_of_unknowns
!	DO j=2,VROW(i,1)+1
!		k = k+1
!		matrix(k)%JCOL = VROW(i,j)		
!	END DO
!
!	matrix(i+1)%SROW = matrix(i)%SROW+VROW(i,1)	
!
!END DO
!
!DEALLOCATE(VROW)


!-----------------------use the sparse matrix to store the Global stiff-----------------------

num_of_nodes	= SIZE(p_basic,2)

x_s = dimensions(1,1)
x_l = dimensions(1,2)
y_s = dimensions(2,1)
y_l = dimensions(2,2)
                 
n_x = nnx
n_y = nny

h_x = (x_l - x_s)/(n_x - 1)
h_y = (y_l - y_s)/(n_y - 1)

h_partition(1) = h_x
h_partition(2) = h_y


DO i=1,SIZE(t_c,2)
  
    ! wsy add for sidg
    h_partition(1) = p_basic(1,t_c(2,i)) - p_basic(1,t_c(1,i))
    h_partition(2) = p_basic(2,t_c(4,i)) - p_basic(2,t_c(1,i))
    
   CALL Generate_Gauss_local(Gauss_coefficient_reference,Gauss_point_reference, &
                             p_basic(1:2,t_c(1,i)),h_partition, Gauss_coefficient_local_Nine, Gauss_point_local_Nine)

   IF (element_index(i) == -1) THEN  !OUTSIDE THE OBJECTS
   
   DO j=1, 4
      
	  DO k=1, 4

      CALL Gauss_integration_local_stiffness(delta, minus_coefficient_function_name,Gauss_coefficient_local_Nine, &
											 Gauss_point_local_Nine, p_basic(1:2,t_c(1,i)), &
											 h_partition,trial_basis_type, j, &
											 trial_derivative_degree_x,trial_derivative_degree_y, test_basis_type, k, &
											 test_derivative_degree_x,test_derivative_degree_y, temp, i);

       El_Stiff(j,k) = temp

	  ENDDO

   ENDDO

   ELSEIF (element_index(i) < -1) THEN  !INSIDE THE OBJECTS

        plus_coefficient_function_name = Global_Beta(ABS(element_index(i)))

        DO j=1, 4
      
            DO k=1, 4

                CALL Gauss_integration_local_stiffness(delta, plus_coefficient_function_name, &
											         Gauss_coefficient_local_Nine, Gauss_point_local_Nine, &
	                                                 p_basic(1:2,t_c(1,i)),h_partition,trial_basis_type, j, &
											         trial_derivative_degree_x,trial_derivative_degree_y, test_basis_type, k, &
											         test_derivative_degree_x,test_derivative_degree_y, temp, i);



                 El_Stiff(j,k) = temp

	        ENDDO
	
	    ENDDO
 
   ELSEIF (element_index(i) > 0) THEN  !THE INTERSECT ELEMENTS

       vertices = p_basic(1:2,t_c(1:4,i))
       information_vector_1 = information_1(:,element_index(i))
       information_vector_2 = information_2(:,element_index(i)) 
   
       DO j= 1, 4

           DO k= 1, 4
     
             CALL Gauss_integration_local_stiffness_IFE(delta, Gauss_coefficient_reference_triangle, &
													    Gauss_point_reference_triangle, &
			                                            vertices,information_vector_1,information_vector_2,trial_basis_type, &
												        j,trial_derivative_degree_x,trial_derivative_degree_y,test_basis_type, &
													    k,test_derivative_degree_x,test_derivative_degree_y,temp, i)

             El_Stiff(j,k) = temp
      
	      ENDDO
      
       ENDDO 

   ENDIF
   !IF(element_index(i)>0)THEN
   !     print*,i
   !     print*,element_index(i)
   !     print*,El_stiff
   !     pause
   !endif
!=======================STORE MATRIX================================

   	DO k=1,n_nodes_in_elem
		IF ( node_type(2,t_c(k,i)) > 0 ) THEN	
		! This is an Unknown node
    		DO j=1,n_nodes_in_elem
    			IF ( node_type(2,t_c(j,i)) > 0 ) THEN	
				! This is an Unknown node
                    IF (size(matrix%K)==1)  Then
                        m=matrix( node_type(2,t_c(k,i)) )%SROW
					IF( matrix(m)%JCOL==node_type(2,t_c(j,i)) )THEN
						matrix(m)%K = matrix(m)%K+El_Stiff(k,j) 
						EXIT
					END IF
					ElSE
					    DO m=matrix( node_type(2,t_c(k,i)) )%SROW,					&
							    matrix( node_type(2,t_c(k,i))+1 )%SROW-1
						    IF( matrix(m)%JCOL==node_type(2,t_c(j,i)) )THEN
							    matrix(m)%K = matrix(m)%K+El_Stiff(k,j) 
							    EXIT
						    END IF
					    END DO
                    ENDIF
                !$ ============= ab.ZWZ 2021/7/9 ==================== \\   
                ELSEIF( node_type(2,t_c(j,i)) < 0 ) THEN
                    ! This is an Known node
                    IF (size(matrix_xt%K)==1)  Then
                        m=matrix_xt( node_type(2,t_c(k,i)) )%SROW
                        IF( matrix_xt(m)%JCOL==node_type(2,t_c(j,i)) )THEN
                            matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff(k,j)
							EXIT
                        END IF
                    ElSE
                        DO m=matrix_xt( node_type(2,t_c(k,i)) )%SROW, &
                            matrix_xt( node_type(2,t_c(k,i))+1 )%SROW-1
                            IF( matrix_xt(m)%JCOL==t_c(j,i) )THEN
                                matrix_xt(m)%K = matrix_xt(m)%K+El_Stiff(k,j)
                                EXIT
                            END IF
                        END DO
                    ENDIF
                !$ ============= ab.ZWZ 2021/7/9 ==================== //  
				END IF
			END DO
		END IF
	END DO


ENDDO

!=======================STORE MATRIX================================


END