SUBROUTINE	generate_load_vector_local_IFE_nonhomogeneous(flux_jump_function_name, Global_Beta, &
                                                        information_1, information_2, p_basic, element_index, t_c, nnx, nny, &
                                                        Gauss_coefficient_reference_1D, Gauss_point_reference_1D,	&
														Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle, &
														test_basis_type, test_derivative_degree_x, test_derivative_degree_y, &
														nonhomogeneous_trial_basis_type, nonhomogeneous_trial_derivative_degree_x, &
														nonhomogeneous_trial_derivative_degree_y, r, node_type)
															
USE IFE_INTERFACE, ONLY: Gauss_integration_of_flux_jump_on_interface, Gauss_integration_local_load_IFE_nonhomogeneous
USE IFE_MAIN_PARAM

IMPLICIT NONE

EXTERNAL										flux_jump_function_name
REAL(8), DIMENSION(:), INTENT(IN)        	::	Global_Beta
INTEGER, DIMENSION(:,:), INTENT(IN)         ::  information_1
REAL(8), DIMENSION(:,:), INTENT(IN)         ::  information_2
REAL(8), DIMENSION(:,:), INTENT(IN)	        ::	p_basic
INTEGER, DIMENSION(:), INTENT(IN)        	::	element_index
INTEGER, DIMENSION(:,:), INTENT(IN)			::	t_c
INTEGER, INTENT(IN)                         ::  nnx,nny
REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_1D
REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_point_reference_1D
REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_triangle
REAL(8), DIMENSION(:,:), INTENT(IN)         ::	Gauss_point_reference_triangle
INTEGER, INTENT(IN)                         ::  nonhomogeneous_trial_basis_type, nonhomogeneous_trial_derivative_degree_x, &
                                                    nonhomogeneous_trial_derivative_degree_y
INTEGER, INTENT(IN)							::	test_basis_type, test_derivative_degree_x, test_derivative_degree_y
REAL(8), DIMENSION(:), POINTER	            ::	r, r_temp
INTEGER, DIMENSION(:,:), INTENT(IN)			::	node_type

INTEGER										::	num_of_nodes, num_of_elements, num_of_unknowns
INTEGER										::	i, j
REAL(8)										::	vertices(2,4)
INTEGER										::	information_vector_1(18)
REAL(8)										::	information_vector_2(8)
REAL(8)										::	q, temp
INTEGER										::	alpha				




num_of_nodes = SIZE(p_basic,2)
num_of_elements = SIZE(t_c,2)
num_of_unknowns = 0

ALLOCATE(r_temp(num_of_nodes))

DO i=1, num_of_nodes
	r_temp(i) = 0.
ENDDO

DO i=1,num_of_elements
	IF(element_index(i) > 0)THEN			!	INTERFACE ELEMENT
		DO j=1,4
			vertices(:,j) = p_basic(:,t_c(j,i))
		ENDDO
		information_vector_1 = information_1(:,element_index(i))
		information_vector_2 = information_2(:,element_index(i))

		CALL Gauss_integration_of_flux_jump_on_interface(flux_jump_function_name, &
		                                                Gauss_coefficient_reference_1D, Gauss_point_reference_1D,	&
															vertices, information_vector_1, information_vector_2, q)
		
		DO alpha=1,4
			CALL Gauss_integration_local_load_IFE_nonhomogeneous(Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle,	&
																	vertices, information_vector_1, information_vector_2, test_basis_type,	&
																	alpha, test_derivative_degree_x, test_derivative_degree_y, nonhomogeneous_trial_basis_type,	&
																	nonhomogeneous_trial_derivative_degree_x, nonhomogeneous_trial_derivative_degree_y, temp, i)
			r_temp(t_c(alpha,i)) = r_temp(t_c(alpha,i)) + q * temp
		ENDDO
	ENDIF
ENDDO


DO i=1,num_of_nodes
	IF(node_type(2,i) > 0)THEN
		num_of_unknowns = num_of_unknowns + 1
	ENDIF
ENDDO

ALLOCATE(r(num_of_unknowns))

DO i=1,num_of_nodes
	IF(node_type(2,i) > 0)THEN
		r(node_type(2,i)) = r_temp(i)
	ENDIF
ENDDO

DEALLOCATE(r_temp)

END SUBROUTINE