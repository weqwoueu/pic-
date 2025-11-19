SUBROUTINE	generate_load_vector_local_IFE_on_interface(interface_coefficient_function_name, information_1, information_2,	&
														p_basic, element_index, t_c, nnx, nny,	&
														Gauss_coefficient_reference_1D, Gauss_point_reference_1D, &
														test_basis_type, r, node_type, node_index)

USE Object_2D
USE Domain_2D
USE IFE_INTERFACE, ONLY: Gauss_integration_local_load_IFE_on_interface, Gauss_integration_local_load_IFE_on_interface_Rectangular
IMPLICIT NONE


EXTERNAL											interface_coefficient_function_name
INTEGER, DIMENSION(:,:), INTENT(IN)				::  information_1
REAL(8), DIMENSION(:,:), INTENT(IN)				::  information_2
REAL(8), DIMENSION(:,:), INTENT(IN)				::	p_basic
INTEGER, DIMENSION(:), INTENT(IN)			   	::	element_index
INTEGER, DIMENSION(:,:), INTENT(IN)				::	t_c
INTEGER, INTENT(IN)                             ::  nnx,nny
REAL(8), DIMENSION(:), INTENT(IN)				::	Gauss_coefficient_reference_1D
REAL(8), DIMENSION(:), INTENT(IN)	            ::	Gauss_point_reference_1D
INTEGER, INTENT(IN)								::	test_basis_type
REAL(8), DIMENSION(:), POINTER					::	r, r_temp
INTEGER, DIMENSION(:,:), INTENT(IN)				::	node_type
INTEGEr, DIMENSION(:), INTENT(IN)				::	node_index


INTEGER										::	num_of_nodes, num_of_elements, num_of_unknowns
INTEGER										::	i, j
REAL(8)										::	vertices(2,4)
INTEGER										::	information_vector_1(18)
REAL(8)										::	information_vector_2(8)
INTEGER										::	alpha
REAL(8)										::	temp
REAL(8)										::	left_lower_point(2), h_partition(2)
INTEGER										::	count, count_x, count_y, object_num

INTEGER                                        ::  N_Boundary, i_boundary, Boundary_Shape, Boundary_Type
REAL(8)                                        ::  x_start, y_start, x_end, y_end

N_Boundary = SIZE(boundaries,1)
num_of_nodes = SIZE(p_basic,2)
num_of_elements = SIZE(t_c,2)
num_of_unknowns = 0

ALLOCATE(r_temp(num_of_nodes))

DO i=1, num_of_nodes
	r_temp(i) = 0.
ENDDO

DO i=1, num_of_elements
	IF(element_index(i) > 0)THEN					! INTERFACE ELEMENT
		DO j=1, 4
			vertices(:,j) = p_basic(:,t_c(j,i))
		ENDDO
		information_vector_1 = information_1(:,element_index(i))
		information_vector_2 = information_2(:,element_index(i))

		DO alpha=1, 4
			CALL Gauss_integration_local_load_IFE_on_interface(interface_coefficient_function_name, &
			                                                Gauss_coefficient_reference_1D, Gauss_point_reference_1D,	&
															vertices, information_vector_1, information_vector_2, &
															test_basis_type, alpha, temp, i)
			r_temp(t_c(alpha,i)) = r_temp(t_c(alpha,i)) + temp
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