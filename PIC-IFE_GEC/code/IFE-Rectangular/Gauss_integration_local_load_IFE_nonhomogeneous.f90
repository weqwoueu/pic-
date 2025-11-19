SUBROUTINE	Gauss_integration_local_load_IFE_nonhomogeneous(Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle,	&
																	vertices, information_vector_1, information_vector_2, test_basis_type,	&
																	test_basis_index, test_derivative_degree_x, &
																	test_derivative_degree_y, nonhomogeneous_trial_basis_type,	&
																	nonhomogeneous_trial_derivative_degree_x, &
																	nonhomogeneous_trial_derivative_degree_y, r, element_Gauss)
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Generate_Gauss_local_triangle
IMPLICIT NONE

!INCLUDE	'Generate_Gauss_local_triangle.inc'

REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_triangle
REAL(8), DIMENSION(:,:), INTENT(IN)         ::	Gauss_point_reference_triangle
REAL(8), INTENT(IN)							::	vertices(2,4)
INTEGER, INTENT(IN)							::	information_vector_1(18)
REAL(8), INTENT(IN)							::	information_vector_2(8)
INTEGER, INTENT(IN)							::	test_basis_type, test_derivative_degree_x, test_derivative_degree_y
INTEGER, INTENT(IN)							::	test_basis_index
INTEGER, INTENT(IN)                         ::  nonhomogeneous_trial_basis_type, nonhomogeneous_trial_derivative_degree_x, &
                                                    nonhomogeneous_trial_derivative_degree_y
REAL(8)										::	r

INTEGER										::	Gpn
INTEGER										::	i
INTEGER										::	interface_element_type
INTEGER										::	pointer_reference_to_local(4)
REAL(8)										::	beta1, beta2
REAL(8)										::	Dx, Ex, Dy, Ey
REAL(8)										::	temp1, temp2, temp3, temp4
REAL(8)										::	r1, r2, r3, r4
REAL(8)										::	vertices_triangle(2,3)


REAL(8),DIMENSION(:),POINTER                ::  Gauss_coefficient_local_triangle
REAL(8),DIMENSION(:,:),POINTER              ::  Gauss_point_local_triangle
!=========LY modification, 2022-7-25=========
Integer :: element_Gauss
!=========LY modification, 2022-7-25=========


Gpn = SIZE(Gauss_coefficient_reference_triangle,1)
interface_element_type = information_vector_1(6)
pointer_reference_to_local = information_vector_1(11:14)
beta1 = information_vector_2(1)
beta2 = information_vector_2(2)
Dx = information_vector_2(3)
Dy = information_vector_2(4)
Ex = information_vector_2(5)
Ey = information_vector_2(6)


r = 0.

IF(interface_element_type == 1)THEN
	vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(1))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
	vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(1))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey

	CALL generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

	DO i=1,Gpn
		CALL retangular_local_basis_IFE_nonhomogeneous(Gauss_point_local_triangle_Nine(i,1), Gauss_point_local_triangle_Nine(i,2),	&
														vertices,information_vector_1,information_vector_2,1,	&
														nonhomogeneous_trial_basis_type,nonhomogeneous_trial_derivative_degree_x,&
														nonhomogeneous_trial_derivative_degree_y, temp1)
		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,1,test_basis_type,test_basis_index, &
										test_derivative_degree_x,test_derivative_degree_y,r1, element_Gauss)

		r = r + Gauss_coefficient_local_triangle_Nine(i) * beta1 * temp1 * r1
	ENDDO


	vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(4))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(4))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = vertices(2,pointer_reference_to_local(3))

	CALL generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

	DO i=1,Gpn
		CALL retangular_local_basis_IFE_nonhomogeneous(Gauss_point_local_triangle_Nine(i,1), Gauss_point_local_triangle_Nine(i,2),	&
														vertices,information_vector_1,information_vector_2,2,	&
														nonhomogeneous_trial_basis_type,nonhomogeneous_trial_derivative_degree_x,&
														nonhomogeneous_trial_derivative_degree_y, temp2)
		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
										test_derivative_degree_x,test_derivative_degree_y,r2, element_Gauss)

		r = r + Gauss_coefficient_local_triangle_Nine(i) * beta2 * temp2 * r2
	ENDDO
		

	vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
	vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(3))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey

	CALL generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

	DO i=1,Gpn
		CALL retangular_local_basis_IFE_nonhomogeneous(Gauss_point_local_triangle_Nine(i,1), Gauss_point_local_triangle_Nine(i,2),	&
														vertices,information_vector_1,information_vector_2,2,	&
														nonhomogeneous_trial_basis_type,nonhomogeneous_trial_derivative_degree_x,&
														nonhomogeneous_trial_derivative_degree_y, temp3)
		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
										test_derivative_degree_x,test_derivative_degree_y,r3, element_Gauss)

		r = r + Gauss_coefficient_local_triangle_Nine(i) * beta2 * temp3 * r3
	ENDDO

	vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(2))
	vertices_triangle(1,2) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(1,3) = Ex
	vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(2))
	vertices_triangle(2,2) = vertices(2,pointer_reference_to_local(3))
	vertices_triangle(2,3) = Ey

	CALL generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

	DO i=1,Gpn
		CALL retangular_local_basis_IFE_nonhomogeneous(Gauss_point_local_triangle_Nine(i,1), Gauss_point_local_triangle_Nine(i,2),	&
														vertices,information_vector_1,information_vector_2,2,	&
														nonhomogeneous_trial_basis_type,nonhomogeneous_trial_derivative_degree_x,&
														nonhomogeneous_trial_derivative_degree_y, temp4)
		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
										test_derivative_degree_x,test_derivative_degree_y,r4, element_Gauss)

		r = r + Gauss_coefficient_local_triangle_Nine(i) * beta2 * temp4 * r4
	ENDDO




ELSEIF(interface_element_type == 2)THEN


	vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(4))
	vertices_triangle(1,2) = vertices(1,pointer_reference_to_local(1))
	vertices_triangle(1,3) = Dx
	vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(4))
	vertices_triangle(2,2) = vertices(2,pointer_reference_to_local(1))
	vertices_triangle(2,3) = Dy

	CALL generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

	DO i=1,Gpn
		CALL retangular_local_basis_IFE_nonhomogeneous(Gauss_point_local_triangle_Nine(i,1), Gauss_point_local_triangle_Nine(i,2),	&
														vertices,information_vector_1,information_vector_2,1,	&
														nonhomogeneous_trial_basis_type,nonhomogeneous_trial_derivative_degree_x,&
														nonhomogeneous_trial_derivative_degree_y, temp1)
		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,1,test_basis_type,test_basis_index, &
										test_derivative_degree_x,test_derivative_degree_y,r1, element_Gauss)

		r = r + Gauss_coefficient_local_triangle_Nine(i) * beta1 * temp1 * r1
	ENDDO



	vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(1))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
	vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(1))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey

	CALL generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

	DO i=1,Gpn
		CALL retangular_local_basis_IFE_nonhomogeneous(Gauss_point_local_triangle_Nine(i,1), Gauss_point_local_triangle_Nine(i,2),	&
														vertices,information_vector_1,information_vector_2,1,	&
														nonhomogeneous_trial_basis_type,nonhomogeneous_trial_derivative_degree_x,&
														nonhomogeneous_trial_derivative_degree_y, temp2)
		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,1,test_basis_type,test_basis_index, &
										test_derivative_degree_x,test_derivative_degree_y,r2, element_Gauss)

		r = r + Gauss_coefficient_local_triangle_Nine(i) * beta1 * temp2 * r2
	ENDDO


	vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(2))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
	vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(2))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey

	CALL generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

	DO i=1,Gpn
		CALL retangular_local_basis_IFE_nonhomogeneous(Gauss_point_local_triangle_Nine(i,1), Gauss_point_local_triangle_Nine(i,2),	&
														vertices,information_vector_1,information_vector_2,2,	&
														nonhomogeneous_trial_basis_type,nonhomogeneous_trial_derivative_degree_x,&
														nonhomogeneous_trial_derivative_degree_y, temp3)
		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
										test_derivative_degree_x,test_derivative_degree_y,r3, element_Gauss)

		r = r + Gauss_coefficient_local_triangle_Nine(i) * beta2 * temp3 * r3
	ENDDO



	vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(1,2) = vertices(1,pointer_reference_to_local(2))
	vertices_triangle(1,3) = Dx
	vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(3))
	vertices_triangle(2,2) = vertices(2,pointer_reference_to_local(2))
	vertices_triangle(2,3) = Dy

	CALL generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

	DO i=1,Gpn
		CALL retangular_local_basis_IFE_nonhomogeneous(Gauss_point_local_triangle_Nine(i,1), Gauss_point_local_triangle_Nine(i,2),	&
														vertices,information_vector_1,information_vector_2,2,	&
														nonhomogeneous_trial_basis_type,nonhomogeneous_trial_derivative_degree_x,&
														nonhomogeneous_trial_derivative_degree_y, temp4)
		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
										test_derivative_degree_x,test_derivative_degree_y,r4, element_Gauss)

		r = r + Gauss_coefficient_local_triangle_Nine(i) * beta2 * temp4 * r4
	ENDDO

ENDIF


END SUBROUTINE