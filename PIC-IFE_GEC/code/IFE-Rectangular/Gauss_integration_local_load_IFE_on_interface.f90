SUBROUTINE	Gauss_integration_local_load_IFE_on_interface(interface_coefficient_function_name, Gauss_coefficient_reference_1D, &
                                                Gauss_point_reference_1D,vertices, information_vector_1, &
                                                information_vector_2, test_basis_type, test_basis_index, r, element_Gauss)

USE Gauss_Data
USE IFE_INTERFACE, ONLY: generate_Gauss_local_1D
IMPLICIT NONE

EXTERNAL										interface_coefficient_function_name
REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_1D
REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_point_reference_1D
REAL(8), INTENT(IN)							::	vertices(2,4)
INTEGER, INTENT(IN)							::	information_vector_1(18)
REAL(8), INTENT(IN)							::	information_vector_2(8)
INTEGER, INTENT(IN)							::	test_basis_type
INTEGER, INTENT(IN)							::	test_basis_index
REAL(8)										::	r

REAL(8)										::	Dx, Dy, Ex, Ey
INTEGER										::	Gpn
REAL(8)										::	slope
REAL(8)										::	lower_bound, upper_bound
REAL(8), DIMENSION(:), POINTER				::	Gauss_coefficient_local_1D
REAL(8), DIMENSION(:), POINTER				::	Gauss_point_local_1D
INTEGER										::	i
REAL(8)										::	temp, temp1
REAL(8)										::	y_Gauss
!=========LY modification, 2022-7-25=========
Integer :: element_Gauss
!=========LY modification, 2022-7-25=========

Gpn = SIZE(Gauss_coefficient_reference_1D,1)
Dx = information_vector_2(3)
Dy = information_vector_2(4)
Ex = information_vector_2(5)
Ey = information_vector_2(6)
slope = (Dy - Ey) / (Dx - Ex)
r = 0.

IF(Dx == Ex)THEN
	lower_bound = MIN(Dy,Ey)
	upper_bound = MAX(Dy,Ey)
	CALL generate_Gauss_local_1D(Gauss_coefficient_reference_1D, Gauss_point_reference_1D,	&
							lower_bound, upper_bound, Gauss_coefficient_local_1D_Four, Gauss_point_local_1D_Four)

	DO i=1,Gpn
		CALL interface_coefficient_function_name(temp, Dx, Gauss_point_local_1D_Four(i))
		CALL retangular_local_basis_IFE(Dx, Gauss_point_local_1D_Four(i), vertices, information_vector_1, information_vector_2,	&
									1, test_basis_type, test_basis_index, 0, 0, temp1, element_Gauss)

		r = r + Gauss_coefficient_local_1D_Four(i) * temp * temp1
	ENDDO
ELSE
	lower_bound = MIN(Dx,Ex)
	upper_bound = MAX(Dx,Ex)
	CALL generate_Gauss_local_1D(Gauss_coefficient_reference_1D, Gauss_point_reference_1D,	&
							lower_bound, upper_bound, Gauss_coefficient_local_1D_Four, Gauss_point_local_1D_Four)

	DO i=1,Gpn
		y_Gauss = slope * (Gauss_point_local_1D_Four(i) - Dx) + Dy
		CALL interface_coefficient_function_name(temp, Gauss_point_local_1D_Four(i), y_Gauss)
		CALL retangular_local_basis_IFE(Gauss_point_local_1D_Four(i), y_Gauss, vertices, information_vector_1, information_vector_2,	&
									1, test_basis_type, test_basis_index, 0, 0, temp1, element_Gauss)


		r = r + Gauss_coefficient_local_1D_Four(i) * SQRT(1 + slope**2) * temp * temp1
	ENDDO

ENDIF

END SUBROUTINE



		

	