SUBROUTINE	Gauss_integration_of_flux_jump_on_interface(flux_jump_function_name, Gauss_coefficient_reference_1D, &
                                                        Gauss_point_reference_1D, vertices, information_vector_1, &
                                                        information_vector_2, r)

USE Gauss_Data
USE IFE_INTERFACE, ONLY: generate_Gauss_local_1D
IMPLICIT NONE

EXTERNAL										flux_jump_function_name
REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_coefficient_reference_1D
REAL(8), DIMENSION(:), INTENT(IN)           ::	Gauss_point_reference_1D
REAL(8), INTENT(IN)							::	vertices(2,4)
INTEGER, INTENT(IN)							::	information_vector_1(18)
REAL(8), INTENT(IN)							::	information_vector_2(8)
REAL(8)										::	r


INTEGER										::	Gpn, i, j
REAL(8)										::	Dx, Dy, Ex, Ey
REAL(8)										::	slope
REAL(8)										::	lower_bound, upper_bound
REAL(8), DIMENSION(:), POINTER				::	Gauss_coefficient_local_1D, Gauss_point_local_1D
REAL(8)										::	temp_flux_jump_function
REAL(8)										::	y_Gauss

Gpn = SIZE(Gauss_coefficient_reference_1D)
Dx = information_vector_2(3)
Dy = information_vector_2(4)
Ex = information_vector_2(5)
Ey = information_vector_2(6)
slope = (Dy - Ey) / (Dx - Ex)


r = 0.

IF(Dx==Ex)THEN
	lower_bound = MIN(Dy, Ey)
	upper_bound = MAX(Dy, Ey)
	
	CALL generate_Gauss_local_1D(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, lower_bound, &
	                            upper_bound, Gauss_coefficient_local_1D_Four, Gauss_point_local_1D_Four)
	
	DO i=1,Gpn
		CALL flux_jump_function_name(temp_flux_jump_function, Dx, Gauss_point_local_1D_Four(i))
		r = r + Gauss_coefficient_local_1D_Four(i) * temp_flux_jump_function
	ENDDO
ELSE
	lower_bound = MIN(Dx, Ex)
	upper_bound = MAX(Dx, Ex)

	CALL generate_Gauss_local_1D(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, lower_bound, &
	                             upper_bound, Gauss_coefficient_local_1D_Four, Gauss_point_local_1D_Four)
	DO i=1,Gpn
		y_Gauss = slope * (Gauss_point_local_1D_Four(i)-Dx)+Dy
		CALL flux_jump_function_name(temp_flux_jump_function, Gauss_point_local_1D_Four(i), y_Gauss)
		r = r + Gauss_coefficient_local_1D_Four(i) * SQRT(1 + slope**2) * temp_flux_jump_function
	ENDDO
ENDIF


END SUBROUTINE
