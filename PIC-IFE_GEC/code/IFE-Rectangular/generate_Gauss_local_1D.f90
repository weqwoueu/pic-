SUBROUTINE generate_Gauss_local_1D(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, lower_bound, &
                                    upper_bound, Gauss_coefficient_local_1D, Gauss_point_local_1D)

IMPLICIT NONE

REAL(8), DIMENSION(:), INTENT(IN)		            ::	Gauss_coefficient_reference_1D
REAL(8), DIMENSION(:), INTENT(IN)		            ::	Gauss_point_reference_1D
REAL(8), INTENT(IN)									::	lower_bound, upper_bound	
REAL(8), DIMENSION(:)								::	Gauss_coefficient_local_1D
REAL(8), DIMENSION(:)								::	Gauss_point_local_1D


Gauss_coefficient_local_1D = (upper_bound-lower_bound) * Gauss_coefficient_reference_1D/2


Gauss_point_local_1D = (upper_bound-lower_bound) * Gauss_point_reference_1D/2+(upper_bound+lower_bound)/2



END SUBROUTINE
