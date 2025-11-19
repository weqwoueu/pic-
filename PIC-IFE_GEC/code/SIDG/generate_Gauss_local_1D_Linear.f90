SUBROUTINE generate_Gauss_local_1D_Linear(Gauss_coefficient_reference_1D, Gauss_point_reference_1D, start_nodes, &
                                    end_nodes, Gauss_coefficient_local_1D_Linear, Gauss_point_local_1D_Linear)

IMPLICIT NONE

REAL(8), DIMENSION(:), INTENT(IN)		            ::	Gauss_coefficient_reference_1D
REAL(8), DIMENSION(:), INTENT(IN)                   ::	Gauss_point_reference_1D
REAL(8), DIMENSION(:), INTENT(IN)                   ::	start_nodes, end_nodes	
REAL(8), DIMENSION(:)								::	Gauss_coefficient_local_1D_Linear
REAL(8), DIMENSION(:, :)                            ::	Gauss_point_local_1D_Linear
REAL(8)                                             ::  upper_bound, lower_bound
INTEGER                                             ::  Gauss_Point_number

Gauss_Point_number = SIZE(Gauss_coefficient_reference_1D)

IF (start_nodes(1) == end_nodes(1)) THEN
    
    upper_bound = MAX(start_nodes(2),end_nodes(2))
    lower_bound = MIN(start_nodes(2),end_nodes(2))
    Gauss_coefficient_local_1D_Linear = (upper_bound-lower_bound)*Gauss_coefficient_reference_1D/2
    Gauss_point_local_1D_Linear(1,:) = start_nodes(1)
    Gauss_point_local_1D_Linear(2,:) = (upper_bound-lower_bound)*Gauss_point_reference_1D/2+(upper_bound+lower_bound)/2
    
ELSEIF (start_nodes(2) == end_nodes(2)) THEN
    
    upper_bound = MAX(start_nodes(1),end_nodes(1))
    lower_bound = MIN(start_nodes(1),end_nodes(1))
    Gauss_coefficient_local_1D_Linear = (upper_bound-lower_bound)*Gauss_coefficient_reference_1D/2
    Gauss_point_local_1D_Linear(1,:) = (upper_bound-lower_bound)*Gauss_point_reference_1D/2+(upper_bound+lower_bound)/2
    Gauss_point_local_1D_Linear(2,:) = start_nodes(2)
    
END IF



END SUBROUTINE
