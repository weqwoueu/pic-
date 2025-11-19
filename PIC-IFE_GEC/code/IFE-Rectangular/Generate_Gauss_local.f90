SUBROUTINE   Generate_Gauss_local( Gauss_coefficient_reference, Gauss_point_reference, &
                                   left_lower_point, h_partition, Gauss_coefficient_local, Gauss_point_local)

!USE IFE_Data

IMPLICIT NONE


REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_reference
REAL(8),DIMENSION(:,:)		                ::    Gauss_point_reference
REAL(8)                           			::    left_lower_point(2)
REAL(8)	                                          h_partition(2)
REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:)				        ::    Gauss_point_local
REAL(8)	                                          Jacobi

Jacobi = h_partition(1)*h_partition(2)/4


Gauss_coefficient_local = Gauss_coefficient_reference*Jacobi

Gauss_point_local(:,1) = h_partition(1)*Gauss_point_reference(:,1)/2.+left_lower_point(1)+h_partition(1)/2.

Gauss_point_local(:,2) = h_partition(2)*Gauss_point_reference(:,2)/2.+left_lower_point(2)+h_partition(2)/2.



END