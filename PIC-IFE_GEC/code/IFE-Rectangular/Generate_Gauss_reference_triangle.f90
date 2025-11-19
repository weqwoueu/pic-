SUBROUTINE Generate_Gauss_reference_triangle( Gauss_point_number, &
                                            Gauss_coefficient_reference_triangle, &
                                            Gauss_point_reference_triangle)


IMPLICIT NONE
        
INTEGER												  Gauss_point_number
REAL(8),DIMENSION(:)							::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:)							::    Gauss_point_reference_triangle
REAL(8)											::    A1(4), A2(9), A3(3)
REAL(8)											::    B1(4,2),B2(9,2), B3(3,2)
INTEGER											::    m



A1(1) = (1.-1/(3**.5))/8.
A1(2) = (1.-1/(3**.5))/8.
A1(3) = (1.+1/(3**.5))/8.
A1(4) = (1.+1/(3**.5))/8.

A2(1) = 64./81.*(1-0)/8.
A2(2) = 100./324.*(1-(0.6**.5))/8.
A2(3) = 100./324.*(1-(0.6**.5))/8.
A2(4) = 100./324.*(1+(0.6**.5))/8.
A2(5) = 100./324.*(1+(0.6**.5))/8.
A2(6) = 40./81.*(1-0)/8.
A2(7) = 40./81.*(1-0)/8.
A2(8) = 40./81.*(1-(0.6**.5))/8.
A2(9) = 40./81.*(1+(0.6**.5))/8.

A3(1) = 1./6.
A3(2) = 1./6.
A3(3) = 1./6.

B1(1,1) = (1./(3.**0.5)+1.)/2.
B1(1,2 )= (1.-1./(3.**0.5))*(1./(3.**0.5)+1.)/4.
B1(2,1) = (1./(3.**0.5)+1.)/2.
B1(2,2 )= (1.-1./(3.**0.5))*(1.-1./(3.**0.5))/4.
B1(3,1) = (1.-1./(3.**0.5))/2.
B1(3,2 )= (1.+1./(3.**0.5))*(1./(3.**0.5)+1.)/4.
B1(4,1) = (1.-1./(3.**0.5))/2.
B1(4,2 )= (1.-1./(3.**0.5))*(1./(3.**0.5)+1.)/4.

B2(1,1)= (1.+0)/2.
B2(1,2)= (1.-0)*(1.+0)/4.
B2(2,1)= (1.+(0.6**.5))/2.
B2(2,2)= (1.+(0.6**.5))*(1.-(0.6**.5))/4.
B2(3,1)= (1.+(0.6**.5))/2.
B2(3,2)= (1.-(0.6**.5))*(1.-(0.6**.5))/4.
B2(4,1)= (1.-(0.6**.5))/2.
B2(4,2)= (1.+(0.6**.5))*(1.+(0.6**.5))/4.
B2(5,1)= (1.-(0.6**.5))/2.
B2(5,2)= (1.+(0.6**.5))*(1.-(0.6**.5))/4.
B2(6,1)= (1.+0)/2.
B2(6,2)= (1.+(0.6**.5))/4.
B2(7,1)= 1./2.
B2(7,2)= (1.-(0.6**.5))/4.
B2(8,1)= (1.+(0.6**.5))/2.
B2(8,2)= (1.-(0.6**.5))/4.
B2(9,1)= (1.-(0.6**.5))/2.
B2(9,2)= (1.+(0.6**.5))/4.

B3(1,1)= 0.5
B3(1,2)= 0
B3(2,1)= 0.5
B3(2,2)= 0.5
B3(3,1)= 0
B3(3,2)= 0.5

m = SIZE(Gauss_coefficient_reference_triangle)

IF (Gauss_point_number==4) THEN
	
    Gauss_coefficient_reference_triangle = A1

    Gauss_point_reference_triangle = B1

ELSEIF (Gauss_point_number==9) THEN
    
	Gauss_coefficient_reference_triangle = A2

	Gauss_point_reference_triangle = B2

ELSEIF (Gauss_point_number==3) THEN

	Gauss_coefficient_reference_triangle = A3
	
	Gauss_point_reference_triangle = B3


ENDIF




END