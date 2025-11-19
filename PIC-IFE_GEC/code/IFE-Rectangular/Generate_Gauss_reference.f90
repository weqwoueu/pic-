SUBROUTINE Generate_Gauss_reference( Gauss_point_number, Gauss_coefficient_reference, Gauss_point_reference)


IMPLICIT NONE
        
INTEGER                                              Gauss_point_number
REAL(8),DIMENSION(:)		                   ::    Gauss_coefficient_reference
REAL(8),DIMENSION(:,:)			               ::    Gauss_point_reference
REAL(8)                                        ::    A1(4), A2(9)
REAL(8)                                        ::    B1(4,2),B2(9,2)
INTEGER										   ::	 m

DATA A1  / 1, 1, 1, 1/
A2(1) = 64./81.
A2(2) = 100./324.
A2(3) = 100./324.
A2(4) = 100./324.
A2(5) = 100./324.
A2(6) = 40./81.
A2(7) = 40./81.
A2(8) = 40./81.
A2(9) = 40./81.


B1(1,1) = 1./(3.**0.5)
B1(1,2 )= 1./(3.**0.5)
B1(2,1) = 1./(3.**0.5)
B1(2,2 )= -1./(3.**0.5)
B1(3,1) = -1./(3.**0.5)
B1(3,2 )= 1./(3.**0.5)
B1(4,1) = -1./(3.**0.5)
B1(4,2 )= -1./(3.**0.5)

B2(1,1)= 0
B2(1,2)= 0
B2(2,1)= (0.6)**0.5
B2(2,2)= (0.6)**0.5
B2(3,1)= (0.6)**0.5
B2(3,2)= -(0.6)**0.5
B2(4,1)= -(0.6)**0.5
B2(4,2)= (0.6)**0.5
B2(5,1)= -(0.6)**0.5
B2(5,2)= -(0.6)**0.5
B2(6,1)= 0
B2(6,2)= (0.6)**0.5
B2(7,1)= 0
B2(7,2)= -(0.6)**0.5
B2(8,1)= (0.6)**0.5
B2(8,2)= 0
B2(9,1)= -(0.6)**0.5
B2(9,2)= 0


m = SIZE(Gauss_coefficient_reference)

IF (Gauss_point_number==4) THEN
	
    Gauss_coefficient_reference = A1
    Gauss_point_reference = B1

ELSEIF (Gauss_point_number==9) THEN
	
	Gauss_coefficient_reference = A2
	Gauss_point_reference = B2

ENDIF




END