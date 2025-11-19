SUBROUTINE	generate_Gauss_reference_1D(Gauss_point_number,	Gauss_coefficient_reference_1D,	Gauss_point_reference_1D)

IMPLICIT NONE

INTEGER, INTENT(IN)						::	Gauss_point_number
REAL(8), DIMENSION(:)					::	Gauss_coefficient_reference_1D, Gauss_point_reference_1D
REAL(8)									::	A1(4), A2(8), A3(2)
REAL(8)									::	B1(4), B2(8), B3(2)


A1(1) = 0.3478548451
A1(2) = 0.3478548451
A1(3) = 0.6521451549
A1(4) = 0.6521451549

A2(1) = 0.1012285363
A2(2) = 0.1012285363
A2(3) = 0.2223810345
A2(4) = 0.2223810345
A2(5) = 0.3137066459
A2(6) = 0.3137066459
A2(7) = 0.3626837834
A2(8) = 0.3626837834

A3(1) = 1.
A3(2) = 1.

B1(1) = 0.8611363116
B1(2) = -0.8611363116
B1(3) = 0.3399810436
B1(4) = -0.3399810436

B2(1) = 0.9602898565
B2(2) = -0.9602898565
B2(3) = 0.7966664774
B2(4) = -0.7966664774
B2(5) = 0.5255324099
B2(6) = -0.5255324099
B2(7) = 0.1834346425
B2(8) = -0.1834346425

B3(1) = -1. / (3**0.5)
B3(2) = 1. / (3**0.5)

IF(Gauss_point_number == 4)THEN
	Gauss_coefficient_reference_1D = A1
	Gauss_point_reference_1D = B1
ELSEIF(Gauss_point_number == 8)THEN
	Gauss_coefficient_reference_1D = A2
	Gauss_point_reference_1D = B2
ELSEIF(Gauss_point_number == 2)THEN
	Gauss_coefficient_reference_1D = A3
	Gauss_point_reference_1D = B3
ENDIF

END	SUBROUTINE
	