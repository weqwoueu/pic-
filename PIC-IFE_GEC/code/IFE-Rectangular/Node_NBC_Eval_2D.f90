SUBROUTINE Node_NBC_Eval_2D(U, xin, yin, Nf, n)

!USE Objects
USE IFE_Data

IMPLICIT NONE

INTEGER									n

REAL(8), DIMENSION(n), INTENT(IN)	::	xin, yin
INTEGER								::	region, i
REAL(8), DIMENSION(n), INTENT(OUT)	::	U
!INTEGER								::	Nf
REAL(8)                             ::  Nf

REAL(8)								beta_minus, beta_plus
REAL(8)                             R, x0, y0, RR
REAL(8)             				x(n), y(n)

beta_minus = Global_Beta(1);
IF (SIZE(Global_Beta)>1) beta_plus  = Global_Beta(2);
x0 = 0.6
y0 = 1.0
R  = 0.5

U = Nf


END