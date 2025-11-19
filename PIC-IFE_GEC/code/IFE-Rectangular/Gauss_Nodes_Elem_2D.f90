SUBROUTINE Gauss_Nodes_Elem_2D(vert, gnodes)

IMPLICIT NONE

REAL(8), INTENT(IN)		::	vert(2,3)
REAL(8), INTENT(OUT)	::	gnodes(2,3)

REAL(8), PARAMETER		::	a = 0.5, b = 0.5

!a = (5.0+3.0*DSQRT(5.0))/20.0
!b = (5.0-DSQRT(5.0))/20.0

gnodes(:,1) = a*vert(:,1) + b*vert(:,2) 
gnodes(:,2) = a*vert(:,2) + b*vert(:,3)
gnodes(:,3) = a*vert(:,3) + b*vert(:,1)


END