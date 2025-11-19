SUBROUTINE Gauss_Nodes_Edge_2D(vert1, vert2, gnodes)

IMPLICIT NONE

REAL(8), INTENT(IN)		::	vert1(2), vert2(2)
REAL(8), INTENT(OUT)	::	gnodes(2,3)

REAL(8)		::	a (3)


a(1) = 0.5*(1.0 - 0.2*DSQRT(15.0D0))
a(2) = 0.5
a(3) = 0.5*(1.0 + 0.2*DSQRT(15.0D0))

gnodes(:,1) = vert1(:) + a(1)*(vert2(:) - vert1(:)) 
gnodes(:,2) = vert1(:) + a(2)*(vert2(:) - vert1(:)) 
gnodes(:,3) = vert1(:) + a(3)*(vert2(:) - vert1(:)) 


END