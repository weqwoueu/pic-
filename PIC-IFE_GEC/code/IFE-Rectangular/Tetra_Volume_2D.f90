SUBROUTINE Tetra_Volume_2D(vert, V)

USE IFE_MAIN_PARAM

IMPLICIT NONE

REAL(8), DIMENSION(2,3), INTENT(IN)	::	vert
REAL(8), INTENT(OUT)				::	V

REAL(8), DIMENSION(:), POINTER	::	a, b !, c

ALLOCATE(a(2),b(2))

a = vert(:,2) - vert(:,1)
b = vert(:,3) - vert(:,1)

!print*, a
!print*, b

CALL Mat_Det_2(RESHAPE((/a,b/),(/2,2/)),V)

V = DABS(V)/2.0

DEALLOCATE(a,b)

END