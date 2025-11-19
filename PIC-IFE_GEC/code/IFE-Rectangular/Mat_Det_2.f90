SUBROUTINE Mat_Det_2(A,D)


IMPLICIT NONE

REAL(8), DIMENSION(2,2), INTENT(IN)	:: A
REAL(8), INTENT(OUT)				:: D

!D	=	A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))- &
!		A(2,1)*(A(1,2)*A(3,3)-A(1,3)*A(3,2))+ &	
!		A(3,1)*(A(1,2)*A(2,3)-A(1,3)*A(2,2))

D	=	A(1,1)*A(2,2)-A(1,2)*A(2,1)

!print*, A(1,:)
!print*, A(2,:)
!print*, D
!pause

END