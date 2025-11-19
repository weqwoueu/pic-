SUBROUTINE	Left_Divide(y, h, g)

!USE IFE_INTERFACE, ONLY: Matrix_Inverse

USE IFE_INTERFACE, ONLY: Matrix_Inverse
IMPLICIT NONE
!INCLUDE	'Matrix_Inverse.inc'

REAL(8), DIMENSION(:,:), INTENT(INOUT)		::	y
REAL(8), DIMENSION(:,:), INTENT(IN)			::	h
REAL(8), DIMENSION(:,:), INTENT(IN)			::	g
INTEGER										::	k, i, j
INTEGER										::	dim1, dim2
REAL(8), DIMENSION(:,:), POINTER			::	h_inverse



k = SIZE(g,1)
dim1 = SIZE(h,1)
dim2 = SIZE(h,2)

ALLOCATE(h_inverse(dim1,dim2))

DO i=1,dim1
	DO j=1,dim2
		h_inverse(i,j) = h(i,j)
	ENDDO
ENDDO

CALL Matrix_Inverse(h_inverse,k)

y = MATMUL(h_inverse,g)

DEALLOCATE(h_inverse)

END SUBROUTINE