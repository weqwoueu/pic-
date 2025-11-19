SUBROUTINE Sparse_Residue(A, X, B, NZ, N, Residue)

USE IFE_MAIN_PARAM

IMPLICIT NONE

INTEGER								NZ, N
TYPE(SPARSE), DIMENSION(NZ)		::	A
REAL(8), DIMENSION(N)			::	X
REAL(8), DIMENSION(N)			::	B
REAL(8)								Residue

INTEGER	i, j
REAL(8) row_sum
REAL*8, DIMENSION(:), POINTER		::	R

ALLOCATE(R(N))

DO i = 1, N
	row_sum = Zero
	DO j = A(i)%SROW, A(i+1)%SROW-1
		row_sum = row_sum + A(j)%K * X(A(j)%JCOL)
	END DO
	R(i) = B(i) - row_sum
END DO

!Residue = Zero
!DO i = 1, N
	!Residue = Residue + R(i)**2
!END DO
Residue = MAXVAL(DABS(R))

DEALLOCATE(R) !$ ab.ZWZ 2021/7/20

END