SUBROUTINE Diagonal_Preconditioner(K_stiff, K_diag)

USE IFE_MAIN_PARAM

IMPLICIT NONE

TYPE(SPARSE), DIMENSION(:), INTENT(IN)	::	K_stiff
REAL(8), DIMENSION(:), POINTER			::	K_diag

INTEGER	i, j, n_unknowns

n_unknowns = MAXVAL(K_stiff(:)%JCOL)

ALLOCATE(K_diag(n_unknowns))

DO i=1,n_unknowns
    IF (n_unknowns ==1) THEN
	j=1
		IF(K_stiff(j)%JCOL==i)THEN
			K_diag(i) = K_stiff(j)%K
			EXIT
		END IF
	ELSE
	DO j=K_stiff(i)%SROW,K_stiff(i+1)%SROW-1
		IF(K_stiff(j)%JCOL==i)THEN
			K_diag(i) = K_stiff(j)%K
			EXIT
		END IF
	END DO
	ENDIF
END DO

END