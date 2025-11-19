SUBROUTINE MY_PCG(ANEL, A, AI, AJ, KNEL, K, KI, KJ, N, R, X, RELTOL, MAXITER, ITER, IERR, MATVEC, MSOLVE)

! MY_PCG:	Preconditioned Conjugate Gradient sparse linear solver.
!	
! PURPOSE:	Solve the system A*x = b, where A is a symmetric, positive
!			definite matrix.	
!	A(ANEL)		REAL(8)		INPUT			System Matrix
!	AI(ANEL)	INTEGER		INPUT			Sparse Structure Columns
!	AJ(ANEL)	INTEGER		INPUT			Sparse Structure Rows
!	K(KNEL)		REAL(8)		INPUT			Preconditioner
!	KI(KNEL)	INTEGER		INPUT			Sparse Structure Columns
!	KJ(KNEL)	INTEGER		INPUT			Sparse Structure Rows
!	R(N)		REAL(8)		INPUT/OUTPUT	RHS/Residue
!	X(N)		REAL(8)		INPUT/OUTPUT	Guess/Solution
!	RELTOL		REAL(8)		INPUT			Relative Tolerance
!	MAXITER		INTEGER		INPUT			Maximum Iterations
!	MATVEC		EXTERNAL	INPUT			Matix-Vector Multiplication Subroutine
!	MSOLVE		EXTERNAL	INPUT			Linear Solver Subroutine (K * w = r)
!	ITER		INTEGER		OUTPUT			Number of Iterations
!	IERR		INTEGER		OUTPUT			Error Index

IMPLICIT NONE

INTEGER				ANEL, KNEL, N
REAL(8), DIMENSION(:)	::	A, K, R, X
INTEGER, DIMENSION(:)	::	AI, AJ, KI, KJ
REAL(8)				RELTOL
INTEGER				MAXITER, ITER, IERR
! VERY IMPORTANT !!
! Procedures which have explicit INTERFACE should not be declared EXTERNAL
!EXTERNAL			MATVEC, MSOLVE

INTERFACE
	SUBROUTINE MATVEC(NEL, A, AI, AJ, N, X, Y)
	INTEGER		NEL, N
	REAL(8), DIMENSION(:)	::	A
	INTEGER, DIMENSION(:)	::	AI, AJ
	REAL(8), DIMENSION(:)	::	X, Y
	END SUBROUTINE

	SUBROUTINE MSOLVE(NEL, A, AI, AJ, N, Y, X)
	INTEGER		NEL, N
	REAL(8), DIMENSION(:)	::	A
	INTEGER, DIMENSION(:)	::	AI, AJ
	REAL(8), DIMENSION(:)	::	X, Y
	END SUBROUTINE
END INTERFACE


REAL(8), DIMENSION(:), POINTER	::	X0, R0, W, P, Q
REAL(8)		alpha, beta, rho, rho0, normR, normRHS, RELTOL2
INTEGER		i, j

RELTOL2 = RELTOL**2

ALLOCATE(X0(N), R0(N), W(N), P(N), Q(N))

normRHS = 0
DO j = 1, N
	normRHS = normRHS + R(j)*R(j)
END DO

DO j = 1, N
	X0(j) = X(j)							! X0 = X
END DO

CALL MATVEC(ANEL, A, AI, AJ, N, X0, R0)
DO j = 1, N
	R0(j) = R(j) - R0(j)					! R(0) = B - A * X(0)
END DO


DO i = 1, MAXITER
	
	! K is not Singular
	CALL MSOLVE(KNEL, K, KI, KJ, N, R0, W)
	
	rho = 0
	DO j = 1, N
		rho = rho + R0(j)*W(j)
	END DO									! rho(i-1) = (R(i))T * W(i-1)

	IF (i == 1) THEN
		DO j = 1, N
			P(j) = W(j)						! P(i) = W(i-1)
		END DO
	ELSE
		IF (rho0 == 0) THEN 
			IERR = 2
			EXIT
		END IF
		! rho0 =/= 0
		beta = rho/rho0
		DO j = 1, N
			P(j) = W(j)+beta*P(j)			! P(i) = W(i-1) + beta(i-1)*P(i-1)
		END DO
	ENDIF
	
	CALL MATVEC(ANEL, A, AI, AJ, N, P, Q)			! Q(i) = A * P(i)
	alpha = 0
	DO j = 1, N
		alpha = alpha + P(j)*Q(j)
	END DO									! alpha(i) = rho(i-1)/(P(i))T * Q(i)

	! alpha =/= 0
	IF (alpha == 0) THEN 
		IERR = 3
		EXIT
	END IF
	alpha = rho/alpha						! alpha(i) = rho(i-1)/(P(i))T * Q(i)

	DO j = 1, N
		X0(j) = X0(j) + alpha*P(j)
		R0(j) = R0(j) - alpha*Q(j)
	END DO									! X(i) = X(i-1) + alpha(i)*P(i)
											! R(i) = R(i-1) - alpha(i)*Q(i)
	rho0 = rho

	normR = 0
	DO j = 1, N
		normR = normR + R0(j)*R0(j)		
	END DO
	IF (normR/normRHS <= RELTOL2) THEN
		EXIT
	END IF

END DO

IF (i >= MAXITER) THEN
	IERR = 1
ELSE
	ITER = i-1
END IF

DO j = 1, N
	X(j) = X0(j)
	R(j) = R0(j)
END DO			

DEALLOCATE(X0, R0, W, P, Q)

END SUBROUTINE
