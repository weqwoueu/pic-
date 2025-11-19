SUBROUTINE Linear_FE_Basis_Coeff_2D(vert, coef)

USE IFE_MAIN_PARAM

IMPLICIT NONE

REAL(8), INTENT(IN)		::	vert(2,4)
REAL(8), INTENT(OUT)	::	coef(4,4)

REAL(8), DIMENSION(:), ALLOCATABLE		::	rhs
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	M, Mtemp
INTEGER										basis_ind, L, i					

ALLOCATE(rhs(4))


ALLOCATE(M(4,4), Mtemp(4,4))
!M(:,1:2)	= TRANSPOSE(vert)
!M(:,3)		= One
M(:,1)		= One
M(:,2:3)	= TRANSPOSE(vert)
!M(:,4)      = TRANSPOSE(vert(1,:)*vert(2,:))
DO i= 1,4
   M(i,4)      = vert(1,i)* vert(2,i)
ENDDO
!N=a+bx+cy+dxy

!print*, vert(1,:)
!print*, vert(2,:)

DO 	basis_ind=1,4
	rhs = Zero
	rhs(basis_ind) = One
	Mtemp = M
	! CALL DLSARG (N, A, LDA, B, IPATH, X)
	! IPATH = 1 means the system A*X = B is solved.
	! IPATH = 2 means the system AT*X = B is solved.

	! PC Version
    ! CALL DLSARG (4, M, 4, rhs, 1, coef(basis_ind,:))
!print*, M(1,:), RHS(1)
!print*, M(2,:), RHS(2)
!print*, M(3,:), RHS(3)
    CALL AGAUS (Mtemp, rhs, 4, coef(basis_ind,:), L)
	IF (L==0) THEN
	  PRINT*, ' Can not Slove the coef, Something is Wrong'
	  STOP
	END IF
!print*, basis_ind, coef(basis_ind,:)
!pause
	! CRAY Version
	! CALL LSARG (4, M, 4, rhs, 1, coef)



END DO

DEALLOCATE(rhs)
DEALLOCATE(M, Mtemp)

END
