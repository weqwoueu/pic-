SUBROUTINE Linear_IFE_Basis_Coeff_2D( vert, intrs, confg_ind, ElBeta, coef)

!USE MSIMSL
USE IFE_MAIN_PARAM

IMPLICIT NONE

REAL(8), INTENT(IN)					::	vert(2,4), intrs(2,2), ElBeta(2) !,tri_sign(2)										
INTEGER									confg_ind
REAL(8), INTENT(OUT)				::	coef(4,8)

REAL(8), DIMENSION(:), ALLOCATABLE		::	rhs, n
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	M, Mtemp
INTEGER										basis_ind, L

ALLOCATE(rhs(8))

ALLOCATE(M(8,8), Mtemp(8,8))
M = Zero

ALLOCATE(n(2))

IF (confg_ind==2) THEN  !cut the diagnal edge
   M(1,2:3) = (vert(:,1));	 M(1,4) =vert(1,1) * vert(2,1);
   M(4,2:3) = (vert(:,4));  M(4,4) =vert(1,4) * vert(2,4);
   M(2,6:7) = (vert(:,2));  M(2,8) =vert(1,2) * vert(2,2); 
   M(3,6:7) = (vert(:,3));  M(3,8) =vert(1,3) * vert(2,3); 
!   M(3,8) =vert(1,3) * vert(2,3);
   M(1,1)   = One;    M(4,1)   = One;	   M(2:3,5)   = One;       !	 M(1:2,4:6) = Zero   

!   M(3,1:3)   = Zero;					 M(3,4:5)   = vert(:,3); M(3,6)	    = One
ELSEIF (confg_ind==1) THEN !cut the NEXT edge

   M(1,2:3) = (vert(:,1));	 M(1,4) =vert(1,1) * vert(2,1);
 !  M(4,2:3) = TRANSPOSE(vert(:,4));  M(4,4) =vert(1,4) * vert(2,4)
!   M(2:4,6:7) = TRANSPOSE(vert(:,2:4));  M(2,8) =vert(1,2) * vert(2,2); M(3,8) =vert(1,3) * vert(2,3); M(4,8) =vert(1,4) * vert(2,4)
   M(4,6:7) = (vert(:,4));  M(4,8) =vert(1,4) * vert(2,4);
   M(2,6:7) = (vert(:,2));  M(2,8) =vert(1,2) * vert(2,2); 
   M(3,6:7) = (vert(:,3));  M(3,8) =vert(1,3) * vert(2,3); 
   M(1,1)   = One;    M(2:4,5)   = One;       !	 M(1:2,4:6) = Zero
 
END IF

M(5:7,1)	  = ONE  
M(5:6,2:3)	  = intrs(:,1:2); M(5,4) = intrs(1,1)*intrs(1,2); M(6,4) = intrs(2,1)*intrs(2,2) 
M(7,2) = (intrs(1,1) +intrs(2,1))/2;  M(7,3) = (intrs(1,2) +intrs(2,2))/2;   
M(7,4)	 =  ((intrs(1,1) +intrs(2,1))/2) *((intrs(1,2) +intrs(2,2))/2);
M(5:7,5)	  = -ONE  
M(5:6,6:7)	  = -intrs(:,1:2); M(5,8) = -intrs(1,1)*intrs(1,2); M(6,8) = -intrs(2,1)*intrs(2,2) 
M(7,6) = -(intrs(1,1) +intrs(2,1))/2;  M(7,7) = -(intrs(1,2) +intrs(2,2))/2;   
M(7,8)	 =  -((intrs(1,1) +intrs(2,1))/2) *((intrs(1,2) +intrs(2,2))/2);
M(8,1) = 0;  M(8,2) = ElBeta(1)*(intrs(2, 2)-intrs(1,2)); M(8,3) = -ElBeta(1)*(intrs(2,1)-intrs(1,1));
M(8,4) = ElBeta(1)*(((intrs(2,2)**2)-(intrs(1,2)**2)-(intrs(2,1)**2)+(intrs(1,1))**2))*.5;
M(8,5) = 0;  M(8,6) = -ElBeta(2)*(intrs(2,2)-intrs(1,2)); M(8,7) = ElBeta(2)*(intrs(2,1)-intrs(1,1));
M(8,8) = -ElBeta(2)*(((intrs(2,2))**2-(intrs(1,2)**2)-(intrs(2,1)**2)+(intrs(1,1))**2))*.5;


IF (intrs(1,1)==intrs(2,1) .OR. intrs(1,2)==intrs(2,2)) THEN
   M(7,1)=0;M(7,2)=0;M(7,3)=0;M(7,4)=1;M(7,5)=0;M(7,6)=0;M(7,7)=0;M(7,8)=-1;
ENDIF

DO basis_ind=1,4
	rhs = Zero
	rhs(basis_ind) = One
	Mtemp = M

!	! CALL DLSARG (N, A, LDA, B, IPATH, X)
!	! IPATH = 1 means the system A*X = B is solved.
!	! IPATH = 2 means the system AT*X = B is solved.
!
!	! PC Version
!	CALL DLSARG (8, M, 8, rhs, 1, coef(basis_ind,:))
!	! CRAY Version
!	!CALL LSARG (4, M, 4, rhs, 1, coef)

    CALL AGAUS (Mtemp, rhs, 8, coef(basis_ind,:), L)
	IF (L==0) THEN
	  PRINT*, ' Can not Slove the coef, Something is Wrong'
	  STOP
	END IF
END DO

DEALLOCATE(rhs)
DEALLOCATE(M, Mtemp)
DEALLOCATE(n)

END
