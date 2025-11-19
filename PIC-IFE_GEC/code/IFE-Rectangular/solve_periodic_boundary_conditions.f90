SUBROUTINE solve_periodic_boundary_conditions(i,j,AL,Q,num_unknows ,A_stiff,f)
          
!Purpose:
!	removal of constraint periodic boundary condtions
!	Akj=Ajk=0,Ajj=1,Aji=Aij=-alpha
!	Aki=Aik=Aki+alpha*Akj
!	Aii=alpha*c+Aii+alpha*Aij, others Akl=Alk
!	Fj=-Q,Fi=c*Q+Fi+alpha*Fj
!
!	c=(1+Ajj)*alpha+Aij
!	xj=alpha*xi+Q
TYPE SPARSE
	REAL(8)	::	K
	INTEGER	::	JCOL, SROW
END TYPE

REAL(8),DIMENSION(:),		POINTER		::	f
REAL(8),DIMENSION(:),		POINTER		::	con_f
REAL(8),DIMENSION(:),		POINTER		::	I_COL,J_COL
REAL(8),DIMENSION(:,:),	POINTER		::	RECORD_I,RECORD_II,RECORD_J
REAL(8),DIMENSION(2,2)                   ::  RECORD_J_2_2

TYPE(SPARSE), DIMENSION(:), POINTER			 ::	A_stiff
TYPE(SPARSE), DIMENSION(:), POINTER			 ::con_matrix
INTEGER num_unknows,m,n,num1,add_num,num2,i,j,num_add_jrow
REAL(8) x,AJJ,AII,AIJ,AL,FJ
REAL(8)	SmallValue

RECORD_J_2_2 = 0.0


!num_unknows=SIZE(A_stiff)
!the i and j col of the matrix 
ALLOCATE(I_COL(num_unknows))
ALLOCATE(J_COL(num_unknows))
SmallValue = 1.0D-5 
I_COL=0
J_COL=0
DO n=1,num_unknows
	DO m=A_stiff(n)%SROW,A_stiff(n+1)%SROW-1
		IF(A_stiff(m)%JCOL==j)THEN
			J_COL(n)=A_stiff(m)%K
		ELSEIF(A_stiff(m)%JCOL==i)THEN
			I_COL(n)=A_stiff(m)%K
		ENDIF
	ENDDO
ENDDO
!1 to j-1 row
DO m=1,j-1
	IF(J_COL(m)/=0 .AND. ABS(I_COL(m))<SmallValue)THEN !I_COL(m)==0
		DO n=A_stiff(m)%SROW,A_stiff(m+1)%SROW-1
			IF(ABS(A_stiff(n)%JCOL-j)<SmallValue)THEN !A_stiff(n)%JCOL==j
				A_stiff(n)%JCOL=i
				A_stiff(n)%K=AL*J_COL(m)
			ENDIF
		ENDDO
	ELSEIF(J_COL(m)/=0 .AND. I_COL(m)/=0)THEN
		DO n=A_stiff(m)%SROW,A_stiff(m+1)%SROW-1
				IF(ABS(A_stiff(n)%JCOL-j)<SmallValue)THEN   !A_stiff(n)%JCOL==j
					x=A_stiff(n)%K
					A_stiff(n)%K=0
				ELSEIF(ABS(A_stiff(n)%JCOL-i)<SmallValue)THEN  !A_stiff(n)%JCOL==i
					A_stiff(n)%K=A_stiff(n)%K+AL*x

			    ENDIF
		ENDDO
	ENDIF
ENDDO

!j row
num_add_jrow=0
DO n=A_stiff(j)%SROW,A_stiff(j+1)%SROW-1

	IF(J_COL(j)/=0 .AND. I_COL(j)/=0)THEN
		num_add_jrow=0
		IF(ABS(A_stiff(n)%JCOL-j)<SmallValue)THEN  !A_stiff(n)%JCOL==j
			AJJ=A_stiff(n)%K
			A_stiff(n)%K=1
		ELSEIF(ABS(A_stiff(n)%JCOL-i)<SmallValue)THEN  !A_stiff(n)%JCOL==i
			A_stiff(n)%K=-AL
		ELSEIF(A_stiff(n)%JCOL/=i .AND. A_stiff(n)%JCOL/=j)THEN
			A_stiff(n)%K=0
		ENDIF
	ELSEIF(ABS(J_COL(j))<SmallVAlue .AND. I_COL(j)/=0)THEN  !J_COL(j)==0
		num_add_jrow=1
!		ALLOCATE(RECORD_J(2,1))
		IF(ABS(A_stiff(n)%JCOL-i)<SmallValue)THEN  !A_stiff(n)%JCOL==i
			A_stiff(n)%K=-AL
		ELSEIF(A_stiff(n)%JCOL/=i .AND. A_stiff(n)%JCOL/=j)THEN
			A_stiff(n)%K=0
		ENDIF
		RECORD_J_2_2(1,1)=j
		RECORD_J_2_2(2,1)=1
		AJJ=0
	ELSEIF(J_COL(j)/=0 .AND. ABS(I_COL(j))<SmallValue)THEN  !I_COL(j)==0
		num_add_jrow=1
!		ALLOCATE(RECORD_J(2,1))
		IF(ABS(A_stiff(n)%JCOL-j)<SmallValue)THEN  !A_stiff(n)%JCOL==j
			AJJ=A_stiff(n)%K
			A_stiff(n)%K=1
		ELSEIF(A_stiff(n)%JCOL/=i .AND. A_stiff(n)%JCOL/=j)THEN
			A_stiff(n)%K=0
		ENDIF
		RECORD_J_2_2(1,1)=i
		RECORD_J_2_2(2,1)=-1
	ELSEIF(ABS(J_COL(j))<SmallValue .AND. ABS(I_COL(j))<SmallValue)THEN  !J_COL(j)==0 .AND. I_COL(j)==0
		num_add_jrow=2
!		ALLOCATE(RECORD_J(2,2))
		AJJ=0
		RECORD_J_2_2(1,1)=j
		RECORD_J_2_2(2,1)=1
		RECORD_J_2_2(1,2)=i
		RECORD_J_2_2(2,2)=-1
		A_stiff(n)%K=0
	ENDIF

ENDDO
!from j+1 to i-1 row
DO m=j+1,i-1
	
	IF(J_COL(m)/=0 .AND. ABS(I_COL(m))<SmallValue)THEN !I_COL(m)==0
		DO n=A_stiff(m)%SROW,A_stiff(m+1)%SROW-1
			IF(ABS(A_stiff(n)%JCOL-j)<SmallValue)THEN  !A_stiff(n)%JCOL==j
				A_stiff(n)%JCOL=i
				A_stiff(n)%K=AL*J_COL(m)
			ENDIF
		ENDDO
	ELSEIF(J_COL(m)/=0 .AND. I_COL(m)/=0)THEN
		DO n=A_stiff(m)%SROW,A_stiff(m+1)%SROW-1
				IF(ABS(A_stiff(n)%JCOL-j)<SmallValue)THEN  !A_stiff(n)%JCOL==j
					x=A_stiff(n)%K
					A_stiff(n)%K=0
				ELSEIF(ABS(A_stiff(n)%JCOL-i)<SmallValue)THEN !A_stiff(n)%JCOL==i
					A_stiff(n)%K=A_stiff(n)%K+AL*x

			    ENDIF
		ENDDO
	ENDIF
ENDDO
ALLOCATE(RECORD_I(2,num_unknows))
RECORD_I=0
! the i row
m=i
num1=0
DO n=1,num_unknows
	IF(n==j)THEN
		AIJ=I_COL(n)
		I_COL(n)=-1
	ELSEIF(n==i) THEN
		AII=I_COL(n)
		I_COL(n)=AII+AL*AL*(1+AJJ)+AL*AIJ+AL*AIJ
	ELSEIF(n/=i .AND. n/=j)THEN
		I_COL(n)=I_COL(n)+AL*J_COL(n)
	ENDIF
ENDDO
DO n=1,num_unknows
	IF(I_COL(n)/=0)THEN
		num1=num1+1
		RECORD_I(1,n)=I_COL(n)
		RECORD_I(2,n)=n
	ENDIF
ENDDO

add_num=num1-(A_stiff(i+1)%SROW-A_stiff(i)%SROW)

ALLOCATE(RECORD_II(2,num1))
num2=0
IF(add_num>0)THEN
	DO m=1,num_unknows
		IF(RECORD_I(1,m)/=0 )THEN
			num2=num2+1
			RECORD_II(:,num2)=RECORD_I(:,m)
		ENDIF
	ENDDO
ENDIF	
num2=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(add_num<0)THEN
	add_num=0
ENDIF

IF(add_num+num_add_jrow>0)THEN

	IF(ABS(num_add_jrow)<SmallValue)THEN  !num_add_jrow==0

		ALLOCATE(con_matrix(size(A_stiff)+add_num))

		DO n=size(A_stiff),A_stiff(i+1)%SROW,-1

			con_matrix(n+add_num)%JCOL=A_stiff(n)%JCOL
			con_matrix(n+add_num)%K=A_stiff(n)%K

		ENDDO
		DO n=A_stiff(i)%SROW,A_stiff(i+1)%SROW-1+add_num
			num2=num2+1

			con_matrix(n)%JCOL=RECORD_II(2,num2)
			con_matrix(n)%K=RECORD_II(1,num2)
		ENDDO
		DO m=1,A_stiff(i)%SROW-1
			con_matrix(m)%JCOL=A_stiff(m)%JCOL
			con_matrix(m)%K=A_stiff(m)%K
		ENDDO
		DO m=1,i
			con_matrix(m)%SROW=A_stiff(m)%SROW
		ENDDO
		DO m=i+1,num_unknows+1
			con_matrix(m)%SROW=A_stiff(m)%SROW+add_num
		ENDDO

	ELSEIF(num_add_jrow>0)THEN
		ALLOCATE(con_matrix(size(A_stiff)+add_num+num_add_jrow))
		! from i row to the end row
		DO n=size(A_stiff),A_stiff(i+1)%SROW,-1
	
			con_matrix(n+add_num+num_add_jrow)%JCOL=A_stiff(n)%JCOL
			con_matrix(n+add_num+num_add_jrow)%K=A_stiff(n)%K
		ENDDO
		! i row
		DO n=A_stiff(i)%SROW+num_add_jrow,A_stiff(i+1)%SROW-1+add_num+num_add_jrow
			num2=num2+1
			con_matrix(n)%JCOL=RECORD_II(2,num2)
			con_matrix(n)%K=RECORD_II(1,num2)
		ENDDO
		!FROM 1 TO J ROW
		DO m=1,A_stiff(j+1)%SROW-1
			con_matrix(m)%JCOL=A_stiff(m)%JCOL
			con_matrix(m)%K=A_stiff(m)%K
		ENDDO
		
		DO m=1,num_add_jrow
			con_matrix(A_stiff(j+1)%SROW-1+m)%JCOL=RECORD_J_2_2(1,m)
			con_matrix(A_stiff(j+1)%SROW-1+m)%K=RECORD_J_2_2(2,m)
		ENDDO
		!j+1 to i-1
		DO m=A_stiff(j+1)%SROW,A_stiff(i)%SROW-1
			con_matrix(m+num_add_jrow)%JCOL=A_stiff(m)%JCOL
			con_matrix(m+num_add_jrow)%K=A_stiff(m)%K
		ENDDO
       !====================GET THE VALUE OF con_matrix(n)%SROW==============
        DO m=1,j
			con_matrix(m)%SROW=A_stiff(m)%SROW
		ENDDO
		DO m=j+1,i
			con_matrix(m)%SROW=A_stiff(m)%SROW+num_add_jrow
		ENDDO
					
		DO m=i+1,num_unknows+1
			con_matrix(m)%SROW=A_stiff(m)%SROW+add_num+num_add_jrow
		ENDDO
	ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSEIF(add_num+num_add_jrow<=0)THEN
	ALLOCATE(con_matrix(size(A_stiff)))
	DO m=1,size(A_stiff)
		con_matrix(m)%SROW=A_stiff(m)%SROW
		con_matrix(m)%JCOL=A_stiff(m)%JCOL
		con_matrix(m)%K=A_stiff(m)%K
	ENDDO
ENDIF


DEALLOCATE(A_stiff)
ALLOCATE(A_stiff(size(con_matrix)))
DO m=1,size(con_matrix)
	A_stiff(m)%SROW=con_matrix(m)%SROW
	A_stiff(m)%JCOL=con_matrix(m)%JCOL
	A_stiff(m)%K=con_matrix(m)%K

ENDDO


ALLOCATE(con_f(size(f)))
DO m=1,size(f)
	IF(ABS(m-j)<SmallValue)THEN  !m==j
		FJ=f(m)
		con_f(m)=-Q
		
	ELSEIF(ABS(m-i)<SmallValue)THEN !m==i
		con_f(m)=f(m)+AL*FJ+((1+AJJ)*AL+AIJ)*Q
	ELSEIf(m/=j .AND. m/=i)THEN
		con_f(m)=f(m)
	ENDIF
ENDDO
DO m=1,size(f)
	f(m)=con_f(m)
ENDDO
DEALLOCATE(RECORD_II)
DEALLOCATE(con_f)
DEALLOCATE(I_COL,J_COL)
DEALLOCATE(con_matrix)
DEALLOCATE(RECORD_I)
END
		

		
	
			
			
		


	

