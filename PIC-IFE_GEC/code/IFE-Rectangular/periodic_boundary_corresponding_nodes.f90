
SUBROUTINE periodic_boundary_corresponding_nodes(nnx,nny,bc_index,p_basic,node_type,ALQ,bc_point_1,bc_point_2,unknow_nodes)

IMPLICIT NONE

TYPE SPARSE

	REAL(8)	::	K
	INTEGER	::	JCOL, SROW
END TYPE

INTEGER nnx,nny
INTEGER, DIMENSION(:),		POINTER		    ::	bc_index
REAL(8), DIMENSION(:,:),	POINTER			::	p_basic
REAL(8), DIMENSION(:,:),	POINTER			::	ALQ
REAL(8), DIMENSION(:,:),    POINTER			::  bc_point_1,bc_point_2
INTEGER	,DIMENSION(:,:),	POINTER			::	unknow_nodes1
REAL(8)	,DIMENSION(:,:),	POINTER			::	unknow_nodes 
REAL(8) ,DIMENSION(:,:),	POINTER			::	check_periodic_boundary
INTEGER, DIMENSION(:,:),   POINTER			::	node_type
REAL xx

REAL(8) ,DIMENSION(:),		POINTER			::	I_COL,J_COL
!REAL(8),DIMENSION(:),		POINTER			::	record_nodes
INTEGER ,DIMENSION(:,:),	POINTER			::	ALQ_index

!INTEGER ,DIMENSION(:,:),	POINTER			::	record_pb
REAL(8)	SmallValue
INTEGER x_min,y_min,x_max,y_max,num_unknows,i,j,m,n,n_boundary
REAL x1_value_min,x1_value_max,y1_value_min,y1_value_max,x2_value_min,x2_value_max,y2_value_min,y2_value_max
INTEGER num_periodic_boundary,num,num_corresponding_boundary
num_unknows=SIZE(p_basic,2)
x_min=MINVAL(p_basic(1,:))
x_max=MAXVAL(p_basic(1,:))
y_min=MINVAL(p_basic(2,:))
y_max=MAXVAL(p_basic(2,:))

n_boundary=SIZE(bc_index)

num_periodic_boundary=0

ALLOCATE(ALQ_index(7,n_boundary))
ALQ_index=-200

DO i=1,n_boundary
	IF(bc_index(i)==-1)THEN
		num_periodic_boundary=num_periodic_boundary+1
		ALQ_index(1,i)=i
		!ALQ(1,num):alpha,ALQ(2,num):Q
		!xj=AL*xi+Q
		ALQ_index(2,i)=ALQ(1,i)
		ALQ_index(3,i)=ALQ(2,i)
		!ALQ_index(4:7,num):the num line's begin and end point
		ALQ_index(4,i)=bc_point_1(1,i)
		ALQ_index(5,i)=bc_point_1(2,i)
		ALQ_index(6,i)=bc_point_2(1,i)
		ALQ_index(7,i)=bc_point_2(2,i)

	ENDIF
ENDDO

!===================Check the boundary conditions===========================
IF(MOD(num_periodic_boundary,2)/=0)THEN
	WRITE(*,*) "Please check the periodic boudary conditions! "
	STOP
ENDIF
!===========================================================================
ALLOCATE(check_periodic_boundary(2,num_periodic_boundary/2))
check_periodic_boundary=-200


!record_pb(1:2,:) are the two periodic boundary conditions lines
!record_pb(3:4,:) are the AL and Q

!ALLOCATE(record_pb(4,num_periodic_boundary/2))
ALLOCATE(unknow_nodes1(4,nnx+nny))
unknow_nodes1=-200
SmallValue = 1.0D-5 
!record_pb=-200

num_corresponding_boundary=0
num=0
DO i=1,n_boundary
	IF(ABS(ALQ_index(4,i)-x_min)<SmallValue .AND. ABS(ALQ_index(6,i)-x_min)<SmallValue)THEN
		

		y1_value_min=MIN(ALQ_index(5,i),ALQ_index(7,i))
		y1_value_max=MAX(ALQ_index(5,i),ALQ_index(7,i))
		num_corresponding_boundary=num_corresponding_boundary+1
		check_periodic_boundary(1,num_corresponding_boundary)=i

		DO j=1,n_boundary

			y2_value_min=MIN(ALQ_index(5,j),ALQ_index(7,j))
			y2_value_max=MAX(ALQ_index(5,j),ALQ_index(7,j))
			!DO n=1,num_unknows
			!IF (p_basic(1,n)==x_min .AND. p_basic(2,n)>y1_value_min .AND. p_basic(2,n)<y1_value_max)THEN
			IF (ABS(ALQ_index(4,j)-x_max)<SmallValue .AND. ABS(ALQ_index(6,j)-x_max)<SmallValue .AND. & 
									ABS(y1_value_min-y2_value_min)<SmallValue .AND. ABS(y1_value_max-y2_value_max)<SmallValue)THEN
							
				check_periodic_boundary(2,num_corresponding_boundary)=j
				!num_corresponding_boundary=num_corresponding_boundary+1
				!record_pb(1,num_corresponding_boundary)=ALQ_index(1,i)
				!record_pb(2,num_corresponding_boundary)=ALQ_index(1,j)
				!record_pb(3,num_corresponding_boundary)=ALQ_index(2,i)
				!record_pb(4,num_corresponding_boundary)=ALQ_index(3,i)
				
				DO n=1,num_unknows
					IF (ABS(p_basic(1,n)-x_min)<SmallValue .AND. p_basic(2,n)>y1_value_min .AND. p_basic(2,n)<y1_value_max)THEN
						
						num=num+1
						unknow_nodes1(1,num)=n

						DO m=1,num_unknows
							IF (ABS(p_basic(1,m)-x_max)<SmallValue .AND. ABS(p_basic(2,n)-p_basic(2,m))<SmallValue)THEN
						
							!num2=num2+1
							!unknow_nodes1(1,num2)=n
							unknow_nodes1(2,num)=m
						
							unknow_nodes1(3,num)=ALQ_index(2,i)
							!unknows_nodes1(3,num2)=alpha
							unknow_nodes1(4,num)=ALQ_index(3,i)
							!unknows_nodes1(4,num2)=Q
							ENDIF
						ENDDO
					ENDIF
				ENDDO
			
			ENDIF
		ENDDO

	ELSEIF(ABS(ALQ_index(5,i)-y_min)<SmallValue .AND. ABS(ALQ_index(7,i)-y_min)<SmallValue)THEN
			
		x1_value_min=MIN(ALQ_index(4,i),ALQ_index(6,i))
		x1_value_max=MAX(ALQ_index(4,i),ALQ_index(6,i))
		num_corresponding_boundary=num_corresponding_boundary+1
		check_periodic_boundary(1,num_corresponding_boundary)=i
			
		DO j=1,SIZE(ALQ_index,2)

			x2_value_min=MIN(ALQ_index(4,j),ALQ_index(6,j))
			x2_value_max=MAX(ALQ_index(4,j),ALQ_index(6,j))
			
			IF (ABS(ALQ_index(5,j)-y_max)<SmallValue .AND. ABS(ALQ_index(7,j)-y_max)<SmallValue .AND. &
													ABS(x1_value_min-x2_value_min)<SmallValue .AND. ABS(x1_value_max-x2_value_max)<SmallValue)THEN
				!record_pb(1,num_corresponding_boundary)=ALQ_index(1,i)
				!record_pb(2,num_corresponding_boundary)=ALQ_index(1,j)
				!record_pb(3,num_corresponding_boundary)=ALQ_index(2,i)
				!record_pb(4,num_corresponding_boundary)=ALQ_index(3,i)
				check_periodic_boundary(2,num_corresponding_boundary)=j

				DO n=1,num_unknows
					IF (ABS(p_basic(2,n)-y_min)<SmallValue .AND. p_basic(1,n)>x1_value_min .AND. p_basic(1,n)<x1_value_max)THEN
						
						num=num+1
						unknow_nodes1(1,num)=n

						DO m=1,num_unknows
							IF (ABS(p_basic(2,m)-y_max)<SmallValue .AND. ABS(p_basic(1,n)-p_basic(1,m))<SmallValue)THEN
						
							!num2=num2+1
							!unknow_nodes1(1,num2)=n
							unknow_nodes1(2,num)=m
							unknow_nodes1(3,num)=ALQ_index(2,i)
							!unknows_nodes1(3,num2)=alpha
							unknow_nodes1(4,num)=ALQ_index(3,i)
							!unknows_nodes1(4,num2)=Q
							ENDIF
						ENDDO
					ENDIF
				ENDDO
			ENDIF
		ENDDO	
	ENDIF
ENDDO
!===================Check the boundary conditions===========================
DO i=1,num
	IF(unknow_nodes1(1,i)<0 .OR. unknow_nodes1(2,i)<0)THEN
		WRITE(*,*) "Please check the periodic boundary conditions or the start and end points of the boundary!"
		STOP
	ENDIF
ENDDO
DO i=1,num_corresponding_boundary
	IF(check_periodic_boundary(2,num_corresponding_boundary)<0 .OR. check_periodic_boundary(1,num_corresponding_boundary)<0)THEN
		WRITE(*,*) "Please check the periodic boundary conditions or the start and end points of the boundary!"
		STOP
	ENDIF
ENDDO

!===========================================================================

IF(num>0)THEN
	ALLOCATE(unknow_nodes(4,num))
	DO i=1,num
		unknow_nodes(:,i)=unknow_nodes1(:,i)
	ENDDO
	
	
ENDIF


DEALLOCATE(ALQ_index,check_periodic_boundary, unknow_nodes1)

END
				
				

