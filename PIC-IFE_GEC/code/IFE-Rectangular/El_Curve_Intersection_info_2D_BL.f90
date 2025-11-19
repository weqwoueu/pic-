SUBROUTINE El_Curve_Intersection_info_2D_BL(vert, object, P_intrs, el_type, el_region, information_1_el, &
                                            information_2_el, node_index_el)

! Purpose:		Intersect a generally defined object with a tetrahedal element BY CHJ


USE IFE_MAIN_PARAM
USE Object_Data_2D
IMPLICIT NONE



REAL(8), DIMENSION(2,4), INTENT(IN)		::	vert

TYPE(ObjectType), INTENT(IN)			    ::	object
REAL(8), DIMENSION(4,6)					::	P_intrs
INTEGER										el_region
INTEGER										n_section,n_nodes

REAL(8), DIMENSION(2,2)					::	xyp

REAL(8)	        ElBound(2,2),LineEnds(2,2)
REAL(8)	        inters_point(2),inters_point_tmp(2)
INTEGER, DIMENSION(:), POINTER	        ::  information_1_el
REAL(8), DIMENSION(:), POINTER	        ::  information_2_el
INTEGER, DIMENSION(:), INTENT(IN)		::	node_index_el
REAL(8)									::	nx, ny, magnitude

REAL(8) x1,y1, x2,y2, diff1, diff2, length,ytemp
REAL(8) XBound(2)

REAL(8)     N1(2), N2(2), N3(2), pm(2), N4(2)

INTEGER e, inters_flag, nint1, nint2, i, j,count, inters_flag_tmp, count2

REAL(8)     edges(4,2)

REAL(8) Norm_2,xi

LOGICAL AtEnd1, AtEnd2,OutOfBound


!CYC add
INTEGER el_type
INTEGER pointer_local_to_reference(4), pointer_reference_to_local(4)
INTEGER	 plus_piece_flag, minus_piece_flag, first_piece_flag, second_piece_flag
INTEGER temp_pointer1(4),temp_pointer2(4),temp_pointer3(4),temp_pointer4(4)
REAL(8)	 beta_minus, beta_plus, beta1, beta2
REAL(8) Dx, Ex, Dy, Ey
INTEGER vert_index(4), edge_index(4)
REAL(8) xy(2)
!CYC add



P_intrs		= Zero

el_region	= Zero

beta_plus = object%Eps
beta_minus = ONE

DO i=1,4
	vert_index(i) = 0
	edge_index(i) = 0
ENDDO

edges = RESHAPE((/	1,	2,	&
					2,	3,	&
					3,	4,	&
					4,	1	/), (/4,2/), ORDER=(/2,1/))

DATA temp_pointer1 /1, 2, 3, 4/
DATA temp_pointer2 /2, 3, 4, 1/
DATA temp_pointer3 /3, 4, 1, 2/
DATA temp_pointer4 /4, 1, 2, 3/


n_nodes	=	object%Wall(1)%nxp



ElBound(1,1)  = MINVAL(vert(1,:))
ElBound(2,1)  = MINVAL(vert(2,:))
ElBound(1,2)  = MAXVAL(vert(1,:))
ElBound(2,2)  = MAXVAL(vert(2,:))

XBound(1)=object%Locations(1,1)
xy=0
IF (object%Wall(1)%Shape==2) THEN
	n_section	=	n_nodes	
	XBound(2)=MAXVAL(object%Wall(1)%node(1:n_nodes,1))
ELSE
	n_section	=	n_nodes-1
	XBound(2)=object%Locations(2,1)
!!	IF(MAXVAL(object%Wall(1)%node(1:n_nodes,1))>object%Locations(2,1))THEN
!!	
!!	    DO i=1,n_nodes
!!	      IF (object%Wall(1)%node(i,1)>=object%Locations(2,1)) THEN     !xmax
!!		    n_section=i
!!		    n_section=n_section-1
!!		    IF (object%Wall(1)%node(i,1)==object%Locations(2,1)) THEN
!!		        ytemp=object%Wall(1)%node(i,2)
!!		    ELSEIF (object%Wall(1)%node(i,1)>object%Locations(2,1)) THEN
!!		        xy(1:2)=object%Wall(1)%node(i,1:2)
!!		        ytemp=(xy(2)-object%Wall(1)%node(i-1,2))/(xy(1)-object%Wall(1)%node(i-1,1))
!!		        ytemp=ytemp*(object%Locations(2,1)-xy(1))+xy(2)
!!		    ENDIF
!!!		    print *,n_section,xy,ytemp,object%Wall(1)%node(i,1:2),object%Locations(2,1)
!!		    EXIT
!!	      ENDIF
!!	    ENDDO
!!		n_section=n_nodes       !n_section+1	
!!!	    xy(1:2)=object%Wall(1)%node(n_nodes-2,1:2)
!!!	    object%Wall(1)%node(n_nodes-2,1)=object%Locations(2,1)
!!!		object%Wall(1)%node(n_nodes-2,2)=(xy(1)-object%Locations(2,1)) &
!!!				        /(xy(1)-object%Wall(1)%node(n_nodes-3,1))
!!!		object%Wall(1)%node(n_nodes-2,2)=object%Wall(1)%node(n_nodes-2,2)*    &
!!!		                    (xy(2)-object%Wall(1)%node(n_nodes-3,2))
!!!		object%Wall(1)%node(n_nodes-2,2)=object%Wall(1)%node(n_nodes-2,2)+xy(2)
!!!	    object%Wall(1)%node(n_nodes-1,1)=object%Locations(2,1)
!!!	    print *,xy,object%Wall(1)%node(n_nodes-2,1:2),XBound(2)
!!	ENDIF
ENDIF

!print *,'pass113'

DO e=1,4

	x1 = vert(1,edges(e,1))
	y1 = vert(2,edges(e,1))
   
	x2 = vert(1,edges(e,2))
	y2 = vert(2,edges(e,2))
 
	LineEnds(1,:) = (/x1, y1/)
	LineEnds(2,:) = (/x2, y2/)

	count= 0
	inters_flag = 0
	inters_point = 0

	DO i = 1, n_section     !-1	
		inters_flag_tmp = 0
		OutOfBound = .FALSE.
!		IF (i==n_section-1) THEN
!			xyp(1,1)=object%Locations(2,1)
!	        xyp(1,2)=object%Wall(1)%node(n_nodes-1,2)
!			xyp(2,1:2)=object%Wall(1)%node(n_nodes-1,1:2)
!	    ELSEIF (i==n_section-2) THEN
!		    xyp(1:2,1)=object%Locations(2,1)
!			xyp(1,2)=ytemp
!			xyp(2,2)=object%Wall(1)%node(n_nodes-1,2)
!		ELSEIF (i==n_section-3) THEN
!		    xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!			xyp(2,1)=object%Locations(2,1)
!			xyp(2,2)=ytemp
!	    ELSEIF (i==n_section) THEN
!	        xyp(1,1:2)=object%Wall(1)%node(n_nodes-1,1:2)
!			xyp(2,1:2)=object%Wall(1)%node(n_nodes,1:2)
!		ELSEIF (i<n_section-3) THEN
		    xyp(1,1:2)=object%Wall(1)%node(i,1:2)
			xyp(2,1:2)=object%Wall(1)%node(i+1,1:2)
!		ENDIF
		DO j = 1, 2
			IF (	ALL(xyp(:,j)<=ElBound(j,1)) .OR.	&
					ALL(xyp(:,j)>=ElBound(j,2)) ) THEN
				OutOfBound = .TRUE.
				EXIT 
			END IF
		END DO
		IF (.NOT.OutOfBound) THEN
			CALL Line_Twopoints_Intersection(LineEnds, xyp, inters_point_tmp, inters_flag_tmp)
			IF (inters_flag_tmp/=0) THEN
				count=count+1
				IF (count>1) THEN
!					count2=count2/2
!					count=count-count2
!					IF (count>1) THEN
						PRINT*, 'Case2: Check Partition and Code, Stop'
 						PAUSE
!					ENDIF
				ELSEIF (count==1) THEN
					inters_flag = inters_flag_tmp
					inters_point = inters_point_tmp
					IF (inters_point(1)<0.5) THEN
					print *,'case1',xyp,inters_point,LineEnds
					endif
					EXIT
				ENDIF
			ENDIF

		ENDIF
	ENDDO

!	IF (object%Wall(1)%Shape==1 ) THEN
!		inters_flag_tmp = 0
!		xyp(1,1:2)=object%Wall(1)%node(n_nodes,1:2)
!		xyp(2,1:2)=object%Wall(1)%node(1,1:2)
!		CALL Line_Twopoints_Intersection(LineEnds, xyp, inters_point_tmp, inters_flag_tmp)
!		IF (inters_flag_tmp/=0) THEN
!			count=count+1
!			IF (count>1) THEN
!				PRINT*, 'Case3: Check Partition and Code, Stop',inters_point_tmp
 !				STOP
!			ELSEIF (count==1) THEN
!				inters_flag = inters_flag_tmp
!				inters_point = inters_point_tmp
!					IF (inters_point(1)<0.5) THEN
!					print *,'case2',xyp,inters_point,LineEnds
!					endif
!			ENDIF
!		ENDIF
!	ENDIF

!	IF (object%Wall(1)%Shape==1 ) THEN
!		inters_flag_tmp = 0
!		IF (object%Wall(1)%Channelwall==3 ) THEN
!			CALL Line_Line_Intersection(LineEnds, XBound(2), object%Locations(2,2), inters_point_tmp, inters_flag_tmp)
!		ELSE
!			CALL Line_Line_Intersection(LineEnds, XBound(2), object%Locations(1,2), inters_point_tmp, inters_flag_tmp)
!		ENDIF
!		IF (inters_flag_tmp/=0) THEN
!			count=count+1
!			IF (count>1) THEN
!				PRINT*, 'Case4: Check Partition and Code, Stop',inters_point_tmp
! 				STOP
!			ELSEIF (count==1) THEN
!				inters_flag = inters_flag_tmp
!				inters_point = inters_point_tmp
!!					IF (inters_point(1)<0.5) THEN
!!					print *,'case3',xyp,inters_point,LineEnds
!!					endif
!			ENDIF
!		ENDIF
!	ENDIF


!
   IF (inters_flag==1) THEN
	  AtEnd1 = .False.
	  AtEnd2 = .False.
	  diff1 = (inters_point(1) - x1)*(inters_point(1) - x1) + (inters_point(2) - y1)*(inters_point(2) - y1)
	  diff2 = (inters_point(1) - x2)*(inters_point(1) - x2) + (inters_point(2) - y2)*(inters_point(2) - y2)
	  length = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)
!print*,	inters_point, inters_flag
!print*, x1,y1,x2,y2
!print*, 'sqrt(diff1)=',sqrt(diff1/length), ' sqrt(diff2)=',sqrt(diff2/length), SmallValue
	  IF (sqrt(diff1/length)<=SmallValue)  AtEnd1 = .True.
	  IF (sqrt(diff2/length)<=SmallValue)  AtEnd2 = .True.
	  IF ((.Not. AtEnd1) .AND. (.Not. AtEnd2)) THEN
        P_intrs(e,3)	   = 1
        P_intrs(e,1:2)	   = inters_point
        P_intrs(e,5:6)	   = edges(e,:)
		edge_index(e)	   = 1
	  ENDIF
	  IF ( AtEnd1 .AND. (.Not. AtEnd2)) THEN
        P_intrs(e,4)		= -1
        P_intrs(e,1)	    = x1
        P_intrs(e,2)	    = y1
        P_intrs(e,5)		= edges(e,1)
		vert_index(e)		= 1 
	  ENDIF
	  IF ((.Not. AtEnd1) .AND. AtEnd2) THEN
        P_intrs(e,4)		= -1
        P_intrs(e,1)	    = x2
        P_intrs(e,2)	    = y2
        P_intrs(e,5)		= edges(e,2)
		IF(e/=4)THEN
			vert_index(e+1)	= 1
		ELSE
			vert_index(e)	= 1
		ENDIF		
	  ENDIF
	  IF (AtEnd1 .AND. AtEnd2) THEN
		PRINT*, ' SmallValue is too big, reset it, stop'
		STOP
	  ENDIF
   ELSEIF (inters_flag==-1) THEN
      P_intrs(e,4)		= -1
      P_intrs(e,1:2)	= inters_point
      P_intrs(e,5)		= edges(e,1)
	  vert_index(e)		= 1
   ELSEIF (inters_flag==-2) THEN
      P_intrs(e,4)		= -1
      P_intrs(e,1:2)	= inters_point
      P_intrs(e,5)		= edges(e,2)
      IF(e/=4)THEN
		vert_index(e+1)	= 1
	  ELSE
		vert_index(1)	= 1
	  ENDIF
   END IF

!print*,'e=',e ,P_intrs(e,1:2), P_intrs(e,3)+P_intrs(e,4)
   
END DO


!Check the ill intersection information
nint1 = SUM(ABS(P_intrs(:,3)))
nint2 = SUM(ABS(P_intrs(:,4)))

IF (nint1==1.AND.nint2==3) THEN
  PRINT*, 'nint1=',nint1, ' nint2=',nint2 
  PRINT*, 'Check Partition and Code, Stop'
  PRINT*,'1',P_intrs(1,:)
  PRINT*,'2',P_intrs(2,:)
  PRINT*,'3',P_intrs(3,:)
  STOP
ENDIF

IF (nint1==2.AND.nint2==1) THEN
  PRINT*, 'nint1=',nint1, ' nint2=',nint2
  PRINT*,vert(1,:)
  PRINT*,vert(2,:) 
  PRINT*, 'Check Partition and Code, Stop'
  PRINT*,'1',P_intrs(1,:)
  PRINT*,'2',P_intrs(2,:)
  PRINT*,'3',P_intrs(3,:)
  STOP
ENDIF

IF (nint1==2.AND.nint2==2) THEN
  PRINT*, 'nint1=',nint1, ' nint2=',nint2 
  PRINT*, 'Check Partition and Code, Stop'
  PRINT*,'1',P_intrs(1,:)
  PRINT*,'2',P_intrs(2,:)
  PRINT*,'3',P_intrs(3,:)
  STOP
ENDIF

IF (nint1==3.AND.nint2==1) THEN
  PRINT*, 'nint1=',nint1, ' nint2=',nint2 
  PRINT*, 'Check Partition and Code, Stop'
  PRINT*,'1',P_intrs(1,:)
  PRINT*,'2',P_intrs(2,:)
  PRINT*,'3',P_intrs(3,:)
  STOP
ENDIF

IF (nint1==3.AND.nint2==0) THEN
  PRINT*, 'nint1=',nint1, ' nint2=',nint2 
  PRINT*, 'Check Partition and Code, Stop'
  PRINT*,'1',P_intrs(1,:)
  PRINT*,'2',P_intrs(2,:)
  PRINT*,'3',P_intrs(3,:)
  STOP
ENDIF

IF (nint1==4.AND.nint2==0) THEN
  PRINT*, 'nint1=',nint1, ' nint2=',nint2 
  PRINT*, 'Check Partition and Code, Stop'
  PRINT*,'1',P_intrs(1,:)
  PRINT*,'2',P_intrs(2,:)
  PRINT*,'3',P_intrs(3,:)
  STOP
ENDIF


! Determine regions of physical property (eps) for non-interface elements
!==========================================================CYC modify for generate information==============================
!==========================================================CHJ modify for MORE THAN ONE REGION==============================
!IF(node_index_el(1)==-1 .AND. node_index_el(2)==-1 .AND. node_index_el(3)==-1 .AND. node_index_el(4)==-2) THEN
IF(node_index_el(1)==-1 .AND. node_index_el(2)==-1 .AND. node_index_el(3)==-1  &
        .AND. node_index_el(4)==object%Regions(1)) THEN
!==========================================================CHJ modify for MORE THAN ONE REGION==============================
!3
	el_type = 1
    pointer_local_to_reference = temp_pointer2
    pointer_reference_to_local = temp_pointer4
	plus_piece_flag = 2
    minus_piece_flag = 1
	first_piece_flag = -1
	second_piece_flag = 1
	beta1 = beta_plus
	beta2 = beta_minus
	Dx = P_intrs(3,1)
	Dy = P_intrs(3,2)
	Ex = P_intrs(4,1)
	Ey = P_intrs(4,2)
	nx = (Ey - Dy) / (Ex - Dx)
	ny = -1
	magnitude = SQRT(nx**2 + ny**2)
	nx = nx / magnitude
	ny = ny / magnitude
!==========================================================CHJ modify for MORE THAN ONE REGION==============================
!ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==-1 .AND. node_index_el(3)==-2 .AND. node_index_el(4)==-1) THEN
ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==-1 .AND.   &
    node_index_el(3)==object%Regions(1) .AND. node_index_el(4)==-1) THEN
!==========================================================CHJ modify for MORE THAN ONE REGION==============================
!5
	el_type = 1
    pointer_local_to_reference = temp_pointer3
    pointer_reference_to_local = temp_pointer3
	plus_piece_flag = 2
	minus_piece_flag = 1
	first_piece_flag = -1
	second_piece_flag = 1
	beta1 = beta_plus
	beta2 = beta_minus
	Dx = P_intrs(2,1)
	Dy = P_intrs(2,2)
	Ex = P_intrs(3,1)
	Ey = P_intrs(3,2)
	nx = (Ey - Dy) / (Ex - Dx)
	ny = -1
	magnitude = SQRT(nx**2 + ny**2)
	nx = nx / magnitude
	ny = ny / magnitude
!ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==-1 .AND. node_index_el(3)==-2 .AND. node_index_el(4)==-2) THEN
ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==-1 .AND. node_index_el(3)==object%Regions(1)   &
        .AND. node_index_el(4)==object%Regions(1)) THEN
	IF(edge_index(1)==0 .AND. edge_index(2)==1 .AND. edge_index(3)==0 .AND. edge_index(4)==1) THEN
		!11
		el_type = 2
        pointer_local_to_reference = temp_pointer2
        pointer_reference_to_local = temp_pointer4
        plus_piece_flag = 2
        minus_piece_flag = 1
        first_piece_flag = -1
        second_piece_flag = 1       
        beta1 = beta_plus
        beta2 = beta_minus
		Dx = P_intrs(2,1)
		Dy = P_intrs(2,2)
		Ex = P_intrs(4,1)
		Ey = P_intrs(4,2)
		nx = (Ey - Dy) / (Ex - Dx)
		ny = -1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ENDIF
	IF(vert_index(1)==0 .AND. vert_index(2)==0 .AND. vert_index(3)==1 .AND. vert_index(4)==0) THEN
		!3
		el_type = 1
        pointer_local_to_reference = temp_pointer2
        pointer_reference_to_local = temp_pointer4
		plus_piece_flag = 2
        minus_piece_flag = 1
		first_piece_flag = -1
		second_piece_flag = 1
		beta1 = beta_plus
		beta2 = beta_minus
		Dx = P_intrs(3,1)
		Dy = P_intrs(3,2)
		Ex = P_intrs(4,1)
		Ey = P_intrs(4,2)
		nx = (Ey - Dy) / (Ex - Dx)
		ny = -1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ENDIF
	IF(vert_index(1)==0 .AND. vert_index(2)==0 .AND. vert_index(3)==0 .AND. vert_index(4)==1) THEN
		!5
		el_type = 1
        pointer_local_to_reference = temp_pointer3
        pointer_reference_to_local = temp_pointer3
		plus_piece_flag = 2
        minus_piece_flag = 1
		first_piece_flag = -1
		second_piece_flag = 1
		beta1 = beta_plus
		beta2 = beta_minus
		Dx = P_intrs(2,1)
		Dy = P_intrs(2,2)
		Ex = P_intrs(4,1)
		Ey = P_intrs(4,2)
		nx = (Ey - Dy) / (Ex - Dx)
		ny = -1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ENDIF
!ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==-2 .AND. node_index_el(3)==-1 .AND. node_index_el(4)==-1) THEN
ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==object%Regions(1) .AND. &
    node_index_el(3)==-1 .AND. node_index_el(4)==-1) THEN
	!7
	el_type = 1
    pointer_local_to_reference = temp_pointer4
    pointer_reference_to_local = temp_pointer2
	plus_piece_flag = 2
    minus_piece_flag = 1
	first_piece_flag = -1
	second_piece_flag = 1
	beta1 = beta_plus
	beta2 = beta_minus
	Dx = P_intrs(1,1)
	Dy = P_intrs(1,2)
	Ex = P_intrs(2,1)
	Ey = P_intrs(2,2)
	nx = -(Ey - Dy) / (Ex - Dx)
    ny = 1
    magnitude = SQRT(nx**2 + ny**2)
    nx = nx / magnitude
    ny = ny / magnitude
!ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==-2 .AND. node_index_el(3)==-2 .AND. node_index_el(4)==-1) THEN
ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==object%Regions(1) .AND. &
        node_index_el(3)==object%Regions(1) .AND. node_index_el(4)==-1) THEN
	IF(edge_index(1)==1 .AND. edge_index(2)==0 .AND. edge_index(3)==1 .AND. edge_index(4)==0) THEN
		!10
		el_type = 2
        pointer_local_to_reference = temp_pointer1
        pointer_reference_to_local = temp_pointer1
        plus_piece_flag = 1 
        minus_piece_flag = 2 
        first_piece_flag = 1 
        second_piece_flag = -1                 
        beta1 = beta_minus
        beta2 = beta_plus
		Dx = P_intrs(3,1)
        Dy = P_intrs(3,2)
		Ex = P_intrs(1,1)
		Ey = P_intrs(1,2)

		IF(Dx > Ex)THEN
			nx = -(Ey - Dy) / (Ex - Dx)
            ny = 1
            magnitude = SQRT(nx**2 + ny**2)
            nx = nx / magnitude
            ny = ny / magnitude                
        ELSEIF(Dx < Ex)THEN
            nx = (Ey - Dy) / (Ex - Dx)
            ny = -1
            magnitude = SQRT(nx**2 + ny**2)
            nx = nx / magnitude
            ny = ny / magnitude  
		ELSEIF(Dx == Ex)THEN
            nx = -1
            ny = 0
        ENDIF

	ENDIF
	IF(vert_index(1)==0 .AND. vert_index(2)==1 .AND. vert_index(3)==0 .AND. vert_index(4)==0) THEN
		!5
		el_type = 1
        pointer_local_to_reference = temp_pointer3
        pointer_reference_to_local = temp_pointer3
		plus_piece_flag = 2
        minus_piece_flag = 1
		first_piece_flag = -1
		second_piece_flag = 1
		beta1 = beta_plus
		beta2 = beta_minus
		Dx = P_intrs(2,1)
		Dy = P_intrs(2,2)
		Ex = P_intrs(3,1)
		Ey = P_intrs(3,2)
		nx = (Ey - Dy) / (Ex - Dx)
		ny = -1
		magnitude = DSQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ENDIF
	IF(vert_index(1)==0 .AND. vert_index(2)==0 .AND. vert_index(3)==1 .AND. vert_index(4)==0) THEN
		!7
		el_type = 1
        pointer_local_to_reference = temp_pointer4
        pointer_reference_to_local = temp_pointer2
		plus_piece_flag = 2
        minus_piece_flag = 1
		first_piece_flag = -1
		second_piece_flag = 1
		beta1 = beta_plus
		beta2 = beta_minus
		Dx = P_intrs(1,1)
		Dy = P_intrs(1,2)
		Ex = P_intrs(2,1)
		Ey = P_intrs(2,2)
		nx = -(Ey - Dy) / (Ex - Dx)
		ny = 1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ENDIF
!ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==-2 .AND. node_index_el(3)==-2 .AND. node_index_el(4)==-2) THEN
ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==object%Regions(1) .AND. &
    node_index_el(3)==object%Regions(1) .AND. node_index_el(4)==object%Regions(1)) THEN
	!2
	el_type = 1
    pointer_local_to_reference = temp_pointer1
    pointer_reference_to_local = temp_pointer1
	plus_piece_flag = 1
    minus_piece_flag = 2
	first_piece_flag = 1
	second_piece_flag = -1		
	beta1 = beta_minus
    beta2 = beta_plus
	Dx = P_intrs(4,1)
	Dy = P_intrs(4,2)
	Ex = P_intrs(1,1)
	Ey = P_intrs(1,2)
	nx = (Ey - Dy) / (Ex - Dx)
	ny = -1
	magnitude = SQRT(nx**2 + ny**2)
	nx = nx / magnitude
	ny = ny / magnitude
!ELSEIF(node_index_el(1)==-2 .AND. node_index_el(2)==-1 .AND. node_index_el(3)==-1 .AND. node_index_el(4)==-1) THEN
ELSEIF(node_index_el(1)==object%Regions(1) .AND. node_index_el(2)==-1 .AND. &
        node_index_el(3)==-1 .AND. node_index_el(4)==-1) THEN
	!1
	el_type = 1
    pointer_local_to_reference = temp_pointer1
    pointer_reference_to_local = temp_pointer1
	plus_piece_flag = 2
    minus_piece_flag = 1
	first_piece_flag = -1
	second_piece_flag = 1
	beta1 = beta_plus
	beta2 = beta_minus
	Dx = P_intrs(4,1)
	Dy = P_intrs(4,2)
	Ex = P_intrs(1,1)
	Ey = P_intrs(1,2)
	nx = -(Ey - Dy) / (Ex-Dx)
	ny = 1
	magnitude = SQRT(nx**2 + ny**2)
	nx = nx / magnitude
	ny = ny / magnitude
!ELSEIF(node_index_el(1)==-2 .AND. node_index_el(2)==-1 .AND. node_index_el(3)==-1 .AND. node_index_el(4)==-2) THEN
ELSEIF(node_index_el(1)==object%Regions(1) .AND. node_index_el(2)==-1 .AND. &
        node_index_el(3)==-1 .AND. node_index_el(4)==object%Regions(1)) THEN
	IF(edge_index(1)==1 .AND. edge_index(2)==0 .AND. edge_index(3)==1 .AND. edge_index(4)==0) THEN
		!9
		el_type = 2
        pointer_local_to_reference = temp_pointer1
        pointer_reference_to_local = temp_pointer1
        plus_piece_flag = 2 
        minus_piece_flag = 1 
        first_piece_flag = -1 
        second_piece_flag = 1                 
        beta1 = beta_plus
        beta2 = beta_minus
		Dx = P_intrs(3,1)
		Dy = P_intrs(3,2)
		Ex = P_intrs(1,1)
		Ey = P_intrs(1,2)
		IF(Dx > Ex)THEN
			nx = (Ey - Dy) / (Ex - Dx)
			ny = -1
			magnitude = SQRT(nx**2 + ny**2)
			nx = nx / magnitude
			ny = ny / magnitude                
		ELSEIF(Dx < Ex)THEN
			nx = -(Ey - Dy) / (Ex - Dx)
			ny = 1
			magnitude = SQRT(nx**2 + ny**2)
			nx = nx / magnitude
			ny = ny / magnitude   
        ELSEIF(Dx == Ex)THEN
			nx = 1
			ny = 0
        ENDIF
	ENDIF
	IF(vert_index(1)==0 .AND. vert_index(2)==0 .AND. vert_index(3)==0 .AND. vert_index(4)==1) THEN
		!1
		el_type = 1
        pointer_local_to_reference = temp_pointer1
        pointer_reference_to_local = temp_pointer1
		plus_piece_flag = 2
        minus_piece_flag = 1
		first_piece_flag = -1
		second_piece_flag = 1
		beta1 = beta_plus
		beta2 = beta_minus
		Dx = P_intrs(4,1)
		Dy = P_intrs(4,2)
		Ex = P_intrs(1,1)
		Ey = P_intrs(1,2)
		nx = -(Ey - Dy) / (Ex-Dx)
		ny = 1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ENDIF
	IF(vert_index(1)==1 .AND. vert_index(2)==0 .AND. vert_index(3)==0 .AND. vert_index(4)==0) THEN
		!3
		el_type = 1
        pointer_local_to_reference = temp_pointer2
        pointer_reference_to_local = temp_pointer4
		plus_piece_flag = 2
        minus_piece_flag = 1
		first_piece_flag = -1
		second_piece_flag = 1
		beta1 = beta_plus
		beta2 = beta_minus
		Dx = P_intrs(3,1)
		Dy = P_intrs(3,2)
		Ex = P_intrs(1,1)
		Ey = P_intrs(1,2)
		nx = (Ey - Dy) / (Ex - Dx)
		ny = -1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ENDIF
!ELSEIF(node_index_el(1)==-2 .AND. node_index_el(2)==-1 .AND. node_index_el(3)==-2 .AND. node_index_el(4)==-2) THEN
ELSEIF(node_index_el(1)==object%Regions(1) .AND. node_index_el(2)==-1 .AND. &
        node_index_el(3)==object%Regions(1) .AND. node_index_el(4)==object%Regions(1)) THEN
		IF(nint2==4)THEN
			!3
			el_type = 1
			pointer_local_to_reference = temp_pointer2
			pointer_reference_to_local = temp_pointer4
			plus_piece_flag = 2
			minus_piece_flag = 1
			first_piece_flag = -1
			second_piece_flag = 1
			beta1 = beta_plus
			beta2 = beta_minus
			Dx = P_intrs(3,1)
			Dy = P_intrs(3,2)
			Ex = P_intrs(1,1)
			Ey = P_intrs(1,2)
			nx = (Ey - Dy) / (Ex - Dx)
			ny = -1
			magnitude = SQRT(nx**2 + ny**2)
			nx = nx / magnitude
			ny = ny / magnitude 
		ELSE
			!8
			el_type = 1
			pointer_local_to_reference = temp_pointer4
			pointer_reference_to_local = temp_pointer2
			plus_piece_flag = 1
			minus_piece_flag = 2
			first_piece_flag = 1
			second_piece_flag = -1
			beta1 = beta_minus
			beta2 = beta_plus
			Dx = P_intrs(1,1)
			Dy = P_intrs(1,2)
			Ex = P_intrs(2,1)
			Ey = P_intrs(2,2)
			nx = (Ey - Dy) / (Ex - Dx)
			ny = -1
			magnitude = SQRT(nx**2 + ny**2)
			nx = nx / magnitude
			ny = ny / magnitude
		ENDIF
!ELSEIF(node_index_el(1)==-2 .AND. node_index_el(2)==-2 .AND. node_index_el(3)==-1 .AND. node_index_el(4)==-1) THEN
ELSEIF(node_index_el(1)==object%Regions(1) .AND. node_index_el(2)==object%Regions(1)   &
             .AND. node_index_el(3)==-1 .AND. node_index_el(4)==-1) THEN
	IF(edge_index(1)==0 .AND. edge_index(2)==1 .AND. edge_index(3)==0 .AND. edge_index(4)==1) THEN
		!12
		el_type = 2
        pointer_local_to_reference = temp_pointer2
        pointer_reference_to_local = temp_pointer4
        plus_piece_flag = 1
        minus_piece_flag = 2
        first_piece_flag = 1
        second_piece_flag = -1
        beta1 = beta_minus       
        beta2 = beta_plus
		Dx = P_intrs(2,1)
		Dy = P_intrs(2,2)
		Ex = P_intrs(4,1)
		Ey = P_intrs(4,2)
		nx = -(Ey - Dy) / (Ex - Dx)
		ny = 1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude

	ENDIF
	IF(vert_index(1)==1 .AND. vert_index(2)==0 .AND. vert_index(3)==0 .AND. vert_index(4)==0) THEN
		!7
		el_type = 1
        pointer_local_to_reference = temp_pointer4
        pointer_reference_to_local = temp_pointer2
		plus_piece_flag = 2
        minus_piece_flag = 1
		first_piece_flag = -1
		second_piece_flag = 1
		beta1 = beta_plus
		beta2 = beta_minus
		Dx = P_intrs(1,1)
		Dy = P_intrs(1,2)
		Ex = P_intrs(2,1)
		Ey = P_intrs(2,2)
		nx = -(Ey - Dy) / (Ex - Dx)
		ny = 1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ENDIF
	IF(vert_index(1)==0 .AND. vert_index(2)==1 .AND. vert_index(3)==0 .AND. vert_index(4)==0) THEN
		!1
		el_type = 1
        pointer_local_to_reference = temp_pointer1
        pointer_reference_to_local = temp_pointer1
		plus_piece_flag = 2
        minus_piece_flag = 1
		first_piece_flag = -1
		second_piece_flag = 1
		beta1 = beta_plus
		beta2 = beta_minus
		Dx = P_intrs(4,1)
		Dy = P_intrs(4,2)
		Ex = P_intrs(2,1)
		Ey = P_intrs(2,2)
		nx = -(Ey - Dy) / (Ex-Dx)
		ny = 1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ENDIF
!ELSEIF(node_index_el(1)==-2 .AND. node_index_el(2)==-2 .AND. node_index_el(3)==-1 .AND. node_index_el(4)==-2) THEN
ELSEIF(node_index_el(1)==object%Regions(1) .AND. node_index_el(2)==object%Regions(1)   &
         .AND. node_index_el(3)==-1 .AND. node_index_el(4)==object%Regions(1)) THEN
	IF(nint2==4) THEN
		!1
		el_type = 1
        pointer_local_to_reference = temp_pointer1
        pointer_reference_to_local = temp_pointer1
		plus_piece_flag = 2
	    minus_piece_flag = 1
		first_piece_flag = -1
		second_piece_flag = 1
		beta1 = beta_plus
		beta2 = beta_minus
		Dx = P_intrs(4,1)
		Dy = P_intrs(4,2)	
		Ex = P_intrs(2,1)
		Ey = P_intrs(2,2)
		nx = -(Ey - Dy) / (Ex-Dx)
		ny = 1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude
	ELSE
		!6
		el_type = 1
		pointer_local_to_reference = temp_pointer3
		pointer_reference_to_local = temp_pointer3
		plus_piece_flag = 1
		minus_piece_flag = 2
		first_piece_flag = 1
		second_piece_flag = -1
		beta1 = beta_minus
		beta2 = beta_plus
		Dx = P_intrs(2,1)
		Dy = P_intrs(2,2)
		Ex = P_intrs(3,1)
		Ey = P_intrs(3,2)
		nx = -(Ey - Dy) / (Ex - Dx)
		ny = 1
		magnitude = SQRT(nx**2 + ny**2)
		nx = nx / magnitude
		ny = ny / magnitude;

	ENDIF
!ELSEIF(node_index_el(1)==-2 .AND. node_index_el(2)==-2 .AND. node_index_el(3)==-2 .AND. node_index_el(4)==-1) THEN
ELSEIF(node_index_el(1)==object%Regions(1) .AND. node_index_el(2)==object%Regions(1)   &
         .AND. node_index_el(3)==object%Regions(1) .AND. node_index_el(4)==-1) THEN
	!4
	el_type = 1
    pointer_local_to_reference = temp_pointer2
    pointer_reference_to_local = temp_pointer4
	plus_piece_flag = 1
    minus_piece_flag = 2
	first_piece_flag = 1
	second_piece_flag = -1
	beta1 = beta_minus
	beta2 = beta_plus
	Dx = P_intrs(3,1)
	Dy = P_intrs(3,2)
	Ex = P_intrs(4,1)
	Ey = P_intrs(4,2)
	nx = -(Ey - Dy) / (Ex - Dx)
	ny = 1
	magnitude = SQRT(nx**2 + ny**2)
	nx = nx / magnitude
	ny = ny / magnitude   
ENDIF

!information_1_el(1:4) = t_basic(:,j)
!information_1_el(5) = j
information_1_el(6) = el_type
information_1_el(7:10) = pointer_local_to_reference
information_1_el(11:14) = pointer_reference_to_local
information_1_el(15) = plus_piece_flag
information_1_el(16) = minus_piece_flag
information_1_el(17) = first_piece_flag
information_1_el(18) = second_piece_flag

information_2_el(1) = beta1
information_2_el(2) = beta2
information_2_el(3) = Dx
information_2_el(4) = Dy
information_2_el(5) = Ex
information_2_el(6) = Ey
information_2_el(7) = nx
information_2_el(8) = ny        



!==========================================================CYC modify for generate information==============================



N1 = vert(:,1)
N2 = vert(:,2)
N3 = vert(:,3)
N4 = vert(:,4)

nint1 = SUM(ABS(P_intrs(:,3)))

el_type = 1	! Default Interface element

!print*, 'nint1=',nint1

IF (object%Wall(1)%Shape==2 ) THEN	! Closed curve

	IF (nint1==0) THEN							! Non-interface element
		pm = (N1+N2+N3+N4)/4.0
		el_type = 0
		count=0	
		count2=0
			
		DO i = 1, n_section
			xyp(1,1:2)=object%Wall(1)%node(i,1:2)
			IF (i<n_nodes) THEN
				xyp(2,1:2)=object%Wall(1)%node(i+1,1:2)
			ELSE
				xyp(2,1:2)=object%Wall(1)%node(1,1:2)
			ENDIF
			IF ((pm(2)-xyp(1,2))*(pm(2)-xyp(2,2))<=0) THEN
				xi=(pm(2)-xyp(2,2))*(xyp(2,1)-xyp(1,1))/(xyp(2,2)-xyp(1,2))
				xi=xi+xyp(2,1)
				IF (xi>=pm(1) ) THEN
					count=count+1
					IF (xi==xyp(2,1) .or. xi==xyp(1,1)) THEN
						count2=count2+1
					ENDIF
				ENDIF
			ENDIF
		ENDDO
		count2=count2/2
		count=count-count2
!print*, 'pm=',pm
		IF	(object%Shape==1) THEN
			IF	(object%Dimensions(2)>0) THEN
				IF (MOD(count,2)==0) THEN
				! Inside Sphere
					el_region = object%Regions(1)
				ELSE	!Outside object
					el_region = min(object%Regions(2),el_region)
				ENDIF
			ELSE
				IF (MOD(count,2)==0) THEN		!Outside object
					el_region = min(object%Regions(2),el_region)
				ELSE	! Inside Sphere
					el_region = object%Regions(1)
				ENDIF
			ENDIF
		ELSEIF	(object%Shape==3) THEN		! box
			IF (MOD(count,2)==0) THEN
			! Outside Sphere
				el_region = min(object%Regions(2),el_region)
			ELSE	!Inside object
				el_region = object%Regions(1)
			ENDIF
		ENDIF
	ELSE
	  	    el_type = 1
		    el_region = object%Regions(1)
	END IF

ELSEIF (object%Wall(1)%Shape==1) THEN	! Polygonal line
	
	IF (nint1==0) THEN							! Non-interface element
		pm = (N1+N2+N3+N4)/4.0
		el_type = 0
		count=0	
		count2=0
			
!!		DO i = 1, n_section
!!!			xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!!!			xyp(2,1:2)=object%Wall(1)%node(i+1,1:2)
!!		    IF (i==n_section-1) THEN
!!			    xyp(1,1)=object%Locations(2,1)
!!	            xyp(1,2)=object%Wall(1)%node(n_nodes-1,2)
!!			    xyp(2,1:2)=object%Wall(1)%node(n_nodes-1,1:2)
!!	        ELSEIF (i==n_section-2) THEN
!!		        xyp(1:2,1)=object%Locations(2,1)
!!			    xyp(1,2)=ytemp
!!			    xyp(2,2)=object%Wall(1)%node(n_nodes-1,2)
!!		    ELSEIF (i==n_section-3) THEN
!!		        xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!!			    xyp(2,1)=object%Locations(2,1)
!!			    xyp(2,2)=ytemp
!!	        ELSEIF (i==n_section) THEN
!!	            xyp(1,1:2)=object%Wall(1)%node(n_nodes-1,1:2)
!!			    xyp(2,1:2)=object%Wall(1)%node(n_nodes,1:2)
!!		    ELSEIF (i<n_section-3) THEN
!!		        xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!!			    xyp(2,1:2)=object%Wall(1)%node(i+1,1:2)
!!		    ENDIF
!!			IF ((pm(2)-xyp(1,2))*(pm(2)-xyp(2,2))<=0) THEN
!!				xi=(pm(2)-xyp(2,2))*(xyp(2,1)-xyp(1,1))/(xyp(2,2)-xyp(1,2))
!!				xi=xi+xyp(2,1)
!!				IF (xi>=pm(1) ) THEN
!!					count=count+1
!!					IF (xi==xyp(2,1) .or. xi==xyp(1,1)) THEN
!!						count2=count2+1
!!					ENDIF
!!				ENDIF
!!			ENDIF
!!		ENDDO
!!		count2=count2/2
!!		count=count-count2
!!		
!!		IF (object%Wall(1)%Channelwall==3) THEN
!!			IF ( pm(2)>object%Locations(2,2)) THEN
!!				! Outside Box
!!				el_region = min(object%Regions(2),el_region)
!!			ELSE
!!				IF (MOD(count,2)==0) THEN
!!				! Outside Box
!!					el_region = min(object%Regions(2),el_region)
!!				ELSE
!!					el_region = object%Regions(1)
!!				ENDIF
!!			ENDIF
!!		ELSEIF (object%Wall(1)%Channelwall==4) THEN
!!			IF ( pm(2)<object%Locations(1,2)) THEN
!!				! Outside Box
!!				el_region = min(object%Regions(2),el_region)
!!			ELSE
!!				IF (MOD(count,2)==0) THEN
!!				! Outside Box
!!					el_region = min(object%Regions(2),el_region)
!!				ELSE
!!					el_region = object%Regions(1)
!!				ENDIF
!!			ENDIF
!!		ENDIF

        DO i = 1, n_section
			xyp(1,1:2)=object%Wall(1)%node(i,1:2)
			xyp(2,1:2)=object%Wall(1)%node(i+1,1:2)

            IF (pm(1)>=xyp(1,1) .AND. pm(1)<=xyp(2,1)) THEN
	            ytemp=(xyp(1,2)-xyp(2,2))/(xyp(1,1)-xyp(2,1))
	            ytemp=ytemp*(pm(1)-xyp(2,1))+xyp(2,2)
	            EXIT
            ENDIF
        ENDDO
	    IF( pm(2)> ytemp) THEN
	        ! Outside Object
		    el_region = min(object%Regions(2),el_region)
		ELSE
		    ! Inside Object
		    el_region = object%Regions(1)
        ENDIF    
	        
	ELSE
	  	    el_type = 1
		    el_region = object%Regions(1)
	END IF
ENDIF

!IF (xy(1)/=0 ) THEN
!
!	object%Wall(1)%node(n_nodes-2,1:2)=xy(1:2)
!	object%Wall(1)%node(n_nodes-1,1)=xy(1)
!
!ENDIF

!DEALLOCATE(N1, N2, N3, N4,pm)
!DEALLOCATE(edges)
!DEALLOCATE(Intrs_pts)      !,ElBound
!DEALLOCATE(LineEnds, inters_point,inters_point_tmp)
!DEALLOCATE(nodes)

END