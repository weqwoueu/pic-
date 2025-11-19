SUBROUTINE El_Object_Intersection_info_2D_BL(vert, object, P_intrs, el_type, el_region, &
                                            information_1_el, information_2_el, node_index_el)

! Purpose:		Intersect a generally defined object with a tetrahedal element
! Last Update:	4/21/2004 10:50 PM

USE IFE_MAIN_PARAM
USE Object_Data_2D
IMPLICIT NONE



REAL(8), DIMENSION(2,4), INTENT(IN)		::	vert

TYPE(ObjectType), INTENT(IN)			::	object
REAL(8), DIMENSION(4,6)					::	P_intrs
INTEGER										el_region

REAL(8), DIMENSION(:,:), ALLOCATABLE	::	Intrs_pts, LineEnds
REAL(8), DIMENSION(:), ALLOCATABLE		::	inters_point
INTEGER, DIMENSION(:), POINTER	        ::  information_1_el
REAL(8), DIMENSION(:), POINTER	        ::  information_2_el
INTEGER, DIMENSION(:), INTENT(IN)		::	node_index_el
REAL(8)									::	nx, ny, magnitude

REAL(8) x1,y1, x2,y2, diff1, diff2, length

REAL(8), DIMENSION(:), ALLOCATABLE	::	N1, N2, N3, pm, N4

INTEGER e, inters_flag, nint1, nint2, i, j

INTEGER, ALLOCATABLE, DIMENSION(:,:)	::	edges

REAL(8) Norm_2

LOGICAL AtEnd1, AtEnd2 


!CYC add
INTEGER el_type
INTEGER pointer_local_to_reference(4), pointer_reference_to_local(4)
INTEGER	plus_piece_flag, minus_piece_flag, first_piece_flag, second_piece_flag
INTEGER temp_pointer1(4),temp_pointer2(4),temp_pointer3(4),temp_pointer4(4)
REAL(8)	beta_minus, beta_plus, beta1, beta2
REAL(8) Dx, Ex, Dy, Ey
INTEGER vert_index(4), edge_index(4)
!=========LY modification, 2022-6-13=========
Integer :: L1, L2
Real(8) :: px1, py1, px2, py2
Real(8) :: lx1, ly1, lx2, ly2
Real(8) :: tmp
Real(8) :: x_c, y_c, x_r, y_r
Integer :: flag1, flag2
Real(8) :: cpx1, cpy1, cpx2, cpy2, cpx3, cpy3
Real(8) :: a1, b1, c1, a2, b2, c2
Real(8) :: xr, yr, radius, xr_location(2)
Real(8) :: x, y
Real(8) :: Fx1, Fy1, Fx2, Fy2
Real(8) :: slope, theta, thetaA, thetaB, thetaC
Real(8) :: temp1, temp2, temp3, temp4
!=========LY modification, 2022-6-13=========


ALLOCATE(N1(2), N2(2), N3(2), pm(2), N4(2))

ALLOCATE(edges(4,2))

ALLOCATE(Intrs_pts(2,2))

Intrs_pts	= Zero
P_intrs		= Zero
el_region	= Zero

!--------IFE-----------
!beta_plus = ONE
!beta_minus = 10
!---------------------

!--------PIC---------
beta_plus = object%Eps
beta_minus = ONE
!-----------------------

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


ALLOCATE(LineEnds(2,2), inters_point(2))



DO e=1,4

	x1 = vert(1,edges(e,1))
	y1 = vert(2,edges(e,1))
   
	x2 = vert(1,edges(e,2))
	y2 = vert(2,edges(e,2))
 
	LineEnds(1,:) = (/x1, y1/)
	LineEnds(2,:) = (/x2, y2/)

	inters_flag = 0

    ! Intersection with objects
	IF	(object%Shape==1 .AND. object%Axis==0) THEN
		CALL Line_Circle_Intersection_BL(	LineEnds, object%Dimensions(1), object%Locations(1,1:2), &
										inters_point, inters_flag)
	ENDIF
	IF	(object%Shape==3) THEN  !LY modification, 2022-6-13, replace 'object%Dimensions, object%Locations' with 'object%Dimensions(1:2), object%Locations(1:2,1:2)'
		CALL Line_Rectangular_Intersection_BL(LineEnds, object%Dimensions(1:2), object%Locations(1:2,1:2), &
										inters_point, inters_flag)
  ENDIF
  
  !=========LY modification, 2022-6-13, quadrangle or triangle, ellipse, arch=========
  IF (object%Shape==4) THEN ! yzh
    DO L1=1,4
      IF (L1==4) THEN
        L2=1
      ELSE
        L2=L1+1
      ENDIF
      px1=object%Locations(L1,1)
      px2=object%Locations(L2,1)
      py1=object%Locations(L1,2)
      py2=object%Locations(L2,2)
      lx1=LineEnds(1,1)
      lx2=LineEnds(2,1)
      ly1=LineEnds(1,2)
      ly2=LineEnds(2,2)
      IF (ABS(lx1-lx2)<SmallValue) THEN  !vertical element edge
        IF (ABS(px1-px2)<SmallValue) THEN  !vertical object edge
          inters_flag=inters_flag+0   !'+0' will be protect the past inters_flag
        ELSEIF (ABS(py1-py2)<SmallValue) THEN  !horizon object edge
          IF (py1<MAX(ly1,ly2) .AND. py1>MIN(ly1,ly2) .AND. lx1<=MAX(px1,px2) .AND. lx1>=MIN(px1,px2)) THEN !intersection
            inters_point(1)=lx1
            inters_point(2)=py1
            inters_flag=1
          ENDIF
        ELSE  !inclined object edge
          tmp=(lx1-px2)*(py1-py2)/(px1-px2)+py2 !may be y ,error prone
          IF (tmp-MAX(ly1,ly2)<SmallValue .AND. tmp-MIN(ly1,ly2)>-SmallValue .AND. &
              lx1-MAX(px1,px2)<SmallValue .AND. lx1-MIN(px1,px2)>-SmallValue) THEN  !intersection
            inters_point(1)=Lx1
            inters_point(2)=tmp
            inters_flag=1
          ENDIF
        ENDIF
      ELSEIF(ABS(ly1-ly2)<SmallValue) THEN   !horizon element edge
        IF (ABS(py1-py2)<SmallValue) THEN    !horizon object edge
          inters_flag=inters_flag+0
        ELSEIF(ABS(px1-px2)<SmallValue) THEN !vertical element edge
          IF (px1<MAX(lx1,lx2) .AND. px1>MIN(lx1,lx2) .AND. ly1<=MAX(py1,py2) .AND. ly1>=MIN(py1,py2)) THEN
            inters_point(1)=px1
            inters_point(2)=ly1
            inters_flag=1
          ENDIF
        ELSE  !inclined object edge
          tmp=(ly1-py2)*(px1-px2)/(py1-py2)+px2
          IF(tmp-MAX(lx1,lx2)<SmallValue .AND. tmp-MIN(lx1,lx2)>-SmallValue .AND. &
             ly1-MAX(py1,py2)<SmallValue .AND. ly1-MIN(py1,py2)>-SmallValue) THEN
            inters_point(1)=tmp
            inters_point(2)=ly1
            inters_flag=1
          ENDIF
        ENDIF
      ELSE
        PRINT*, 'Check EL_OBJECT_INTERSECTION_INFO_2D_BL,SHAPE=4,not a true edge'
        STOP
      ENDIF
    ENDDO
  ENDIF !finish
  
  If (object%Shape==5) Then   !ellipse, LY modification for new parameter input way, 2022-5-30
    Call Line_Ellipse_Intersection_BL(LineEnds, object%Dimensions(1:2), object%Locations(1:2,1:2), &
                    inters_point, inters_flag)
  End If
  
  If (object%Shape==6) Then   !Arch
    lx1=LineEnds(1,1)
    lx2=LineEnds(2,1)
    ly1=LineEnds(1,2)
    ly2=LineEnds(2,2)

    !=========LY modification for New parameter input way, 2022-5-31=========    
    cpx1=object%Locations(1,1)  !X coordinate of the first node  A
    cpy1=object%Locations(1,2)  !Y coordinate of the first node  A
    cpx2=object%Locations(2,1)  !X coordinate of the second node B
    cpy2=object%Locations(2,2)  !Y coordinate of the second node B
    cpx3=object%Locations(3,1)  !X coordinate of the third node  C
    cpy3=object%Locations(3,2)  !Y coordinate of the third node  C
    !=========LY modification for New parameter input way, 2022-5-31=========
    
    a1=2*(cpx2-cpx1)
    b1=2*(cpy2-cpy1)
    c1=cpx2*cpx2+cpy2*cpy2-cpx1*cpx1-cpy1*cpy1
    a2=2*(cpx3-cpx2)
    b2=2*(cpy3-cpy2)
    c2=cpx3*cpx3+cpy3*cpy3-cpx2*cpx2-cpy2*cpy2
    
    xr=((c1*b2)-(c2*b1))/((a1*b2)-(a2*b1))      !Centre x of the circle of the arch
    yr=((a1*c2)-(a2*c1))/((a1*b2)-(a2*b1))      !Centre y of the circle of the arch
    radius=SQRT((yr-cpy1)*(yr-cpy1)+(xr-cpx1)*(xr-cpx1))    !Radius of the circle of the arch
    xr_location(1)=xr
    xr_location(2)=yr
    
    CALL Line_Circle_Intersection_BL(	LineEnds, radius, xr_location, &  !intersection with circular curve?
                                      inters_point, inters_flag)
    
    x=inters_point(1)
    y=inters_point(2)
    
    !=========LY modification for New parameter input way, 2022-5-31=========
    If (ABS(cpx1-cpx3) < SmallValue) Then   !cpx1==cpx3: left or right
      If ((cpx2-cpx1) > SmallValue) Then  !cpx2 > cpx1 : right arch, Bx > Ax
        If (x-cpx1 > -SmallValue) Then
          inters_flag = inters_flag + 0
        Else
          inters_flag = 0
        End If
      Elseif ((cpx2-cpx1) < -SmallValue) Then   !cpx2 < cpx1 : left arch, Bx < Ax
        If (cpx1-x > -SmallValue) Then
          inters_flag = inters_flag + 0
        Else
          inters_flag = 0
        End If
      End If
      
    Elseif (ABS(cpy1-cpy3) < SmallValue) Then   !cpy1==cpy3: bottom or top
      If ((cpy2-cpy1) > SmallValue) Then  !cpy2 > cpy1 : top arch, By > Ay
        If (y-cpy1 > -SmallValue) Then
          inters_flag = inters_flag + 0
        Else
          inters_flag = 0
        End If
      Elseif ((cpy2-cpy1) < -SmallValue) Then   !cpy2 < cpy1 : bottom arch, By < Ay
        If (cpy1-y > -SmallValue) Then
          inters_flag = inters_flag + 0
        Else
          inters_flag = 0
        End If
      End If
      
    Else    !cpx1/=cpx3 and cpy1/=cpy3: inclined arch.
      If (inters_flag==0) Then !no intersect weith curve
        lx1=LineEnds(1,1)
        ly1=LineEnds(1,2)
        lx2=LineEnds(2,1)
        ly2=LineEnds(2,2)
        If (ABS(lx1-lx2) < SmallValue) Then  !vertical edge.
          tmp=(lx1-cpx1)*(cpy3-cpy1)/(cpx3-cpx1)+cpy1
          If (lx1-MIN(cpx1,cpx3)>-SmallValue .AND. MAX(cpx1,cpx3)-lx1>-SmallValue .AND. tmp-MIN(ly1,ly2)>-SmallValue .AND.&
              MAX(ly1,ly2)-tmp>-SmallValue) Then
            inters_flag=1
            inters_point(1)=lx1
            inters_point(2)=tmp
          End If
        Elseif (ABS(ly1-ly2) < SmallValue) Then  !horizon edge.
          tmp=(ly1-cpy1)*(cpx3-cpx1)/(cpy3-cpy1)+cpx1
          If (ly1-MIN(cpy1,cpy3)>-SmallValue .AND. MAX(cpy1,cpy3)-ly1>-SmallValue .AND. tmp-MIN(lx1,lx2)>-SmallValue .AND.&
              MAX(lx1,lx2)-tmp>-SmallValue) Then
            inters_flag=1
            inters_point(1)=tmp
            inters_point(2)=ly1
          End If
        End If
        
      Else

        tmp=(y-cpy1)*(cpx3-cpx1)/(cpy3-cpy1)+cpx1
        If (cpx2 - ((cpy2-cpy1)*(cpx3-cpx1)/(cpy3-cpy1)+cpx1) > SmallValue) Then   !right arch.
          If (x-tmp>-SmallValue) Then
            inters_flag=inters_flag+0
          Else
            inters_flag=0
          End If
        Elseif (cpx2 - ((cpy2-cpy1)*(cpx3-cpx1)/(cpy3-cpy1)+cpx1) < -SmallValue) Then   !left arch.
          If (tmp-x>-SmallValue) Then
            inters_flag=inters_flag+0
          Else
            inters_flag=0
          End If
        End If
      End If
    End If
    !=========LY modification for New parameter input way, 2022-5-31=========
  End If
  !=========LY modification, 2022-6-13, quadrangle or triangle, ellipse, arch=========
  
!	IF		(object%Shape==1 .AND. object%Axis==1) THEN
!		! Cylinder // x-axis
!		CALL Line_Cylinder_x_Intersection(	LineEnds, object%Dimensions(1), object%Locations(1,2:3), &
!											object%Locations(1:2,1), inters_point, inters_flag)
!	ELSEIF	(object%Shape==1 .AND. object%Axis==2) THEN
!		! Cylinder // y-axis
!		CALL Line_Cylinder_y_Intersection(	LineEnds, object%Dimensions(1), object%Locations(1,(/1,3/)), &
!											object%Locations(1:2,2), inters_point, inters_flag)
!	ELSEIF	(object%Shape==1 .AND. object%Axis==3) THEN
!		! Cylinder // z-axis
!		CALL Line_Cylinder_z_Intersection(	LineEnds, object%Dimensions(1), object%Locations(1,(/1,2/)), &
!											object%Locations(1:2,3), inters_point, inters_flag)
!	ELSEIF	(object%Shape==1 .AND. object%Axis==4) THEN
!		CALL Line_Cylinder_Intersection(	LineEnds, object%Dimensions(1), object%Locations(1,1:3), &
!											object%Locations(2,1:3), inters_point, inters_flag)
!	ELSEIF	(object%Shape==2 .AND. object%Axis==3) THEN
!!		CALL Line_Sphere_z_Intersection(	LineEnds, object%Radius(1),				&
!											(/object%CenterLine,object%Ends(1)/),	&
!											object%Ends(2), inters_point, inters_flag)
!		CALL Line_Sphere_x_Intersection(	LineEnds, object%Dimensions(1), object%Locations(1,1:3), &
!											object%Locations(2,1), inters_point, inters_flag)
!	ELSEIF	(object%Shape==2 .AND. object%Axis==0) THEN
!!		CALL Line_Sphere_z_Intersection(	LineEnds, object%Radius(1),				&
!!											(/object%CenterLine,object%Ends(1)/),	&
!!											object%Ends(2), inters_point, inters_flag)
!		CALL Line_Sphere_Intersection(	LineEnds, object%Dimensions(1), object%Locations(1,1:3), &
!										inters_point, inters_flag)
!
!
!		!IF (inters_flag/=0) THEN
!		!	PRINT*, SNGL(LineEnds(1,:))
!		!	PRINT*, SNGL(LineEnds(2,:))
!		!	PRINT*, SNGL(inters_point)
!		!END IF
!	END IF

   IF (inters_flag==1) THEN
	  AtEnd1 = .False.
	  AtEnd2 = .False.
	  diff1 = (inters_point(1) - x1)*(inters_point(1) - x1) + (inters_point(2) - y1)*(inters_point(2) - y1)
	  diff2 = (inters_point(1) - x2)*(inters_point(1) - x2) + (inters_point(2) - y2)*(inters_point(2) - y2)
	  length = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)

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
			vert_index(1)	= 1   !LY modification for replace 'vert_index(e)=1' with 'vert_index(1)=1', 2022-6-13
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

   
END DO


!Check the ill intersection information
nint1 = SUM(ABS(P_intrs(:,3)))  !★统计四条边中交点在边上的数目
nint2 = SUM(ABS(P_intrs(:,4)))  !★统计四条边中交点在边顶点上的数目

IF (nint1==1.AND.nint2==3) THEN             !★交点在边上有1条，交点在边顶点上有3条，停止程序并检查
  PRINT*, 'nint1=',nint1, ' nint2=',nint2 
  PRINT*, 'Check Partition and Code, Stop'
  PRINT*,'1',P_intrs(1,:)
  PRINT*,'2',P_intrs(2,:)
  PRINT*,'3',P_intrs(3,:)
  STOP
ENDIF

IF (nint1==2.AND.nint2==1) THEN             !★交点在边上有2条，交点在边顶点上有1条，停止程序并检查
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
!LY modification for Test Multi-Object-Fix-Potential, 2022-7-25, replacing 'node_index_el==-2' with 'node_index_el==object%Regions(1)'
!==========================================================CYC modify for generate information==============================
IF(node_index_el(1)==-1 .AND. node_index_el(2)==-1 .AND. &
   node_index_el(3)==-1 .AND. node_index_el(4)== object%Regions(1)) THEN
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
ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==-1 .AND. &
       node_index_el(3)==object%Regions(1) .AND. node_index_el(4)==-1) THEN
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
ELSEIF(node_index_el(1)==-1 .AND. node_index_el(2)==-1 .AND. &
       node_index_el(3)==object%Regions(1) .AND. node_index_el(4)==object%Regions(1)) THEN
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
ELSEIF(node_index_el(1)==object%Regions(1) .AND. node_index_el(2)==object%Regions(1) .AND. &
       node_index_el(3)==-1 .AND. node_index_el(4)==-1) THEN
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
ELSEIF(node_index_el(1)==object%Regions(1) .AND. node_index_el(2)==object%Regions(1) .AND. &
       node_index_el(3)==-1 .AND. node_index_el(4)==object%Regions(1)) THEN
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
ELSEIF(node_index_el(1)==object%Regions(1) .AND. node_index_el(2)==object%Regions(1) .AND. &
       node_index_el(3)==object%Regions(1) .AND. node_index_el(4)==-1) THEN
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

IF (object%Shape==1 .AND. object%Axis==0) THEN	! Circle

	IF (nint1==0) THEN							! Non-interface element
		pm = (N1+N2+N3)/3.0
		el_type = 0

!print*, 'pm=',pm

		IF ((Norm_2(N1-object%Locations(1,:))<=object%Dimensions(1)).AND.&
		    (Norm_2(N2-object%Locations(1,:))<=object%Dimensions(1)).AND.&
		    (Norm_2(N3-object%Locations(1,:))<=object%Dimensions(1)).AND.&
		    (Norm_2(N4-object%Locations(1,:))<=object%Dimensions(1))) THEN
			! Inside Sphere
			el_region = object%Regions(1)
		ELSE
			! Outside Sphere
			el_region = object%Regions(2)
		ENDIF
	ELSE
	  	    el_type = 1
		    el_region = object%Regions(1)
	END IF

ELSEIF (object%Shape==3) THEN	! Box

	IF (nint1==0) THEN							! Non-interface element

		pm = (N1+N2+N3+N4)/4.0
		el_type = 0

		IF (pm(1)>object%Locations(1,1)	.AND. pm(1)<object%Locations(2,1)	.AND.	& 
				pm(2)>object%Locations(1,2)	.AND. pm(2)<object%Locations(2,2)) THEN
			! Inside Box
			el_region = object%Regions(1)
		ELSE
			! Outside Box
			el_region = object%Regions(2)
		ENDIF
  END IF
  
ELSEIF (object%Shape==4 .OR. object%Shape==5 .OR. object%Shape==6) THEN  !LY modification, 2022-6-1, quadrangle or triangle, ellipse and arch
  
  !=========LY modification, 2022-6-2=========  
  IF (DABS(P_intrs(1,1)-P_intrs(2,1))<SmallValue .AND. DABS(P_intrs(1,2)-P_intrs(2,2))<SmallValue .AND.&
          (P_intrs(1,1)/=0 .OR. P_intrs(1,2)/=0)) THEN
    nint2=nint2-1
  ENDIF
  IF (DABS(P_intrs(2,1)-P_intrs(3,1))<SmallValue .AND. DABS(P_intrs(2,2)-P_intrs(3,2))<SmallValue .AND.&
          (P_intrs(2,1)/=0 .OR. P_intrs(2,2)/=0)) THEN
    nint2=nint2-1
  ENDIF
  IF (DABS(P_intrs(3,1)-P_intrs(4,1))<SmallValue .AND. DABS(P_intrs(3,2)-P_intrs(4,2))<SmallValue .AND.&
          (P_intrs(3,1)/=0 .OR. P_intrs(3,2)/=0)) THEN
    nint2=nint2-1
  ENDIF
  IF (DABS(P_intrs(4,1)-P_intrs(1,1))<SmallValue .AND. DABS(P_intrs(4,2)-P_intrs(1,2))<SmallValue .AND.&
          (P_intrs(4,1)/=0 .OR. P_intrs(4,2)/=0)) THEN
    nint2=nint2-1
  ENDIF
  !=========LY modification, 2022-6-2=========
  
  !=========YZH edition=========
  !IF (DABS(P_intrs(1,1)-P_intrs(2,1))<SmallValue .AND. DABS(P_intrs(1,2)-P_intrs(2,2))<SmallValue .AND.&
  !    P_intrs(1,1)+P_intrs(1,2)/=0) THEN
  !    nint2=nint2-1
  !ENDIF
  !IF (DABS(P_intrs(2,1)-P_intrs(3,1))<SmallValue .AND. DABS(P_intrs(2,2)-P_intrs(3,2))<SmallValue .AND.&
  !    P_intrs(2,1)+P_intrs(2,2)/=0) THEN
  !    nint2=nint2-1
  !ENDIF
  !IF (DABS(P_intrs(3,1)-P_intrs(4,1))<SmallValue .AND. DABS(P_intrs(3,2)-P_intrs(4,2))<SmallValue .AND.&
  !    P_intrs(3,1)+P_intrs(3,2)/=0) THEN
  !    nint2=nint2-1
  !ENDIF
  !IF (DABS(P_intrs(4,1)-P_intrs(1,1))<SmallValue .AND. DABS(P_intrs(4,2)-P_intrs(1,2))<SmallValue .AND.&
  !    P_intrs(4,1)+P_intrs(4,2)/=0) THEN
  !    nint2=nint2-1
  !ENDIF
  !=========YZH edition=========

  IF (nint1+nint2==0 .OR. nint1+nint2==1) THEN							! Non-interface element

    !only four points are all -2,then it is a inside element

    el_type = 0

    !IF (node_index_el(1)+node_index_el(2)+node_index_el(3)+node_index_el(4)==object%Regions(1)*4) THEN
    IF (node_index_el(1)<=object%Regions(1) .AND. node_index_el(2)<=object%Regions(1) .AND.&
        node_index_el(3)<=object%Regions(1) .AND. node_index_el(4)<=object%Regions(1)) THEN   !LY modifiction for Test Multi-Object-Fix-Potential, 2022-7-25
      ! Inside Sphere
      el_region = object%Regions(1)
    ELSE
      ! Outside Sphere
      el_region = object%Regions(2)
    ENDIF

    information_1_el(6)=0 !yzh ,for pingjie judge
    
  END IF  !finish
  
ENDIF


!
!IF (object%Shape==1 .AND. object%Axis>=1 .AND. object%Axis<=3)	THEN	! Cylinder // one of the axes
!
!	IF		(object%Axis==1) THEN
!		coord1 = 1
!		coord2 = 2
!		coord3 = 3
!	ELSEIF	(object%Axis==2) THEN
!		coord1 = 2
!		coord2 = 3
!		coord3 = 1
!	ELSEIF	(object%Axis==3) THEN
!		coord1 = 3
!		coord2 = 1
!		coord3 = 2
!	ENDIF
!
!
!	IF (nint1==0) THEN							! Non-interface element
!
!		pm = (N1+N2+N3+N4)/4
!		el_type = 0
!		IF (	pm(coord1)>object%Locations(1,coord1) .AND.	& 
!				pm(coord1)<object%Locations(2,coord1) .AND.	&
!				Norm_3(pm((/coord2,coord3/))-object%Locations(1,(/coord2,coord3/)))<object%Dimensions(1) ) THEN
!			! Inside Cylinder
!			el_region = object%Regions(1)
!		ELSE
!			! Outside Cylinder
!			el_region = object%Regions(2)
!		ENDIF
!		
!
!	END IF
!
!ELSEIF (object%Shape==1 .AND. object%Axis>3) THEN
!
!	r  = object%Dimensions(1)
!	Xb = object%Locations(1,:)
!	Xt = object%Locations(2,:)
!
!	! From here till the next ELSEIF the code works for a general cylinder
!	Xtb = Xt - Xb
!	CALL Vec_Dot_3(Xtb, Xtb, XtbDotXtb)
!
!	IF (nint1==0) THEN							! Non-interface element
!
!		pm  = (N1+N2+N3+N4)/4
!		el_type = 0
!
!		Xp  = pm
!		Xpb = Xp - Xb
!		CALL Vec_Dot_3(Xpb, Xtb, XpbDotXtb)
!		muc = -XpbDotXtb/XtbDotXtb
!		Xc  = Xb + muc*(Xt - Xb)
!		Xpt = Xp - Xt
!		CALL Vec_Dot_3(Xpt, Xtb, XptDotXtb)
!		Xpc = Xp - Xc
!		CALL Vec_Dot_3(Xpc, Xpc, XpcDotXpc)
!
!		IF (XpbDotXtb>=-MZero .AND. XptDotXtb<=MZero .AND. XpcDotXpc<=r**2) THEN
!			! Inside the object
!			el_region = object%Regions(1)
!		ELSE
!			! Outside Cylinder
!			el_region = object%Regions(2)
!		ENDIF
!		
!
!	END IF
!
!ELSEIF (object%Shape==2 .AND. object%Axis>0) THEN	! Sphere
!
!	IF		(object%Axis==1) THEN
!		coord1 = 1
!		coord2 = 2
!		coord3 = 3
!	ELSEIF	(object%Axis==2) THEN
!		coord1 = 2
!		coord2 = 3
!		coord3 = 1
!	ELSEIF	(object%Axis==3) THEN
!		coord1 = 3
!		coord2 = 1
!		coord3 = 2
!	ENDIF
!
!	IF (nint1==0) THEN							! Non-interface element
!
!		pm = (N1+N2+N3+N4)/4
!		el_type = 0
!
!		IF (	pm(coord1)>object%Locations(2,1)							.AND.	& 
!				pm(coord1)<object%Locations(1,coord1)+object%Dimensions(1)	.AND.	&
!				Norm_3( pm-object%Locations(1,:))<object%Dimensions(1) ) THEN
!			! Inside Sphere
!			el_region = object%Regions(1)
!		ELSE
!			! Outside Sphere
!			el_region = object%Regions(2)
!		ENDIF
!		
!
!	END IF
!
!ELSEIF (object%Shape==2 .AND. object%Axis==0) THEN	! Sphere
!
!	IF (nint1==0) THEN							! Non-interface element
!
!		pm = (N1+N2+N3+N4)/4
!		el_type = 0
!
!		IF (	Norm_3( pm-object%Locations(1,:))<object%Dimensions(1) ) THEN
!			! Inside Sphere
!			el_region = object%Regions(1)
!		ELSE
!			! Outside Sphere
!			el_region = object%Regions(2)
!		ENDIF
!		
!
!	END IF
!
!ELSEIF (object%Shape==3) THEN	! Box
!
!	IF (nint1==0) THEN							! Non-interface element
!
!		pm = (N1+N2+N3+N4)/4
!		el_type = 0
!
!		IF (	pm(1)>object%Locations(1,1)	.AND. pm(1)<object%Locations(2,1)	.AND.	& 
!				pm(2)>object%Locations(1,2)	.AND. pm(2)<object%Locations(2,2)	.AND.	&
!				pm(3)>object%Locations(1,3)	.AND. pm(3)<object%Locations(2,3)	) THEN
!			! Inside Sphere
!			el_region = object%Regions(1)
!		ELSE
!			! Outside Sphere
!			el_region = object%Regions(2)
!		ENDIF
!
!		
!
!	END IF
!
!ENDIF
!
!print*,Intrs_pts
!stop

DEALLOCATE(N1, N2, N3, pm, N4)
DEALLOCATE(edges)
DEALLOCATE(Intrs_pts)
DEALLOCATE(LineEnds, inters_point)


END