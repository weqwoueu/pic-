SUBROUTINE Update_Locate_Node_2D	(	p_node, object, node_region,Phi_p	)

! Purpose:		Locate node with respect to an immersed element
! Last Update:	4/21/2004 07:09 PM

USE IFE_MAIN_PARAM
USE Object_Data_2D
!USE Wall_2D    !LY modification, 2022-6-13

IMPLICIT NONE


REAL(8), DIMENSION(2), INTENT(IN)	::	p_node
TYPE(ObjectType), INTENT(IN)		::	object
INTEGER									i,node_region, coord1, coord2		
INTEGER									n_section,n_nodes,count,count2

REAL(8), DIMENSION(2)				::	X1, X2, Xb, Xt, Xc, Xtb, Xpb, Xpt, Xpc
REAL(8)									y_node,Phi_p,xi,ytemp
REAL(8), DIMENSION(2,2)				::	xyp

REAL(8) Norm_2 
!=========LY modification, 2022-6-13=========
Real(8) :: x, y
Real(8) :: ax, ay, bx, by, cx, cy, dx, dy
Real(8) :: S123, S134, S234, S124
Real(8) :: Spab, Spbc, Spcd, Spda, Spbd, Sq, Sp
Real(8) :: mark_abc, mark_bcd, mark_cda, mark_dab
Real(8) :: x_c, y_c, x_r, y_r
Real(8) :: x_a, y_a, x_b, y_b
Real(8) :: Larch, harch, Rarch
Real(8) :: node1(3), node2(3), node3(3)
Real(8) :: center_arch(3), Cz
Real(8) :: Mxarch, Miarch
Real(8) :: yL
Integer :: flag1, flag2
Real(8) :: cpx1, cpy1, cpx2, cpy2, cpx3, cpy3
Real(8) :: a1, b1, c1, a2, b2, c2, xr, yr, radius
Real(8) :: tmp
Real(8) :: Fx1, Fy1, Fx2, Fy2
Real(8) :: vecACx, vecACy, vecAZx, vecAZy, vecANx, vecANy
Real(8) :: vecAN_len, vecAC_len, vecAZ_len
Real(8) :: vecAN_dot_vecAC, vecAN_dot_vecAZ
Real(8) :: slope, theta, thetaA, thetaB, thetaC
Real(8) :: temp
!=========LY modification, 2022-6-13=========

IF		(object%Axis==1) THEN
		coord1 = 1
		coord2 = 2
ELSEIF	(object%Axis==2) THEN
		coord1 = 2
		coord2 = 1
ENDIF

IF (object%Erosion==0) THEN
	IF (object%Shape==1 .AND. object%Axis==0) THEN	! Circle

		IF	(object%Dimensions(2)>0) THEN
			IF (	Norm_2( p_node-object%Locations(1,:))<object%Dimensions(2) ) THEN
        !=========LY modification, 2022-6-13=========
        !Inside Sphere
        node_region = object%Regions(1)
				!Phi_p=0
        !=========LY modification, 2022-6-13=========
      ELSE
        !=========LY modification, 2022-6-13=========
        !Outside Sphere
        node_region = min(object%Regions(2),node_region)
				!Phi_p=object%Phi
        !=========LY modification, 2022-6-13=========
			ENDIF	
		ELSE
			IF (	Norm_2( p_node-object%Locations(1,:))<=object%Dimensions(1) ) THEN
			! Inside Sphere
				node_region = object%Regions(1)
				!Phi_p=object%Phi
			ELSE
			! Outside Sphere
				node_region = min(object%Regions(2),node_region)
        !!!!*************** bjw add 2019.9.12 ****************
        !Phi_p=object%Phi
				!Phi_p=0
			ENDIF
		ENDIF	

	ELSEIF (object%Shape==3) THEN	! Box
	
		IF (p_node(1)>(object%Locations(1,1)-SmallValue).AND. p_node(1)<(object%Locations(2,1)+SmallValue).AND.	& 
				p_node(2)>(object%Locations(1,2)-SmallValue).AND. p_node(2)<(object%Locations(2,2)+SmallValue)) THEN
			! Inside 
			node_region = object%Regions(1)
			!Phi_p=object%Phi
		ELSE
			! Outside 
			node_region = min(object%Regions(2),node_region)
!			Phi_p=0
    ENDIF
    
  Elseif (object%Shape==4) Then  !LY modification for quadrangle and triangle, 2022-6-13
    
    x = p_node(1)
    y = p_node(2)
    ax = object%Locations(1,1)
    ay = object%Locations(1,2)
    bx = object%Locations(2,1)
    by = object%Locations(2,2)
    cx = object%Locations(3,1)
    cy = object%Locations(3,2)
    dx = object%Locations(4,1)
    dy = object%Locations(4,2)
    
    !Firstly, we need to calculate some triangle area.
    !Note: In Concave quadrilateral(or ao quadrilateral), the inner node is 1 node, the bottom nodes is 2 and 4 node, the top node is 3 node.
    S123 = DABS(0.5*(ax*by + ay*cx + bx*cy - cx*by - cy*ax - bx*ay))  !left triangle of tu quadrilateral
	  S134 = DABS(0.5*(ax*cy + ay*dx + cx*dy - dx*cy - dy*ax - cx*ay))  !right triangle of tu quadrilateral
	  S234 = DABS(0.5*((cx - bx)*(dy - by) - (dx - bx)*(cy - by)))  !big triangle of ao quadrilateral
	  S124 = DABS(0.5*((bx - ax)*(dy - ay) - (dx - ax)*(by - ay)))  !small triangle of ao quadrilateral(it will be minus)
    
    !the triangle area between of mesh node and tu quadrilateral node.
    Spab = DABS(0.5*((ax - x)*(by - y) - (bx - x)*(ay - y)))  !P-ab triangle
	  Spbc = DABS(0.5*((bx - x)*(cy - y) - (cx - x)*(by - y)))  !P-bc triangle
	  Spcd = DABS(0.5*((cx - x)*(dy - y) - (dx - x)*(cy - y)))  !P-cd triangle
	  Spda = DABS(0.5*((dx - x)*(ay - y) - (ax - x)*(dy - y)))  !P-da triangle
    
    !the triangle area between of mesh node and ao quadrilateral bottom node(2 and 4 node).
    Spbd = DABS(0.5*((bx - x)*(dy - y) - (dx - x)*(by - y)))  !P-bd triangle
    
    !the area of tu quadrilateral.
    Sq = S123 + S134  !Sq=S123+S134
    
    !the area of four triangle from mesh node and tu quadrilateral node.
    Sp = Spab + Spbc + Spcd + Spda
    
    !Secondly, we need to judge the type of quadrilateral: tu or ao?
    mark_abc = ax*by + ay*cx + bx*cy - cx*by - cy*ax - bx*ay  !abc vector product
	  mark_bcd = bx*cy + by*dx + cx*dy - dx*cy - dy*bx - cx*by  !bcd vector product
	  mark_cda = cx*dy + cy*ax + dx*ay - ax*dy - ay*cx - dx*cy  !cda vector product
	  mark_dab = dx*ay + dy*bx + ax*by - bx*ay - by*dx - ax*dy  !dab vector product
    
    !Finally, we can judge the type of quadrilateral and the node location.
    IF ((mark_abc > -SmallValue .AND. mark_bcd > -SmallValue .AND. mark_cda > -SmallValue .AND. mark_dab > -SmallValue) .OR. &
        (mark_abc < SmallValue .AND. mark_bcd < SmallValue .AND. mark_cda < SmallValue .AND. mark_dab < SmallValue)) THEN   !tu quadrilateral
      IF (DABS(Sp - Sq) < SmallValue) THEN  
        node_region = object%Regions(1)   !node locate inside tu quadrilateral
      ELSE 
        node_region = object%Regions(2)   !node locate outside tu quadrilateral
      ENDIF

    ELSE  !ao quadrilateral
      IF (DABS(Spbd + Spbc + Spcd - S234)<SmallValue) THEN  !node locate inside big triangle of ao quadrilateral
        IF ((DABS(Spbd + Spbc + Spcd - S234)<SmallValue) .AND. (DABS(Spda) > SmallValue) .AND. (DABS(Spab) > SmallValue)) THEN
          node_region = object%Regions(2)  !node locate inside small triangle of ao quadrilateral.
        ELSE 
          node_region = object%Regions(1)  !node locate inside the triangle between big triangle minus small triangle of ao quadrilateral.
        ENDIF
      ELSE
        node_region = object%Regions(2)   !node locate outside big triangle of ao quadrilateral
      ENDIF
    ENDIF
    !=========LY modification, 2022-6-6=========
    !In Triangle or Quardangle, we consider that the endpoints and nodes on Triangle or Quardangle line is inner node.
    !=========LY modification, 2022-6-6=========
    
  Elseif (object%Shape==5) Then  !LY modification for ellipse, 2022-6-13
    
    x = p_node(1)
    y = p_node(2)
    !=========LY modification, 2022-5-30=========
    Fx1 = object%Locations(1,1)   !Focus1 x of the ellipse
    Fy1 = object%Locations(1,2)   !Focus1 y of the ellipse
    Fx2 = object%Locations(2,1)   !Focus2 x of the ellipse
    Fy2 = object%Locations(2,2)   !Focus2 y of the ellipse
    
    x_c = (Fx1+Fx2)/2.0           !Centre x of the ellipse
    y_c = (Fy1+Fy2)/2.0           !Centre y of the ellipse
    x_r = object%Dimensions(1)    !Radius x of the ellipse
    y_r = object%Dimensions(2)    !Radius y of the ellipse
    !=========LY modification, 2022-5-30=========
    
    !=========LY modification, 2022-6-1=========
    If (ABS(Fy1-Fy2) < 1.0D-8) Then    !horizontal ellipse, focus on x axis.
      If ((x-x_c)*(x-x_c)/(x_r*x_r) + (y-y_c)*(y-y_c)/(y_r*y_r) - 1.0 < -SmallValue) Then
        !Inside Sphere
        node_region = object%Regions(1)
      Else
        !Outside Sphere
        node_region = min(object%Regions(2), node_region)
      End If
      
    Elseif (ABS(Fx1-Fx2) < 1.0D-8) Then    !vertical ellipse, focus on y axis.
      If ((x-x_c)*(x-x_c)/(x_r*x_r) + (y-y_c)*(y-y_c)/(y_r*y_r) - 1.0 < -SmallValue) Then
        !Inside Sphere
        node_region = object%Regions(1)
      Else
        !Outside Sphere
        node_region = min(object%Regions(2), node_region)
      End If
      
    Else    !inclined ellipse
      slope = (Fy2-Fy1) / (Fx2-Fx1)
      If (ABS(slope) < 1.0D-8) Then   !There denote slope=0, but we already deal with the condition in above.
        Write(*,*) 'The slope = ', slope
        Stop
      End If
      
      If (slope > 1.0D-8) Then
        theta = DATAN(slope)
      Elseif (slope < -1.0D-8) Then
        theta = DATAN(slope) + pi
      End If
      thetaA = (((DCOS(theta))**2)/(x_r**2) + ((DSIN(theta))**2)/(y_r**2))
      thetaB = (1/(x_r**2) - 1/(y_r**2)) * DSIN(2*theta)
      thetaC = (((DSIN(theta))**2)/(x_r**2) + ((DCOS(theta))**2)/(y_r**2))
      
      temp = thetaA*x*x + thetaB*x*y + thetaC*y*y - (2*thetaA*x_c+thetaB*y_c)*x - (thetaB*x_c+2*thetaC*y_c)*y + &
              (thetaA*x_c*x_c+thetaB*x_c*y_c+thetaC*y_c*y_c)
      
      If (temp - 1.0 < -1.0D-8) Then
        !Inside Sphere
        node_region = object%Regions(1)
      Else
        !Outside Sphere
        node_region = min(object%Regions(2), node_region)
      End If
    End If
    !=========LY modification, 2022-6-1=========
    !=========LY modification, 2022-6-6=========
    !In Ellipse, we consider that the nodes on Ellipse curve is outer node.
    !=========LY modification, 2022-6-6=========
    
  Elseif (object%Shape==6) Then  !LY modification for arch, 2022-6-13

    !=========LY modification for Yao edition, 2022-5-12=========
    !=========LY modification for New parameter input way, 2022-5-31=========
    cpx1=object%Locations(1,1)  !X coordinate of the first node  A
    cpy1=object%Locations(1,2)  !Y coordinate of the first node  A
    cpx2=object%Locations(2,1)  !X coordinate of the second node B
    cpy2=object%Locations(2,2)  !Y coordinate of the second node B
    cpx3=object%Locations(3,1)  !X coordinate of the third node  C
    cpy3=object%Locations(3,2)  !Y coordinate of the third node  C
    
    vecACx = cpx3 - cpx1    !vector AC x component
    vecACy = cpy3 - cpy1    !vector AC y component
    vecAZx = cpy3 - cpy1    !vector AZ x component
    vecAZy = cpx1 - cpx3    !vector AZ y component
    !=========LY modification for New parameter input way, 2022-5-31=========
    
    a1=2*(cpx2-cpx1)
    b1=2*(cpy2-cpy1)
    c1=cpx2*cpx2+cpy2*cpy2-cpx1*cpx1-cpy1*cpy1
    a2=2*(cpx3-cpx2)
    b2=2*(cpy3-cpy2)
    c2=cpx3*cpx3+cpy3*cpy3-cpx2*cpx2-cpy2*cpy2

    xr=((c1*b2)-(c2*b1))/((a1*b2)-(a2*b1))
    yr=((a1*c2)-(a2*c1))/((a1*b2)-(a2*b1))
    radius=DSQRT((yr-cpy1)*(yr-cpy1)+(xr-cpx1)*(xr-cpx1))

    x=p_node(1)
    y=p_node(2)

    IF(DSQRT((x-xr)*(x-xr)+(y-yr)*(y-yr))-radius<SmallValue)THEN

      !=========LY modification for New parameter input way, 2022-5-31=========
      vecANx = x - cpx1   !vector AN x component
      vecANy = y - cpy1   !vector AN y component
      vecAN_len = DSQRT(vecANx**2 + vecANy**2)    !vector AN length
      vecAC_len = DSQRT(vecACx**2 + vecACy**2)    !vector AC length
      vecAZ_len = DSQRT(vecAZx**2 + vecAZy**2)    !vector AZ length
      vecAN_dot_vecAC = vecANx * vecACx + vecANy * vecACy   !dot product between vector AN and AC
      vecAN_dot_vecAZ = vecANx * vecAZx + vecANy * vecAZy   !dot product between vector AN and AZ

      !Node locate inside arch only cos¦È1 = (AN¡¤AC)/(|AN|*|AC|)>0 and cos¦È2 = (AN¡¤AZ)/(|AN|*|AZ|)>0.
      !At same time, we consider node locate inside arch when node is A, B and C.
      If (vecAN_len < SmallValue .OR. ABS(vecAN_dot_vecAZ) < SmallValue) Then
        !|AN|=0, that is the node N is node A. We consider the node locate inside arch.
        !AN¡¤AZ=0, that is the node N inside AC. We consider the node locate inside arch.
        node_region = object%Regions(1)  
      Else
        If ((vecAN_dot_vecAC/vecAN_len/vecAC_len > SmallValue) .AND. &
            (vecAN_dot_vecAZ/vecAN_len/vecAZ_len > SmallValue)) Then
          node_region = object%Regions(1)    !Node locate inside arch.
        Else
          node_region = object%Regions(2)    !Node locate outside arch.
        End If
      End If
      !=========LY modification for New parameter input way, 2022-5-31=========

    ELSE
      node_region = object%Regions(2) !outside
    ENDIF
    !=========LY modification for Yao edition, 2022-5-12=========
    !=========LY modification, 2022-6-6=========
    !In here, we consider that the endpoints and node on arch line and curve is inner node.
    !=========LY modification, 2022-6-6=========
    
	ELSEIF (object%Shape==11) THEN	! Line, Plate parallel to X axis(delta ==0) or Z axis(delta ==1)

			X1 = object%Locations(1,:)
			X2 = object%Locations(2,:)

!		IF (	(X2(1)-p_node(1))*(p_node(1)-X1(1)) >= -SmallValue						.AND.	&
!				(X2(2)-p_node(2))*(p_node(2)-X1(2)) >= -SmallValue						.AND.	&
!				DABS((p_node(coord2)-X1(coord2))*(X2(coord1)-X1(coord1)))<= SmallValue	) THEN
		IF (	(X2(1)-p_node(1))*(p_node(1)-X1(1)) >= 0						.AND.	&
				(X2(2)-p_node(2))*(p_node(2)-X1(2)) >= 0						.AND.	&
				DABS((p_node(coord2)-X1(coord2))*(X2(coord1)-X1(coord1)))<= 0	) THEN
			! Inside 
			node_region = object%Regions(1)
			Phi_p=object%Phi

		ELSE
			! Outside 
			node_region = object%Regions(2)
			Phi_p=0
		ENDIF
	
	ENDIF
ELSE
  !=========LY modification, 2022-6-13=========
  !There denote we use erosion, error.
  WRITE(*,*) 'object%Erosion = ', object%Erosion
  STOP
  !=========LY modification, 2022-6-13=========
!	n_nodes	=	object%Wall(1)%nxp
!
!	IF (object%Wall(1)%Shape==2) THEN
!		n_section	=	n_nodes	
!	ELSE
!	        DO i=1,n_nodes      !object%Wall(1)%nxp
!			  IF (object%Wall(1)%node(i,1)>=object%Locations(2,1)) THEN
!				n_section=i
!				IF (object%Wall(1)%node(i,1)==object%Locations(2,1)) THEN
!					ytemp=object%Wall(1)%node(i,2)
!				ELSEIF (object%Wall(1)%node(i,1)>object%Locations(2,1)) THEN
!					xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!					ytemp=(xyp(1,2)-object%Wall(1)%node(i-1,2))/(xyp(1,1)-object%Wall(1)%node(i-1,1))
!					ytemp=ytemp*(object%Locations(2,1)-xyp(1,1))+xyp(1,2)
!				ENDIF
!!				print *,n_section,object%Wall(1)%node(i,1),ytemp
!				EXIT
!			  ENDIF
!			ENDDO
!			
!		n_section	=	n_section-1     !n_nodes-1
!	ENDIF
!
!	IF (object%Wall(1)%Shape==2 ) THEN	! Closed curve
!
!		count=0	
!		count2=0
!			
!		DO i = 1, n_section
!			xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!			IF (i<n_nodes) THEN
!				xyp(2,1:2)=object%Wall(1)%node(i+1,1:2)
!			ELSE
!				xyp(2,1:2)=object%Wall(1)%node(1,1:2)
!			ENDIF
!			IF ((p_node(2)-xyp(1,2))*(p_node(2)-xyp(2,2))<=0) THEN
!				xi=(p_node(2)-xyp(2,2))*(xyp(2,1)-xyp(1,1))/(xyp(2,2)-xyp(1,2))
!				xi=xi+xyp(2,1)
!				IF (xi>=p_node(1) ) THEN
!					count=count+1
!					IF (xi==xyp(2,1) .or. xi==xyp(1,1)) THEN
!						count2=count2+1
!					ENDIF
!				ENDIF
!			ENDIF
!		ENDDO
!		count2=count2/2
!		count=count-count2
!	!print*, 'pm=',pm
!		IF	(object%Shape==1) THEN
!			IF	(object%Dimensions(2)>0) THEN
!				IF (MOD(count,2)==0) THEN
!					! Inside Sphere
!					node_region = object%Regions(1)
!					Phi_p=object%Phi
!				ELSE
!					! Outside Sphere
!					node_region = min(object%Regions(2),node_region)
!!					Phi_p=0
!				ENDIF
!			ELSE
!				IF (MOD(count,2)==0) THEN		!Outside object
!					node_region = min(object%Regions(2),node_region)
!!					Phi_p=0
!				ELSE	! Inside Sphere
!					node_region = object%Regions(1)
!					Phi_p=object%Phi
!				ENDIF
!			ENDIF
!		ELSEIF	(object%Shape==3) THEN		! box
!			IF (MOD(count,2)==0) THEN
!			! Outside Sphere
!				node_region = min(object%Regions(2),node_region)
!!				Phi_p=0
!			ELSE	!Inside object
!				node_region = object%Regions(1)
!				Phi_p=object%Phi
!			ENDIF
!		ENDIF
!
!	ELSEIF (object%Wall(1)%Shape==1) THEN	! Polygonal line
!		
!		count=0	
!		count2=0
!			  	
!!		DO i = 1, n_section
!!			xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!!			IF (i<n_nodes) THEN
!!				xyp(2,1:2)=object%Wall(1)%node(i+1,1:2)
!!			ELSE
!!				xyp(2,1:2)=object%Wall(1)%node(1,1:2)
!!			ENDIF
!!
!!		
!!			IF ((p_node(2)-xyp(1,2))*(p_node(2)-xyp(2,2))<=0) THEN
!!				xi=(p_node(2)-xyp(2,2))*(xyp(2,1)-xyp(1,1))/(xyp(2,2)-xyp(1,2))
!!				xi=xi+xyp(2,1)
!!				IF (xi>=p_node(1) ) THEN
!!					count=count+1
!!					IF (xi==xyp(2,1) .or. xi==xyp(1,1)) THEN
!!						count2=count2+1
!!					ENDIF
!!				ENDIF
!!			ENDIF
!!		ENDDO
!!!		DO i = 1, n_section+2
!!!			
!!!			IF (i<n_section) THEN
!!!			    xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!!!				xyp(2,1:2)=object%Wall(1)%node(i+1,1:2)
!!!!			ELSEIF (i==n_section-1) THEN
!!!!			    xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!!!!				xyp(2,1)=object%Locations(2,1)
!!!!				xyp(2,2)=ytemp
!!!			ELSEIF (i==n_section) THEN
!!!			    xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!!!				xyp(2,1)=object%Locations(2,1)
!!!				xyp(2,2)=ytemp
!!!			ELSEIF (i==n_section+1) THEN
!!!			    xyp(1:2,1)=object%Locations(2,1)
!!!			    xyp(1,2) = ytemp
!!!			    xyp(2,2)=object%Wall(1)%node(n_nodes-1,2)
!!!			ELSEIF (i==n_section+2) THEN
!!!			    xyp(1:2,2)=object%Wall(1)%node(n_nodes,2)
!!!			    xyp(1,1)=object%Locations(2,1)
!!!			    xyp(2,1)=object%Wall(1)%node(n_nodes,1)
!!!!			ELSEIF (i==n_section+3) THEN
!!!!			    xyp(1,1:2)=object%Wall(1)%node(n_nodes,1:2)
!!!!			    xyp(2,1:2)=object%Wall(1)%node(1,1:2)
!!!			ENDIF
!!!
!!!		
!!!			IF ((p_node(2)-xyp(1,2))*(p_node(2)-xyp(2,2))<=0) THEN
!!!				xi=(p_node(2)-xyp(2,2))*(xyp(2,1)-xyp(1,1))/(xyp(2,2)-xyp(1,2))
!!!				xi=xi+xyp(2,1)
!!!				IF (xi>=p_node(1) ) THEN
!!!					count=count+1
!!!					IF (xi==xyp(2,1) .or. xi==xyp(1,1)) THEN
!!!						count2=count2+1
!!!					ENDIF
!!!				ENDIF
!!!			ENDIF
!!!		ENDDO
!!!		count2=count2/2
!!!		count=count-count2
!!!
!!!		IF (object%Wall(1)%Channelwall==3) THEN
!!!			IF ( p_node(2)>object%Locations(2,2)) THEN
!!!				! Outside Box
!!!				node_region = min(object%Regions(2),node_region)
!!!!				Phi_p=0
!!!			ELSE
!!!				IF (MOD(count,2)==0) THEN
!!!				! Outside Box
!!!					node_region = min(object%Regions(2),node_region)
!!!!					Phi_p=0
!!!				ELSE
!!!					node_region = object%Regions(1)
!!!					Phi_p=object%Phi
!!!				ENDIF
!!!			ENDIF
!!!		ELSEIF (object%Wall(1)%Channelwall==4) THEN
!!!			IF ( p_node(2)<object%Locations(1,2)) THEN
!!!				! Outside Box
!!!				node_region = object%Regions(2)
!!!				Phi_p=0
!!!			ELSE
!!!				IF (MOD(count,2)==0) THEN
!!!				! Outside Box
!!!					node_region = min(object%Regions(2),node_region)
!!!!					Phi_p=0
!!!				ELSE
!!!					node_region = object%Regions(1)
!!!					Phi_p=object%Phi
!!!				ENDIF
!!!			ENDIF
!!!		ENDIF
!
!!!!!!!!!!!!!!!!!!!!chj 2014-12-10!!!!!!!!!!!!!!!!!!!!!!!!!
!        IF ( p_node(2)>object%Locations(2,2)) THEN
!            ! Outside Box
!			node_region = min(object%Regions(2),node_region)
!		ELSE
!			xi = p_node(1)
!    	
!	        DO i = 1, n_section 
!	            xyp(1,1:2)=object%Wall(1)%node(i,1:2)
!				xyp(2,1:2)=object%Wall(1)%node(i+1,1:2)
!
!	            IF (xi>=xyp(1,1) .AND. xi<=xyp(2,1)) THEN
!		            ytemp=(xyp(1,2)-xyp(2,2))/(xyp(1,1)-xyp(2,1))
!		            ytemp=ytemp*(xi-xyp(2,1))+xyp(2,2)
!		            EXIT
!	            ENDIF
!	        ENDDO
!    	
!	        IF( p_node(2)> ytemp) THEN
!		        ! Outside Object
!			    node_region = min(object%Regions(2),node_region)
!			ELSE
!			    ! Inside Object
!			    node_region = object%Regions(1)
!				Phi_p=object%Phi
!	        ENDIF
!				
!		ENDIF
!
!!!!!!!!!!!!!!!!!!!!chj 2014-12-10!!!!!!!!!!!!!!!!!!!!!!!!!
!        
!	ENDIF

ENDIF

END