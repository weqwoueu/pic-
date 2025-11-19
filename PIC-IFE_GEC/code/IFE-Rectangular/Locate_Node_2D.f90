SUBROUTINE Locate_Node_2D	(	p_node, object, node_region		)

! Purpose:		Locate node with respect to an immersed element
! Last Update:	4/21/2004 07:09 PM

USE IFE_MAIN_PARAM
USE Object_Data_2D

IMPLICIT NONE

REAL(8), DIMENSION(2), INTENT(IN)	::	p_node
TYPE(ObjectType), INTENT(IN)		::	object
INTEGER									node_region, coord1, coord2


REAL(8), DIMENSION(2)				::	X1, X2, Xb, Xt, Xp, Xc, Xtb, Xpb, Xpt, Xpc
REAL(8)									r, XpbDotXtb, XtbDotXtb, XptDotXtb,	XpcDotXpc, muc

!INTEGER	coord1, coord2, coord3

REAL(8) Norm_2 

IF		(object%Axis==1) THEN
		coord1 = 1
		coord2 = 2
ELSEIF	(object%Axis==2) THEN
		coord1 = 2
		coord2 = 1
ENDIF


IF (object%Shape==1 .AND. object%Axis==0) THEN	! Circle

!print*, object%Regions(:)
!print*, p_node(1), p_node(2)
!print*, object%Locations(1,:)
!print*, p_node-object%Locations(1,:)
!print*

	IF (	Norm_2( p_node-object%Locations(1,:))<=(object%Dimensions(1)+SmallValue) ) THEN
		! Inside Sphere
		node_region = object%Regions(1)
	ELSE
		! Outside Sphere
	    node_region = object%Regions(2)
	ENDIF	

ELSEIF (object%Shape==3) THEN	! Box
	
	IF (	p_node(1)>(object%Locations(1,1)-SmallValue).AND. p_node(1)<(object%Locations(2,1)+SmallValue).AND.	& 
			p_node(2)>(object%Locations(1,2)-SmallValue).AND. p_node(2)<(object%Locations(2,2)+SmallValue)) THEN
		! Inside 
		node_region = object%Regions(1)
	ELSE
		! Outside 
		node_region = object%Regions(2)
	ENDIF

ELSEIF (object%Shape==11) THEN	! Line, Plate parallel to X axis(delta ==0) or Z axis(delta ==1)

		X1 = object%Locations(1,:)
		X2 = object%Locations(2,:)

!	IF (	(X2(coord1)-p_node(coord1))*(p_node(coord1)-X1(coord1)) > -SmallValue .AND.	&
!			ABS(X2(coord2)-p_node(coord2)) < SmallValue) THEN
	IF (	(X2(1)-p_node(1))*(p_node(1)-X1(1)) >= -SmallValue						.AND.	&
			(X2(2)-p_node(2))*(p_node(2)-X1(2)) >= -SmallValue						.AND.	&
			DABS((p_node(coord2)-X1(coord2))*(X2(coord1)-X1(coord1)))<= SmallValue	) THEN

		! Inside 
		node_region = object%Regions(1)
	ELSE
		! Outside 
		node_region = object%Regions(2)
	ENDIF
	
ENDIF



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
!	IF (	p_node(coord1)>object%Locations(1,coord1) .AND.	& 
!			p_node(coord1)<object%Locations(2,coord1) .AND.	&
!			Norm_3(p_node((/coord2,coord3/))-object%Locations(1,(/coord2,coord3/)))<object%Dimensions(1) ) THEN
!		! Inside Cylinder
!		node_region = object%Regions(1)
!	ELSE
!		! Outside Cylinder
!		node_region = object%Regions(2)
!	ENDIF
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
!	Xp  = p_node
!	Xpb = Xp - Xb
!	CALL Vec_Dot_3(Xpb, Xtb, XpbDotXtb)
!	muc = -XpbDotXtb/XtbDotXtb
!	Xc  = Xb + muc*(Xt - Xb)
!	Xpt = Xp - Xt
!	CALL Vec_Dot_3(Xpt, Xtb, XptDotXtb)
!	Xpc = Xp - Xc
!	CALL Vec_Dot_3(Xpc, Xpc, XpcDotXpc)
!
!	IF (XpbDotXtb>=-MZero .AND. XptDotXtb<=MZero .AND. XpcDotXpc<=r**2) THEN
!	! Inside the object
!		node_region = object%Regions(1)
!	ELSE
!		! Outside Cylinder
!		node_region = object%Regions(2)
!	ENDIF
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
!
!	IF (	p_node(coord1)>object%Locations(2,1)							.AND.	& 
!			p_node(coord1)<object%Locations(1,coord1)+object%Dimensions(1)	.AND.	&
!			Norm_3( p_node-object%Locations(1,:))<object%Dimensions(1) ) THEN
!		! Inside Sphere
!		node_region = object%Regions(1)
!	ELSE
!		! Outside Sphere
!		node_region = object%Regions(2)
!	ENDIF
!
!ELSEIF (object%Shape==2 .AND. object%Axis==0) THEN	! Sphere
!
!	IF (	Norm_3( p_node-object%Locations(1,:))<object%Dimensions(1) ) THEN
!		! Inside Sphere
!		node_region = object%Regions(1)
!	ELSE
!		! Outside Sphere
!		node_region = object%Regions(2)
!	ENDIF	
!	
!ELSEIF (object%Shape==3) THEN	! Box
!
!	
!	IF (	p_node(1)>object%Locations(1,1)	.AND. p_node(1)<object%Locations(2,1)	.AND.	& 
!			p_node(2)>object%Locations(1,2)	.AND. p_node(2)<object%Locations(2,2)	.AND.	&
!			p_node(3)>object%Locations(1,3)	.AND. p_node(3)<object%Locations(2,3)	) THEN
!		! Inside Sphere
!		node_region = object%Regions(1)
!	ELSE
!		! Outside Sphere
!		node_region = object%Regions(2)
!	ENDIF
!
!ENDIF


END