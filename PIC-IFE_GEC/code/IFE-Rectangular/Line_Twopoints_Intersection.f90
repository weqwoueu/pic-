SUBROUTINE Line_Twopoints_Intersection(line_ends,	xyp,	inters_point, inters_flag)

! Purpose:		Intersect a Box but only the lines parallel to x-axis.
! Last Update:	10/3/2003 07:05 AM
									
USE PIC_MAIN_PARAM_2D
USE Domain_2D

IMPLICIT NONE



!INTEGER						nnxp
REAL(8), INTENT(IN)		::	line_ends(2,2), &
							xyp(2,2)
REAL(8), INTENT(OUT)	::	inters_point(2)
INTEGER, INTENT(OUT)	::	inters_flag
!REAL(8)						LowBound(2),UpBound(2)

REAL(8)	x1,y1, x2,y2, x21,y21, xi,yi,y1i,y2i,k1,k2

x1	= line_ends(1,1)
y1	= line_ends(1,2)

x2	= line_ends(2,1)
y2	= line_ends(2,2)

x21 = x2-x1
y21 = y2-y1

inters_flag 	= 0
inters_point 	= (/Zero,Zero/)



IF (x1 == x2) THEN
	xi = x1
	k1 = (xyp(2,2) - xyp(1,2)) / (xyp(2,1) - xyp(1,1))
	yi = xyp(1,2) + k1 * (xi - xyp(1,1))
	IF ((y1-yi)*(y2-yi)<=0 ) THEN
		IF ((xyp(1,2)-yi)*(xyp(2,2)-yi)<=0 .and. (xyp(1,1)-xi)*(xyp(2,1)-xi)<=0  ) THEN		! Intersection point between Box ends
			IF (y1==yi ) THEN	! ObjLowBound(1)-xi = 0 : Intersection point at first END IF			
				inters_flag		= -1
			ELSEIF (y2==yi ) THEN	! Intersection point at second END IF
				inters_flag		= -2
			ELSE
				inters_flag		= 1
			ENDIF
			inters_point	= (/xi,yi/)
		ENDIF
	ENDIF
ENDIF


IF (x1/=x2 .AND. y2==y1) THEN
	yi=y1
	k2=(xyp(2,1)-xyp(1,1))/(xyp(2,2)-xyp(1,2))
	xi=(y1-xyp(1,2))*k2+xyp(1,1)
	IF ((x1-xi)*(x2-xi)<=0 ) THEN
		IF ((xyp(1,1)-xi)*(xyp(2,1)-xi)<=0 .and. (xyp(1,2)-yi)*(xyp(2,2)-yi)<=0 ) THEN		! Intersection point between Box ends
			IF (x1==xi ) THEN	! ObjLowBound(1)-xi = 0 : Intersection point at first END IF			
				inters_flag		= -1
			ELSEIF (x2==xi ) THEN	! Intersection point at second END IF
				inters_flag		= -2
			ELSE
				inters_flag		= 1
			ENDIF
			inters_point	= (/xi,yi/)
		ENDIF
	ENDIF
ENDIF




END