SUBROUTINE	Line_Circle_Intersection_BL(line_ends, Circle_radius, Circle_center, inters_point, inters_flag)


USE	IFE_MAIN_PARAM
IMPLICIT NONE

REAL(8), INTENT(IN)		::	line_ends(2,2), Circle_radius, Circle_center(2)
REAL(8), INTENT(OUT)	::	inters_point(2)
INTEGER, INTENT(OUT)	::	inters_flag


REAL(8), DIMENSION(2)	::	root1, root2
REAL(8)					::	r, xc,yc, x1,y1, x2,y2
REAL(8)					::	rx1,ry1, rx2,ry2
INTEGER					::	root1_flag, root2_flag
REAL(8)					::	delta
INTEGER					::	n_root

inters_point = (/0.,0./)
inters_flag = 0

r = Circle_radius

xc = Circle_center(1)
yc = Circle_center(2)

x1 = line_ends(1,1)
y1 = line_ends(1,2)

x2 = line_ends(2,1)
y2 = line_ends(2,2)

root1 = 0.
root2 = 0.

root1_flag = 0
root2_flag = 0

IF(ABS(x1-x2) < SmallValue) THEN
	delta = r*r-(x1-xc)*(x1-xc)
	IF(delta < -SmallValue) THEN
		n_root = 0
	ELSEIF(delta > SmallValue) THEN
		n_root = 2
		root1(1) = x1
		root1(2) = yc+SQRT(delta)
		root2(1) = x1
		root2(2) = yc-SQRT(delta)
	ELSEIF(ABS(delta) < SmallValue) THEN
		n_root = 1
		root1(1) = x1
		root1(2) = yc
		root2(1) = x1
		root2(2) = yc
	ENDIF
ENDIF

IF(x1/=x2 .AND. ABS(y1-y2) < SmallValue) THEN
	delta = r*r-(y1-yc)*(y1-yc)
	IF(delta < -SmallValue) THEN
		n_root = 0
	ELSEIF(delta > SmallValue) THEN
		n_root = 2
		root1(1) = xc+SQRT(delta)
		root1(2) = y1
		root2(1) = xc-SQRT(delta)
		root2(2) = y1
	ELSEIF(ABS(delta) < SmallValue) THEN
		n_root = 1
		root1(1) = xc
		root1(2) = y1
		root2(1) = xc
		root2(2) = y1
	ENDIF
ENDIF

IF(n_root/=0) THEN
	rx1 = x1-root1(1)
	ry1 = y1-root1(2)
	rx2 = x2-root1(1)
	ry2 = y2-root1(2)
	r = rx1*rx2+ry1*ry2

	IF(r < -SmallValue) THEN
		root1_flag = 1
	ELSEIF(ABS(r) < SmallValue) THEN
		IF((rx1*rx1+ry1*ry1) > (rx2*rx2+ry2*ry2)) THEN
			root1_flag = -2
		ELSEIF((rx1*rx1+ry1*ry1) < (rx2*rx2+ry2*ry2))THEN
			root1_flag = -1
		ELSE
			PRINT*,'STOP'
			STOP
		ENDIF
!		IF(ABS(rx1)<SmallValue .AND. ABS(ry1)<SmallValue) THEN
!			root1_flag = -1
!		ENDIF
!		IF(ABS(rx2)<SmallValue .AND. ABS(ry2)<SmallValue) THEN
!			root1_flag = -2
!		ENDIF
!		PRINT*,'STOP, Maybe, it is a bug'
!		STOP
	ELSEIF(r>MZero) THEN
		root1_flag = 0
	ENDIF

	rx1=x1-root2(1)
	ry1=y1-root2(2)
	rx2=x2-root2(1)
	ry2=y2-root2(2)
	r=rx1*rx2+ry1*ry2

	IF(r<-SmallValue) THEN
		root2_flag = 1
	ELSEIF(ABS(r) < SmallValue) THEN
		IF((rx1*rx1+ry1*ry1) > (rx2*rx2+ry2*ry2)) THEN
			root2_flag = -2
		ELSEIF((rx1*rx1+ry1*ry1) < (rx2*rx2+ry2*ry2))THEN
			root2_flag = -1
		ELSE
			PRINT*,'STOP'
			STOP
		ENDIF
!		IF(ABS(rx1)<SmallValue .AND. ABS(ry1)<SmallValue) THEN
!			root2_flag = -1
!		ENDIF
!		IF(ABS(rx2)<SmallValue .AND. ABS(ry2)<SmallValue) THEN
!			root2_flag = -2
!		ENDIF
!		PRINT*,'STOP, Maybe, it is a bug'
!		STOP
	ELSEIF(r>SmallValue) THEN
		root2_flag = 0
	ENDIF
ELSEIF(n_root==0) THEN
	root1_flag = 0
	root2_flag = 0
ENDIF


IF(root1_flag/=0 .and. root2_flag==0) THEN
            inters_flag 	= root1_flag
            inters_point 	= root1
ELSEIF(root1_flag==-1 .and. root2_flag==-1) THEN
            inters_flag 	= root1_flag
            inters_point 	= root1
ELSEIF(root1_flag==-2 .and. root2_flag==-2) THEN
            inters_flag 	= root1_flag
            inters_point 	= root1
ELSEIF(root1_flag==0 .and. root2_flag/=0) THEN
            inters_flag 	= root2_flag
            inters_point 	= root2
ELSEIF((root1_flag==-1 .and. root2_flag==-2).OR.(root1_flag==-2 .and. root2_flag==-1).OR.(root1_flag==0 .and. root2_flag==0)) THEN
            inters_flag 	= 0
ELSE
      PRINT*,'n_root=',n_root
      PRINT*,'root1=',root1, ' root1_flag=',root1_flag 
      PRINT*,'root2=',root2, ' root2_flag=',root2_flag 
	  PRINT*,'LineEnd1=',x1,y1
	  PRINT*,'LineEnd2=',x2,y2
	  PRINT*,'Xc=',xc,' Yc=',yc,' R=',rx1,ry1,rx2,ry2
	  PRINT*,'Check the Mesh, STOP'
	  STOP
ENDIF


END SUBROUTINE
