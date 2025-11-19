SUBROUTINE	Line_Rectangular_Intersection_BL(line_ends, Rectangular_Size, Corner_Points, inters_point, inters_flag)


USE	IFE_MAIN_PARAM
IMPLICIT NONE

REAL(8), INTENT(IN)		::	line_ends(2,2), Rectangular_Size(2), Corner_Points(2,2)
REAL(8), INTENT(OUT)	::	inters_point(2)
INTEGER, INTENT(OUT)	::	inters_flag

REAL(8)					::	x1, y1, x2, y2
REAL(8)					::	x_min, y_min, x_max, y_max
REAL(8)					::	object_x_min, object_x_max, object_y_min, object_y_max
INTEGER					::	n_root
REAL(8)					::	delta1, delta2, delta_min, delta_max

inters_point = (/0.,0./)
inters_flag = 0

x1 = line_ends(1,1)
y1 = line_ends(1,2)
x2 = line_ends(2,1)
y2 = line_ends(2,2)

object_x_min = Corner_Points(1,1)
object_x_max = Corner_Points(2,1)
object_y_min = Corner_Points(1,2)
object_y_max = Corner_Points(2,2)


IF(ABS(y1-y2) < SmallValue) THEN                           !$ horizontal edge
	x_min = MIN(x1, x2)                     !$ x minimum of the edge
	x_max = MAX(x1, x2)                     !$ x maximum of the edge
	delta1 = y1 - object_y_min              !$ the distance of edge first point from bottom of object
	delta2 = y1 - object_y_max
	IF(delta1 < -SmallValue .OR. delta2 > SmallValue) THEN
		n_root = 0
	ELSE
		delta_min = x_min - object_x_min    
		delta_max = x_max - object_x_max
		IF(delta_min < -SmallValue) THEN
			n_root = 1
			inters_point(1) = object_x_min
			inters_point(2) = y1
			inters_flag = 1
		ELSEIF(delta_max > SmallValue) THEN
			n_root = 1
			inters_point(1) = object_x_max
			inters_point(2) = y1
			inters_flag = 1
		ELSE
			n_root = 0
		ENDIF
	ENDIF
ENDIF

IF(ABS(x1-x2) < SmallValue) THEN                           !$ vertical edge
	y_min = MIN(y1, y2)
	y_max = MAX(y1, y2)
	delta1 = x1 - object_x_min
	delta2 = x1 - object_x_max
	IF(delta1 < -SmallValue .OR. delta2 > SmallValue) THEN
		n_root = 0
	ELSE
		delta_min = y_min - object_y_min
		delta_max = y_max - object_y_max
		IF(delta_min < -SmallValue) THEN
			n_root = 1
			inters_point(1) = x1
			inters_point(2) = object_y_min
			inters_flag = 1
		ELSEIF(delta_max > SmallValue) THEN
			n_root = 1
			inters_point(1) = x1
			inters_point(2) = object_y_max
			inters_flag = 1
		ELSE
			n_root = 0
		ENDIF
	ENDIF
ENDIF			



END SUBROUTINE
