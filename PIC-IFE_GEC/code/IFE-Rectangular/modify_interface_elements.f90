SUBROUTINE	modify_interface_elements(h, t_c, p_basic, element_index, information_1, information_2,	&
										modified_element_index, modified_information_1, modified_information_2)


IMPLICIT NONE

REAL(8), INTENT(IN)								::	h(2)
INTEGER, DIMENSION(:,:), INTENT(IN)				::	t_c
REAL(8), DIMENSION(:,:), INTENT(IN)				::	p_basic
INTEGER, DIMENSION(:), INTENT(IN)				::	element_index
INTEGER, DIMENSION(:,:), INTENT(IN)				::  information_1
REAL(8), DIMENSION(:,:), INTENT(IN)				::  information_2
INTEGER, DIMENSION(:), POINTER					::	modified_element_index
INTEGER, DIMENSION(:,:), POINTER				::  modified_information_1
REAL(8), DIMENSION(:,:), POINTER				::  modified_information_2


INTEGER											::	num_of_elements
INTEGER											::	n, count, ien
REAL(8)											::	vertices(2,4)
INTEGER											::	i, j
INTEGER											::	interface_element_type
INTEGER											::	pointer_reference_to_local(4), plus_piece_flag
REAL(8)											::	Dx, Dy, Ex, Ey
REAL(8)											::	x1, x2, x3
REAL(8)											::	y1, y2, y3
REAL(8)											::	area, area1, area2

num_of_elements = SIZE(t_c,2)
modified_element_index = element_index
count = 1
!ALLOCATE(modified_information_1(SIZE(information_1,1),SIZE(inforamtion_1,2)))
!ALLOCATE(modified_information_2(SIZE(information_2,1),SIZE(inforamtion_2,2)))

DO n=1, num_of_elements
	ien = element_index(n)
	IF(ien > 0)THEN
		DO i=1,4
			vertices(:,i) = p_basic(:,t_c(i,n))
		ENDDO
		interface_element_type = information_1(6,ien)
		pointer_reference_to_local = information_1(11:14,ien)
		plus_piece_flag = information_1(15,ien)
		Dx = information_2(3,ien)
		Dy = information_2(4,ien)
		Ex = information_2(5,ien)
		Ey = information_2(6,ien)

		IF(interface_element_type == 1)THEN
			x1 = Dx
			y1 = Dy
			x2 = Ex
			y2 = Ey
			x3 = vertices(1, pointer_reference_to_local(1))
			y3 = vertices(2, pointer_reference_to_local(1))
			area = ABS(0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)))
			IF(area < h(1)**2 * h(2)**2)THEN
				IF(plus_piece_flag == 1)THEN
					modified_element_index(n) = -2
				ELSEIF(plus_piece_flag == 2)THEN
					modified_element_index(n) = -1
				ENDIF
			ELSE
				modified_information_1(:,count) = information_1(:,ien)
                modified_information_2(:,count) = information_2(:,ien)
                count=count+1
			ENDIF
		
		ELSEIF(interface_element_type == 2)THEN
			x1 = Dx
			y1 = Dy
			x2 = Ex
			y2 = Ey
			x3 = vertices(1, pointer_reference_to_local(1))
			y3 = vertices(2, pointer_reference_to_local(1))
			area1 = ABS(0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)))
            x1 = Dx
            y1 = Dy
            x2 = vertices(1,pointer_reference_to_local(4))
            y2 = vertices(2,pointer_reference_to_local(4))
            x3 = vertices(1,pointer_reference_to_local(1))
            y3 = vertices(2,pointer_reference_to_local(1))
			area2 = ABS(0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)))
			area = area1 + area2
			IF(area < h(1)**2 * h(2)**2)THEN
				IF(plus_piece_flag == 1)THEN
					modified_element_index(n) = -2
				ELSEIF(plus_piece_flag == 2)THEN
					modified_element_index(n) = -1
				ENDIF
			ELSE
				x1 = Dx
				y1 = Dy
				x2 = Ex
				y2 = Ey
				x3 = vertices(1, pointer_reference_to_local(2))
				y3 = vertices(2, pointer_reference_to_local(2))			
				area1 = ABS(0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)))
				x1 = Dx
				y1 = Dy
				x2 = vertices(1,pointer_reference_to_local(2))
				y2 = vertices(2,pointer_reference_to_local(2))
				x3 = vertices(1,pointer_reference_to_local(3))
				y3 = vertices(2,pointer_reference_to_local(3))
				area2 = ABS(0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)))
				area = area1 + area2
				IF(area < h(1)**2 * h(2)**2)THEN
					IF(plus_piece_flag == 1)THEN
						modified_element_index(n) = -1
					ELSEIF(plus_piece_flag == 2)THEN
						modified_element_index(n) = -2
					ENDIF
				ELSE
					modified_information_1(:,count) = information_1(:,ien)
					modified_information_2(:,count) = information_2(:,ien)
					count=count+1
				ENDIF

			ENDIF

		ENDIF

	ENDIF

ENDDO

count = 1

DO j=1, num_of_elements
	IF(modified_element_index(j) > 0)THEN
		modified_element_index(j) = count
		count = count +1
	ENDIF
ENDDO
	

END SUBROUTINE


					


				
				
