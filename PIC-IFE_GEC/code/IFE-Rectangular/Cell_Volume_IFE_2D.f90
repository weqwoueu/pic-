SUBROUTINE Cell_Volume_IFE_2D(vert, information_1, information_2, node_index_el, area)

IMPLICIT NONE

REAL(8), DIMENSION(2,4), INTENT(IN)              :: vert
INTEGER, DIMENSION(18), INTENT(IN)               :: information_1
REAL(8), DIMENSION(8), INTENT(IN)                :: information_2
INTEGER, DIMENSION(4), INTENT(IN)                ::	node_index_el
REAL(8)                                                area

REAL(8)                            xmin, ymin, xmax, ymax
REAL(8)                            Dx, Dy, Ex, Ey
INTEGER                            el_type, n, num, n_temp
REAL(8)                            x_temp, y_temp
INTEGER                            i, j
REAL(8)                            area_rec, area_temp
REAL(8), DIMENSION(2,2)      ::  vert_two
REAL(8)                        ::  height, wide1, wide2

n = 0
num = 0

xmin = MINVAL(vert(1,:))
ymin = MINVAL(vert(2,:))
xmax = MAXVAL(vert(1,:))
ymax = MAXVAL(vert(2,:))

area_rec = (xmax - xmin) * (ymax - ymin)

Dx = information_2(3)
Dy = information_2(4)
Ex = information_2(5)
Ey = information_2(6)

el_type = information_1(6)

IF(el_type == 1) THEN              ! cut two adjacent edges
    n_temp = 0
    DO j=1, 4
        IF(node_index_el(j) == -1) THEN
            n_temp = n_temp + 1
        ENDIF
    ENDDO
    IF(n_temp == 1) THEN
        DO j=1, 4
            IF(node_index_el(j) == -1) THEN
                n = j
            ENDIF
        ENDDO
        x_temp = vert(1,n)
        y_temp = vert(2,n)
        area = ABS(0.5 * ((Ex - Dx) * (y_temp - Dy) - (x_temp - Dx) * (Ey - Dy)))
    ELSEIF(n_temp == 3) THEN
        DO j=1, 4
            IF(node_index_el(j) < -1) THEN
                n = j
            ENDIF
        ENDDO      
        x_temp = vert(1,n)
        y_temp = vert(2,n)
        area_temp = ABS(0.5 * ((Ex - Dx) * (y_temp - Dy) - (x_temp - Dx) * (Ey - Dy)))
        area = area_rec - area_temp
    ENDIF    
ELSEIF(el_type == 2) THEN         ! cut two opposite edges
    DO j=1, 4
        IF(node_index_el(j) == -1) THEN
            num = num+1
            vert_two(:,num) = vert(:,j)
        ENDIF
    ENDDO
    IF(Dx==xmin .OR. Dx==xmax) THEN    ! cut two edges: xmin and xmax
        height = ABS(vert_two(1,1) - vert_two(1,2))
        IF(Dx==vert_two(1,1)) THEN
            wide1 = ABS(Dy - vert_two(2,1))
            wide2 = ABS(Ey - vert_two(2,2))
        ENDIF
        IF(Dx==vert_two(1,2)) THEN
            wide1 = ABS(Dy - vert_two(2,2))
            wide2 = ABS(Ey - vert_two(2,1))
        ENDIF
        area = 0.5*(wide1 + wide2)*height
    ENDIF
    IF(Dy==ymin .OR. Dy==ymax) THEN
        height = ABS(vert_two(2,1)-vert_two(2,2))
        IF(Dy==vert_two(2,1)) THEN
            wide1 = ABS(Dx - vert_two(1,1))
            wide2 = ABS(Ex - vert_two(1,2))
        ENDIF
        IF(Dy==vert_two(2,2)) THEN
            wide1 = ABS(Dx - vert_two(1,2))
            wide2 = ABS(Ex - vert_two(1,1))
        ENDIF
        area = 0.5*(wide1 + wide2)*height
    ENDIF
ELSE
    PRINT*, 'wrong element type, check!'
    stop
ENDIF    

END SUBROUTINE