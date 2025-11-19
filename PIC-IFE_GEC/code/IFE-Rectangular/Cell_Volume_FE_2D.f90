SUBROUTINE Cell_Volume_FE_2D(vert, area)

IMPLICIT NONE

REAL(8), DIMENSION(2,4), INTENT(IN)              ::  vert
REAL(8)                                            ::  area
REAL(8)                                            ::  xmin, xmax, ymin, ymax

xmin = MINVAL(vert(1,:))
xmax = MAXVAL(vert(1,:))
ymin = MINVAL(vert(2,:))
ymax = MAXVAL(vert(2,:))

area = (xmax - xmin) * (ymax - ymin)

END SUBROUTINE