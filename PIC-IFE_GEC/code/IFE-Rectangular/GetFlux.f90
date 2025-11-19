SUBROUTINE	GetFlux(r, x, y)

USE Object_2D
USE Domain_2D
USE IFE_Data
IMPLICIT NONE

REAL(8), INTENT(IN)				::	x, y
REAL(8)							::	r

REAL(8)                            ::  hx_temp, hy_temp
REAL(8)                            ::  xmin, ymin
INTEGER                            ::  nx_temp, ny_temp, ipre, jpre
INTEGER                            ::  num_of_element, num_of_interface_element, intr_num_inter_ele
INTEGER                            ::  i,j

REAL(8)                            ::  Dx, Dy, Ex, Ey
REAL(8)                            ::  length, length_D, length_E


hx_temp = hx(1)
hy_temp = hx(2)

OPEN(10, ACTION='READ', FILE='normalize.inp')!
READ(10,*) xmin, ymin
READ(10,*) nx_temp, ny_temp
CLOSE(10)

jpre = INT((y - ymin) / hy_temp) + 1
ipre = INT((x - xmin) / hx_temp) + 1
num_of_element = (ipre-1) * (ny_temp-1) + jpre
num_of_interface_element = SIZE(information_2,2)


DO i=1, num_of_interface_element
    IF(information_1(5,i) == num_of_element) THEN
        intr_num_inter_ele = i
        EXIT
    ENDIF
ENDDO

Dx = information_2(3,intr_num_inter_ele)
Dy = information_2(4,intr_num_inter_ele)
Ex = information_2(5,intr_num_inter_ele)
Ey = information_2(6,intr_num_inter_ele)

length = SQRT((Dx-Ex)**2+(Dy-Ey)**2)
length_D = SQRT((x-Dx)**2+(y-Dy)**2)
length_E = SQRT((x-Ex)**2+(y-Ey)**2)

r = collectq(5,intr_num_inter_ele) + (collectq(6,intr_num_inter_ele)    &
     - collectq(5,intr_num_inter_ele)) * (length_D / length)



END SUBROUTINE
