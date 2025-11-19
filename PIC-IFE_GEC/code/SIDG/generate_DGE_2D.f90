SUBROUTINE generate_DGE_2D(DGP, DGT, DGE, dimensions)
!----------------2021/05/06 wsy add-----------------------------------------------
!生成DG-FE中的边界边
!第一列和第二列表示边界边起点和终点的全局坐标
!第三列和第四列表示边界边起点和终点的局部编号
!第五列和第二六表示边界边所在网格元和邻边（如果是外边界，邻边默认为0）
!这里还没有进行网格加密
IMPLICIT NONE

REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions
REAL(8), DIMENSION(:,:), POINTER 	::	DGP
INTEGER, DIMENSION(:,:), POINTER	::	DGT, DGE
INTEGER ::  number_of_edges, number_of_element, number_of_row_element, n
REAL(8) ::  xmax, xmin, ymax, ymin, hx, hy

number_of_edges = SIZE(DGP, 2)
ALLOCATE(DGE(7, number_of_edges))
number_of_element = SIZE(DGT, 2)
xmin = dimensions(1, 1)
ymin = dimensions(2, 1)
xmax = dimensions(1, 2)
ymax = dimensions(2, 2)
hx = DGP(1, DGT(2, 1)) - DGP(1, DGT(1, 1))
hy = DGP(2, DGT(4, 1)) - DGP(2, DGT(1, 1))
number_of_row_element = (ymax - ymin) /  hy

DO n = 1, number_of_edges
    
    IF (MOD(n,4) == 0) THEN
        DGE(1,n) = DGT(4, INT(n/4))
        DGE(2,n) = DGT(1, INT(n/4))
        DGE(3,n) = 4
        DGE(4,n) = 1
        DGE(5,n) = INT(n/4)
        IF (DGE(5,n) <= number_of_row_element) THEN
            DGE(6,n) = 0
        ELSE
            DGE(6,n) = (DGE(5,n) - number_of_row_element) * 4 - 2
        END IF
     ELSE
        DGE(1,n) = DGT(MOD(n,4),INT(n/4) + 1)
        DGE(2,n) = DGT(MOD(n,4) + 1,INT(n/4) + 1)
        DGE(3,n) = MOD(n,4)
        DGE(4,n) = MOD(n,4) + 1
        DGE(5,n) = INT(n/4) + 1
        IF (MOD(n,4) == 1) THEN
            IF (MOD(DGE(5,n),number_of_row_element) == 1) THEN
                DGE(6,n) = 0
            ELSE
                DGE(6,n) = (DGE(5,n) - 1)*4 - 1
            END IF
        ELSEIF (MOD(n,4) == 2) THEN
            IF (DGE(5,n) > number_of_element - number_of_row_element) THEN
                DGE(6,n) = 0
            ELSE
                DGE(6,n) = (DGE(5,n) + number_of_row_element)*4
            END IF
        ELSEIF (MOD(n,4) == 3) THEN
            IF (MOD(DGE(5,n),number_of_row_element) == 0) THEN
                DGE(6,n) = 0
            ELSE
                DGE(6,n) = (DGE(5,n) + 1)*4 - 3
            END IF
        END IF
    
     END IF
     
     ! 因为局部网格加密的需要，要放上坐标信息(第七列后面会用上)
     DGE(7,n) = 0
     
END DO


END