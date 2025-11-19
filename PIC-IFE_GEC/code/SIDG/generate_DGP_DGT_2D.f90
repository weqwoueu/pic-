SUBROUTINE generate_DGP_DGT_2D(dimensions, n_nodes, P, T, DGP, DGT)

!----------------------- wsy add DGP and DGT 2021/7/27 -------------------------------------
    
IMPLICIT NONE

REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions
INTEGER, DIMENSION(2), INTENT(IN)		::	n_nodes

REAL(8), DIMENSION(:,:), POINTER 	::	DGP, P
INTEGER, DIMENSION(:,:), POINTER	::	DGT, T

INTEGER								::	n_x, n_y, count, i, n_grid, DG_n_nodes

                 
n_x = n_nodes(1);
n_y = n_nodes(2);

n_grid=  (n_x-1)*(n_y-1);
DG_n_nodes = n_grid * 4;

ALLOCATE(DGT(4,n_grid))
DGT (:,:) = 0

DO count = 1, n_grid
    DGT(1, count) =  4 * count - 3;
    DGT(2, count) =  4 * count - 2;
    DGT(3, count) =  4 * count - 1;
    DGT(4, count) =  4 * count;
END DO

ALLOCATE(DGP(2,4 * n_grid))
DGP (:,:) = 0.

DO count = 1, n_grid
    DO i = 1, 4
        
        DGP(1:2, DGT(i,count)) = P(1:2,T(i,count));
        
    END DO
END DO


END