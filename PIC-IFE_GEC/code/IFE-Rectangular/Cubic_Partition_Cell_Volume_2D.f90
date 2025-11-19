!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: Cubic_Partition_Cell_Volume_2D.f90                 C
!
!  Purpose: Get the coordinate of volume cell and index
!                                                                      C
!  Author: Yuchuan Chu                                Date: 13-Dec-12  C
!  Comments:                                                           C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
SUBROUTINE Cubic_Partition_Cell_Volume_2D(dimensions, n_nodes, P, T)

IMPLICIT NONE

REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions
INTEGER, DIMENSION(2), INTENT(IN)		::	n_nodes

REAL(8), DIMENSION(:,:), POINTER 	::	P
INTEGER, DIMENSION(:,:), POINTER 	::	T
INTEGER, DIMENSION(:,:), POINTER    ::  Q  

REAL(8)								::	x_s, x_l, y_s, y_l, h_x, h_y
REAL(8), DIMENSION(:), POINTER		::	x, y
INTEGER								::	n_x, n_y, count, i, j, n, row, column, n_grid
INTEGER                                ::  nnx, nny


x_s = dimensions(1,1);
x_l = dimensions(1,2);
y_s = dimensions(2,1);
y_l = dimensions(2,2);
                 
n_x = n_nodes(1);
n_y = n_nodes(2);
nnx = n_x - 1
nny = n_y - 1

n_grid=  (n_x-1)*(n_y-1)
                 
h_x = (x_l - x_s)/(2*(nnx - 1));
h_y = (y_l - y_s)/(2*(nny - 1));

ALLOCATE(x(n_x), y(n_y))

DO count=1,n_x
    IF(count==1 .OR. count==2) THEN
        x(count) = x_s + h_x*(count-1)
    ELSEIF(count==n_x) THEN
        x(count) = x_l
    ELSE
        x(count) = x(2) + 2*h_x*(count-2)
    ENDIF
END DO
DO count=1,n_y
    IF(count==1 .OR. count==2) THEN
        y(count) = y_s + h_y*(count-1)
    ELSEIF(count==n_y) THEN
        y(count) = y_l
    ELSE
        y(count) = y(2) + 2*h_y*(count-2)
    ENDIF
END DO

ALLOCATE(P(2,n_x*n_y))
P (:,:) = 0.

count = 1;

DO i = 1,n_x
     DO j = 1,n_y 
        P(1,count) = x(i)
        P(2,count) = y(j)
        count = count + 1
    END DO
END DO
ALLOCATE(Q(n_x,n_y))
ALLOCATE(T(4, n_grid))
T(:,:)=0
Q(:,:)=0

DO i=1,n_x
   DO j=1,n_y
         Q(i,j)=(i-1)*n_y+j
   END DO
END DO

DO n=1,(n_x-1)*(n_y-1)
   IF (MOD(n,(n_y-1))==0) THEN
        row = n_y-1
		column = n/(n_y-1)
   ELSE
        row = mod(n,n_y-1)
		column = n/(n_y-1)+1
   ENDIF
    T(1,n) = Q(column,row)
    T(2,n) = Q(column+1,row)
    T(3,n) = Q(column+1,row+1)
    T(4,n) = Q(column,row+1)
ENDDO
END
   

              


