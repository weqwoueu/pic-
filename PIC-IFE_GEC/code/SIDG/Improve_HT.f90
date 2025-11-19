Subroutine Improve_HT(HT, HP, HE, repeat_refinement)

USE Domain_2D

IMPLICIT NONE
REAL(8), DIMENSION(:,:), POINTER    ::	HP
INTEGER, DIMENSION(:,:), POINTER    ::	HT
INTEGER, DIMENSION(:,:), POINTER    ::  HE
INTEGER                             ::  repeat_refinement

INTEGER, ALLOCATABLE                ::  NodeMap(:,:,:)

INTEGER                             ::  i, j , k, l, count, self, bro
REAL                                ::  dx, dy
INTEGER, ALLOCATABLE                ::  NodeMap_temp(:,:), HT_map(:) , HE_TEMP(:,:)
REAL(8), ALLOCATABLE                ::  HP_temp(:,:)

ALLOCATE(NodeMap((2**repeat_refinement)*(nx - 1) + 1, (2**repeat_refinement)*(ny -1) + 1, 5))
NodeMap = 0

dxmin = minval(HP(1,:))
dymin = minval(HP(2,:))

dx = hx(1) / (2**repeat_refinement)
dy = hx(2) / (2**repeat_refinement)

Do i = 1, SIZE(HP, 2)
    
    j = ((HP(1, i)- dxmin) / dx ) + 1
    k = ((HP(2, i)- dymin) / dy ) + 1
    
    NodeMap(j, k, 1) = NodeMap(j, k, 1) + 1
    l = NodeMap(j, k, 1) + 1
    
    IF (l > 5) THEN
        write(*, *) 'error occur in Improve_HT'
        stop
    END IF
    
    NodeMap(j, k, l) = i 
END Do


Count = 0
DO i = 1, (2**repeat_refinement)*(nx - 1) + 1
    DO j = 1, (2**repeat_refinement)*(ny - 1) + 1
        
        IF (NodeMap(i, j, 1) /= 0) THEN
            Count = Count + 1
        END IF
        
    END DO
END DO


ALLOCATE(NodeMap_temp(5, count))
ALLOCATE(HP_temp(2, count))
NodeMap_temp(:, :) = 0

count = 0
DO i = 1, (2**repeat_refinement)*(nx - 1) + 1
    DO j = 1, (2**repeat_refinement)*(ny - 1) + 1

        IF (NodeMap(i, j, 1) /= 0) THEN
            count = count + 1
            NodeMap_temp(:, count) = NodeMap(i, j, :)
        END IF  
        
    END DO
END DO

DO i = 1, SIZE(NodeMap_temp, 2)
        l = NodeMap_temp(2, i)
        HP_temp(:, i) = HP(:, l)
END DO

ALLOCATE(HT_map(size(HP,2)))
HT_map = 0
DO i = 1, SIZE(NodeMap_temp, 2)
    DO j = 2, NodeMap_temp(1, i) + 1
        l = NodeMap_temp(j, i)
        HT_map(l) = i
    END DO
END DO

DO i = 1, size(HT, 2)
    DO j = 1, 4
        HT(j, i) = HT_map(HT(j, i))
    END DO
END DO

Deallocate(HP)
Allocate(HP(2, size(HP_temp,2)))

DO i = 1, size(HP_temp,2)
    HP(:,i) =  HP_temp(:,i) 
END DO

count = 0
DO i = 1, SIZE(HE, 2)
    self = HE(5, i)
    bro = HE(6, i)
    
    count = count + 1
    HE(1, i) = HT_map(HE(1, i))
    HE(2, i) = HT_map(HE(2, i))
    
    IF (bro /= 0) THEN
        
    
        IF (HT(6, self) /= 0 .And. HT(6, self) /=  HT(6, bro)) THEN
            HE(8, i) = 1
        ELSE
            HE(8, i) = 0
        END IF
    ELSE
        HE(8, i) = 0
    END IF
    
END DO


!Allocate(HE_temp(8, count))
!count = 0
!DO i = 1, SIZE(HE, 2)
!    IF (HE(8, i) == 1) THEN
!        count = count + 1
!        HE_temp(:, count) = HE(:, i)
!    END IF
!END DO
!
!Deallocate(HE)
!ALLOCATE(HE(8, count))
!DO i = 1, size(HE_temp, 2)
!    HE(:, i) = HE_temp(:, i)
!END DO

Deallocate(NodeMap, NodeMap_temp, HP_temp, HT_map)

End Subroutine Improve_HT