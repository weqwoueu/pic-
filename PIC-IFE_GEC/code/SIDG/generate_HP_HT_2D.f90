SUBROUTINE generate_HP_HT_2D(P, T, DGT_IFE_partition, HP, HT, P_average)
! wsy add SIDG 2021/10/27
USE IFE_MAIN_PARAM

IMPLICIT NONE

REAL(8), DIMENSION(:,:), POINTER 	      ::	P, HP, HP_temp, P_average, P_average_temp
INTEGER, DIMENSION(:,:), POINTER	      ::	T, HT
INTEGER, DIMENSION(:,:), POINTER        ::  DGT, DGT_IFE_partition

LOGICAL                                 ::  flag
INTEGER                                 ::  max_node, i, j, k, count, num
INTEGER, DIMENSION(:), ALLOCATABLE      ::  hashset, hashmap
REAL(8)                                 ::  delta

flag = .false.
max_node = MAXVAL(T)
ALLOCATE(HT(5, SIZE(T,2)))
ALLOCATE(HP_temp(2, SIZE(P,2) * 8)) ! enough big
ALLOCATE(P_average_temp(7, SIZE(P,2)))   !LY Add for Average,2021-11-22

P_average_temp = 0

DO i = 1, SIZE(P,2)
  HP_temp(:, i) = P(:, i)
  P_average_temp(1:2, i) = P(:, i)
END DO

DO i = 1, SIZE(T, 2)
    
    IF (DGT_IFE_partition(5, i) > 0) THEN
        
        max_node = max_node + 1
        HP_temp(:, max_node) = P(:, T(1, i))
        HT(1, i) = max_node
        P_average_temp(3, T(1, i)) = P_average_temp(3, T(1, i)) + 1.0
        P_average_temp(INT((3.0+P_average_temp(3, T(1, i)))), T(1, i)) = max_node
        
        max_node = max_node + 1
        HP_temp(:, max_node) = P(:, T(2, i))
        HT(2, i) = max_node
        P_average_temp(3, T(2, i)) = P_average_temp(3, T(2, i)) + 1.0
        P_average_temp(INT((3.0+P_average_temp(3, T(2, i)))), T(2, i)) = max_node
        
        max_node = max_node + 1
        HP_temp(:, max_node) = P(:, T(3, i))
        HT(3, i) = max_node
        P_average_temp(3, T(3, i)) = P_average_temp(3, T(3, i)) + 1.0
        P_average_temp(INT((3.0+P_average_temp(3, T(3, i)))), T(3, i)) = max_node
        
        max_node = max_node + 1
        HP_temp(:, max_node) = P(:, T(4, i))
        HT(4, i) = max_node
        P_average_temp(3, T(4, i)) = P_average_temp(3, T(4, i)) + 1.0
        P_average_temp(INT((3.0+P_average_temp(3, T(4, i)))), T(4, i)) = max_node
        
        HT(5, i) = 1
        
    ELSE
        
        HT(1:4, i) = T(:, i)
        HT(5, i) = 0
        
    END IF
    
END DO

ALLOCATE(hashset(max_node))
hashset = 0

DO i = 1, SIZE(HT, 2)
    DO j = 1, 4
        
        hashset(HT(j, i)) = hashset(HT(j, i)) + 1
        
    END DO
END DO

ALLOCATE(hashmap(max_node))
hashmap = 0
count = 0

DO i = 1, max_node
    
    IF (hashset(i) > 0) THEN
        count = count + 1
        hashmap(i) = count
    ELSE
        flag = .TRUE.
    END IF
    
END DO

!LY Add for Average,2021-11-22
DO j = 1, SIZE(P_average_temp, 2)

  delta = P_average_temp(3, j) - 4.0
  IF (delta < -SmallValue) THEN
    P_average_temp(3, j) = P_average_temp(3, j) + 1.0
    P_average_temp(INT((3.0+P_average_temp(3, j))), j) = j
  END IF
  
END DO

DO i = 1, max_node
  DO j = 1, SIZE(P_average_temp, 2) 
  
    IF (i==P_average_temp(4,j)) THEN
      P_average_temp(4,j) = hashmap(i)
    ELSEIF (i==P_average_temp(5,j)) THEN
      P_average_temp(5,j) = hashmap(i)
    ELSEIF (i==P_average_temp(6,j)) THEN
      P_average_temp(6,j) = hashmap(i)
    ELSEIF (i==P_average_temp(7,j)) THEN
      P_average_temp(7,j) = hashmap(i)
    END IF
    
  END DO
END DO

DO j = 1, SIZE(P_average_temp, 2)
  num = 0
  DO i = 4, 7
  
    IF (P_average_temp(i, j) /= 0) THEN
      num = num + 1  
    ENDIF
    P_average_temp(3, j) = num
    
  END DO
END DO

ALLOCATE(P_average(6, SIZE(P_average_temp, 2)))
DO i = 1, SIZE(P_average_temp, 2)
  P_average(1:2, i) = P_average_temp(1:2, i)
  P_average(3:6, i) = P_average_temp(4:7, i)
END DO
!LY Add for Average,2021-11-22

IF (flag .EQV. .true.) THEN
    
    DO i = 1, SIZE(HT, 2)
        DO j = 1, 4
            HT(j, i) =  hashmap(HT(j, i))
        END DO
    END DO
    
    ALLOCATE(HP(2, count))
    
    DO i = 1, max_node
        IF (hashmap(i) > 0) THEN
            HP(:, hashmap(i)) = HP_temp(:, i)
        END IF
    END DO
ELSE

    ALLOCATE(HP(2, max_node))
    DO i = 1, max_node
        HP(:, i) = HP_temp(:, i)
    END DO
    
END IF

DEALLOCATE(HP_temp, P_average_temp)

END