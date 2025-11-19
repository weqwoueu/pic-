SUBROUTINE generate_Adaptive_DG_P_average(HP, HP_average)
!-------------------------LY Add for P average,2021-11-24-------------------------------
USE qsort_c_module_2D
USE ShellSort_2D

IMPLICIT NONE
REAL(8), DIMENSION(:,:), INTENT(IN)   :: HP
REAL(8), DIMENSION(:,:), POINTER  :: HP_average, HP_average_temp, HP_temp

INTEGER :: i, j, k, count1, count2

ALLOCATE(HP_temp(3, SIZE(HP,2)))
HP_temp = 0.0

DO j = 1, SIZE(HP, 2)
  HP_temp(1:2, j) = HP(1:2, j)
  HP_temp(3, j) = j
END DO

!LY REVISE, 2021-12-17, replace QuickSort from ShellSort
!CALL QsortC(HP_temp)
CALL ShellSort(HP_temp, 1)

i = 1
j = 1

DO WHILE (i <= SIZE(HP_temp,2))
  DO WHILE (ABS(HP_temp(1,i)-HP_temp(1,j))<1.0D-12)
    j = j + 1
    IF (j>SIZE(HP_temp,2)) THEN
      EXIT
    END IF
  END DO
  IF (j-1 <= SIZE(HP_temp,2)) THEN
    !CALL QsortC(HP_temp(2:3, i:j-1))
    CALL ShellSort(HP_temp(2:3, i:j-1), 1)
    i = j
  END IF
END DO

ALLOCATE(HP_average_temp(6,SIZE(HP_temp,2)))
HP_average_temp = 0.0

i = 1
j = 1
k = 1
count1 = 1
count2 = 0

DO WHILE (i <= SIZE(HP_temp,2))
  DO WHILE (ABS(HP_temp(1,i)-HP_temp(1,j))<1.0D-12)
    j = j + 1
    IF (j > SIZE(HP_temp, 2)) THEN
      EXIT
    END IF
  END DO

  count2 = 0
  DO WHILE (k <= j-1)
    IF (ABS(HP_temp(2, i)-HP_temp(2, k))<1.0D-12) THEN
      HP_average_temp(1:2, count1) = HP_temp(1:2, i)
      HP_average_temp(3+count2, count1) = HP_temp(3, k)
      count2 = count2 + 1
    ELSE
      count2 = 0
      count1 = count1 + 1
      HP_average_temp(1:2, count1) = HP_temp(1:2, k)
      HP_average_temp(3+count2, count1) = HP_temp(3, k)
      i = k
      count2 = count2 + 1
    END IF
    k = k + 1
  END DO
  count1 = count1 + 1

  i = j
  k = j
END DO

count1 = count1 - 1
ALLOCATE(HP_average(6, count1))
HP_average = 0.0

DO i = 1, 6
  DO j = 1, count1
    HP_average(i, j) = HP_average_temp(i, j)
  END DO
END DO

DEALLOCATE(HP_temp, HP_average_temp)

!P_average:
!P_average(1:2, i): the coordinates of i coordinate node
!P_average(3:6, i): the corresponding mesh node index of i coordinate node

END SUBROUTINE generate_Adaptive_DG_P_average