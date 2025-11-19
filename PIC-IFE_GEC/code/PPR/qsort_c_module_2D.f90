!2D Quick Sort
!LY Add for P_average, 2021-11-24
MODULE qsort_c_module_2D
  IMPLICIT NONE
  PUBLIC :: QsortC
  PRIVATE :: Partition
  
  CONTAINS
    
    RECURSIVE SUBROUTINE QsortC(A)
      REAL(8), INTENT(IN OUT), DIMENSION(:,:) :: A
      INTEGER :: iq
      IF (SIZE(A,2) > 1) THEN
        CALL Partition(A, iq)
        CALL QsortC(A(:, :iq-1))
        CALL QsortC(A(:, iq:))
      END IF
    END SUBROUTINE QsortC
    
    SUBROUTINE Partition(A, marker)
      REAL(8), INTENT(IN OUT), DIMENSION(:,:) :: A
      INTEGER, INTENT(OUT) :: marker
      INTEGER :: i, j
      REAL(8),DIMENSION(SIZE(A,1)) :: temp
      REAL(8) :: x     !pivot point
      x = A(1,1)
      i = 1
      j = SIZE(A,2)
      temp = 0.0
      DO
        DO WHILE (A(1,j) >= x .And. j > i)
            j = j - 1
        END DO
        DO WHILE (A(1,i) <= x .And. j > i)
            i = i + 1
        END DO
        IF (i < j) THEN
          ! exchange A(i) and A(j)
          temp = A(:,i)
          A(:,i) = A(:,j)
          A(:,j) = temp
        ELSEIF (i == j) THEN
          marker = i + 1
          temp = A(:,i)
          A(:,i) = A(:,1)
          A(:,1) = temp
          RETURN
        ELSE
          temp = A(:,i)
          A(:,i) = A(:,1)
          A(:,1) = temp
          marker = i
          RETURN
        END IF
      END DO
    END SUBROUTINE Partition
  
END MODULE qsort_c_module_2D