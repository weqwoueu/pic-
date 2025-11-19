SUBROUTINE generate_E_Dirichlet_2D(HP, HT, E_Dirichlet)
USE IFE_Boundary
!=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========
Use IFE_Data, Only: E_DirichletValue
!=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========
IMPLICIT NONE

REAL(8), DIMENSION(:,:), POINTER 	    ::	HP
INTEGER, DIMENSION(:,:), POINTER	    ::	HT, E_Dirichlet, E_Dirichlet_temp

REAL(8)                                 ::  top, bottom, left, right, node1(2), node2(2)
INTEGER                                 ::  count, i, j, N_boundary, k
Integer                                 ::  IN1, IN2

!=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========
Real(8), Dimension(:), Allocatable :: E_DirichletValueTemp
!=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========

left = MINVAL(HP(1, :))
right = MAXVAL(HP(1, :))
top = MAXVAL(HP(2, :))
bottom = MINVAL(HP(2, :))
count = 0
N_boundary = size(bc_index)

ALLOCATE(E_Dirichlet_temp(7, size(HP,2)))
!LY Add for Boundary Condition, 2021-11-18
!E_Dirichlet_temp(6, :)��1--Dirichlet, 2--Neumann, 3--Robin.

!=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========
Allocate(E_DirichletValueTemp(Size(HP,2)))
!=========LY modification for dealing with Floating Dirichlet boundary condition, 2022-7-25=========

DO i = 1, size(HT, 2)
  IF (HT(5, i) > 0) THEN
    DO j = 1, 4
      IF (j /= 4) THEN
        !IF ( (HP(1, HT(j, i)) == left .AND. HP(1, HT(j + 1, i)) == left) .OR. (HP(1, HT(j, i)) == right .AND. &
        !  HP(1, HT(j + 1, i)) == right) .OR. ((abs(HP(2, HT(j, i)) - top) < 10D-6) .AND. &
        !  (abs(HP(2, HT(j + 1, i)) - top) <  10D-6 ))  .OR. &
        !  (HP(2, HT(j, i)) == bottom .AND. HP(2, HT(j + 1, i)) == bottom) ) THEN
         IF ( (ABS(HP(1, HT(j, i))-left)<1.0D-12 .AND. ABS(HP(1, HT(j + 1, i))-left)<1.0D-12) .OR. &
              (ABS(HP(1, HT(j, i))-right)<1.0D-12 .AND. ABS(HP(1, HT(j + 1, i))-right)<1.0D-12) .OR. &
              (ABS(HP(2, HT(j, i))-top)<1.0D-12 .AND. ABS(HP(2, HT(j + 1, i))-top)<1.0D-12) .OR. &
              (ABS(HP(2, HT(j, i))-bottom)<1.0D-12 .AND. ABS(HP(2, HT(j + 1, i))-bottom)<1.0D-12) ) THEN

          count = count + 1
          E_Dirichlet_temp(1, count) = HT(j, i)
          E_Dirichlet_temp(2, count) = HT(j + 1, i)
          E_Dirichlet_temp(5, count) = i
          IF (j == 1) THEN                        !bottom edge
            E_Dirichlet_temp(3, count) = 0
            E_Dirichlet_temp(4, count) = -1

          ELSEIF (j == 2) THEN                    !right edge
            E_Dirichlet_temp(3, count) = 1
            E_Dirichlet_temp(4, count) = 0

          ELSEIF (j == 3) THEN                    !top edge
            E_Dirichlet_temp(3, count) = 0
            E_Dirichlet_temp(4, count) = 1

          END IF
          
          node1(1:2) = HP(1:2, HT(j, i))
          node2(1:2) = HP(1:2, HT(j + 1, i))
          
            CALL Edge_in_Boundary_2D(node1, node2, IN1, IN2)
          
            !If (IN1 /= IN2) Then
            !pause
            !write(*,*) 'IN1 /= IN2'
            !End If
          
            IF ( IN1 /= -10 .AND. IN2 /= -10) THEN
                E_Dirichlet_temp(6, count) = bc_index(IN1)      !WSY REVISE 2022 7 24
                E_Dirichlet_temp(7, count) = bc_value(IN1)
                E_DirichletValueTemp(count) = bc_value(IN1)   !LY modification, 2022-7-25
            ENDIF
          
        End If

      ELSE

        IF ( (ABS(HP(1, HT(4, i))-left)<1.0D-12 .AND. ABS(HP(1, HT(1, i))-left)<1.0D-12) .OR. &
             (ABS(HP(1, HT(4, i))-right)<1.0D-12 .AND. ABS(HP(1, HT(1, i))-right)<1.0D-12) ) then

          count = count + 1                       !left edge
          E_Dirichlet_temp(1, count) = HT(4, i)
          E_Dirichlet_temp(2, count) = HT(1, i)
          E_Dirichlet_temp(5, count) = i
          E_Dirichlet_temp(3, count) = -1
          E_Dirichlet_temp(4, count) = 0

            
          node1(1:2) = HP(1:2, HT(4, i))
          node2(1:2) = HP(1:2, HT(1, i))
          
          CALL Edge_in_Boundary_2D(node1, node2, IN1, IN2)
          
            If (IN1 /= IN2) Then
                PRINT*, 'Case3: Check Partition and Code, Stop'
                pause
                write(*,*) 'IN1 /= IN2'
            End If
          
            IF ( IN1 /= -10 .AND. IN2 /= -10) THEN
                E_Dirichlet_temp(6, count) = bc_index(IN1)                        !WSY REVISE 2022 7 24
                E_Dirichlet_temp(7, count) = bc_value(IN1)
                E_DirichletValueTemp(count) = bc_value(IN1)   !LY modification, 2022-7-25
            ENDIF
          
        END IF
        
      END IF
        
        
      
    END DO
  END IF
END DO


ALLOCATE(E_Dirichlet(7, count))
Allocate(E_DirichletValue(count))

DO i = 1, count
    E_Dirichlet(:, i) = E_Dirichlet_temp(:, i)
    E_DirichletValue(i) = E_DirichletValueTemp(i)
END DO

DEALLOCATE(E_Dirichlet_temp, E_DirichletValueTemp)

END