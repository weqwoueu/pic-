SUBROUTINE generate_DG_node_flag(P, T, P_average, P_flag)
!DATE: 2021-12-3
!FUNCTION: LY REVISE FOR NODE FLAG TO DIFFERENT CG OR DG NODE

IMPLICIT NONE
REAL(8), DIMENSION(:,:), POINTER :: P, P_flag, P_average
INTEGER, DIMENSION(:,:), POINTER :: T

INTEGER :: num_of_element, n, i, num

ALLOCATE(P_flag(5,SIZE(P,2)))
P_flag = 0
DO n = 1, SIZE(P,2)
  P_flag(1:2, n) = P(1:2, n)
END DO

num_of_element = SIZE(T,2)

!P_flag(3,:)=1 denote the node is DG node.
DO n = 1, num_of_element
  IF (T(5,n) /= 0) THEN
    P_flag(3, T(1:4,n)) = 1
    !P_flag(5, T(1:4,n)) = n
  END IF
  P_flag(5, T(1:4,n)) = n
END DO

DO n = 1, SIZE(P_average,2)
  num = 0
  DO i = 3, 6
    IF (P_average(i, n) /= 0) THEN 
      num = num + 1
    END IF
  END DO
  
  !If any coordinate node only have two mesh node on P_average, meantime all this mesh node is DG, then the coordinate node is "hand node". 
  IF (num==2) THEN
    IF (INT(P_flag(3,INT(P_average(3,n))))==1 .AND. INT(P_flag(3,INT(P_average(4,n))))==1) THEN
      P_flag(4,INT(P_average(3,n))) = 1
      P_flag(4,INT(P_average(4,n))) = 1
    END IF
  END IF
END DO

!P_flag:
!P_flag(1:2, i): the coordinates of i node
!P_flag(3, i): the flag of CG or DG, CG--0, DG--1
!P_flag(4, i): the flag of 'hand' node, no--0, yes--1
!P_flag(5, i): the element index of DG node
END SUBROUTINE generate_DG_node_flag