SUBROUTINE Edge_in_Boundary_2D(Node1, Node2, IN1, IN2)

USE IFE_Boundary

IMPLICIT NONE

REAL(8), INTENT(IN)		::	Node1(2), Node2(2)
INTEGER, INTENT(OUT)	::	IN1, IN2

LOGICAL	  IN_LINE1, IN_LINE2
LOGICAL	  IN_LINE1_1, IN_LINE1_2
LOGICAL	  IN_LINE2_1, IN_LINE2_2
INTEGER   n_bound, i
REAL(8)  :: P1(2), P2(2), Cross1, Cross2, SmallV 


SmallV = 1.0D-5

n_bound		=   SIZE(bc_index)

IN1 = -10
IN2 = -10

IN_LINE1 = .FALSE.
IN_LINE2 = .FALSE.


DO i = 1, n_bound

   P1 =	Node1 - bc_point_1(:, i)
   P2 =	Node1 - bc_point_2(:, i)  
   Cross1 = P1(1)*P2(2) - P1(2)*P2(1)

   P1 =	Node2 - bc_point_1(:, i)
   P2 =	Node2 - bc_point_2(:, i)  
   Cross2 = P1(1)*P2(2) - P1(2)*P2(1)

   IF (IN1 == -10) CALL Node_in_Boundary_2D(Node1, bc_point_1(:, i), bc_point_2(:, i), IN_LINE1)
   IF (IN1 == -10 .AND. IN_LINE1.AND.(ABS(Cross1)<=SmallV).AND.(ABS(Cross2)<=SmallV))   THEN
      IN1 = i
   ENDIF

   IF (IN2 == -10) CALL Node_in_Boundary_2D(Node2, bc_point_1(:, i), bc_point_2(:, i), IN_LINE2)
   IF (IN2 == -10 .AND. IN_LINE2.AND.(ABS(Cross1)<=SmallV).AND.(ABS(Cross2)<=SmallV))   THEN
      IN2 = i
   ENDIF

    IF (IN1 /= -10 .AND. IN2 /= -10)  THEN
	   IF (IN1 == IN2) EXIT
	   IF (IN1 /= IN2) THEN
  		  IN_LINE1_1 = .FALSE.
          IN_LINE1_2 = .FALSE.
		  CALL Node_in_Boundary_2D(Node1, bc_point_1(:, IN1), bc_point_2(:, IN1), IN_LINE1_1)
		  CALL Node_in_Boundary_2D(Node2, bc_point_1(:, IN1), bc_point_2(:, IN1), IN_LINE1_2)
		  IN_LINE2_1 = .FALSE.
          IN_LINE2_2 = .FALSE.
		  CALL Node_in_Boundary_2D(Node1, bc_point_1(:, IN2), bc_point_2(:, IN2), IN_LINE2_1)
		  CALL Node_in_Boundary_2D(Node2, bc_point_1(:, IN2), bc_point_2(:, IN2), IN_LINE2_2)
		  IF(IN_LINE1_1.AND.IN_LINE1_2.AND.IN_LINE2_1.AND.(.NOT.IN_LINE2_2)) IN2=IN1
		  IF(IN_LINE1_1.AND.IN_LINE1_2.AND.(.NOT.IN_LINE2_1).AND.IN_LINE2_2) IN2=IN1
		  IF(IN_LINE1_1.AND.(.NOT.IN_LINE1_2).AND.IN_LINE2_1.AND.IN_LINE2_2) IN1=IN2
		  IF((.NOT.IN_LINE1_1).AND.IN_LINE1_2.AND.IN_LINE2_1.AND.IN_LINE2_2) IN1=IN2
		  EXIT
	   ENDIF
	ENDIF



ENDDO



END