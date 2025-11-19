SUBROUTINE Boundary_Condition_info_2D  ( n_grid, t_basic, p_basic, e_basic)


USE IFE_MAIN_PARAM
USE IFE_Boundary

IMPLICIT NONE

INTEGER, DIMENSION(:,:), POINTER			::	t_basic
INTEGER, DIMENSION(:,:), POINTER	        ::	e_basic, e_basic_tmp
REAL(8), DIMENSION(:,:), POINTER			::	p_basic


INTEGER		n_elements, n_nodes, e, eg, n_grid
INTEGER     el_edge, edge_count
INTEGER 	IN1, IN2



n_elements	=	n_grid
n_nodes		=	SIZE(p_basic,2)


el_edge     = 4

ALLOCATE(e_basic_tmp(5, n_elements*el_edge))
e_basic_tmp = -10


edge_count = 0

DO e=1,n_elements

  DO eg=1, el_edge

    IF( eg /= el_edge) THEN

      CALL Edge_in_Boundary_2D(p_basic(:,t_basic(eg,e)), p_basic(:,t_basic(eg+1,e)), IN1, IN2)


      IF ( IN1 /= -10 .AND. IN2 /= -10) THEN
        edge_count = edge_count + 1
        e_basic_tmp(1, edge_count) = t_basic(eg,e) !eg
        e_basic_tmp(2, edge_count) = t_basic(eg+1,e) !eg+1
        e_basic_tmp(3, edge_count) = e
        e_basic_tmp(4, edge_count) = IN1
        e_basic_tmp(5, edge_count) = IN2
      ENDIF
      ! 		      IF(IN1 /= -10 .AND. IN2 == -10) THEN
      !			     PRINT*, 'It is single Node on the Bounday', IN1, 'Please Check it carefully'
      !				 PRINT*, 'Nodes Coordinate'
      !				 PRINT*, p_basic(:,t_basic(eg,e)), p_basic(:,t_basic(eg+1,e))
      !				 PRINT*, 'Bounday Information'
      !				 PRINT*, bc_point_1(:, IN1), bc_point_2(:, IN1)
      !				 PRINT*
      !		      ENDIF
      ! 		      IF(IN1 == -10 .AND. IN2 /= -10) THEN
      !			     PRINT*, 'It is single Node on the Bounday', IN2, 'Please Check it carefully'
      !				 PRINT*, 'Nodes Coordinate'
      !				 PRINT*, p_basic(:,t_basic(eg,e)), p_basic(:,t_basic(eg+1,e))
      !				 PRINT*, 'Bounday Information'
      !				 PRINT*, bc_point_1(:, IN2), bc_point_2(:, IN2)
      !				 PRINT*
      !		      ENDIF
    ELSE


      CALL Edge_in_Boundary_2D(p_basic(:,t_basic(eg,e)), p_basic(:,t_basic(1,e)), IN1, IN2)


      IF ( IN1 /= -10 .AND. IN2 /= -10) THEN
        edge_count = edge_count + 1
        e_basic_tmp(1, edge_count) = t_basic(eg,e) !eg
        e_basic_tmp(2, edge_count) = t_basic(1,e) !1
        e_basic_tmp(3, edge_count) = e
        e_basic_tmp(4, edge_count) = IN1
        e_basic_tmp(5, edge_count) = IN2
      ENDIF
      ! 		      IF(IN1 /= -10 .AND. IN2 == -10) THEN
      !			     PRINT*, 'It is single Node on the Bounday', IN1, 'Please Check it carefully'
      !				 PRINT*, 'Nodes Coordinate'
      !				 PRINT*, p_basic(:,t_basic(eg,e)), p_basic(:,t_basic(1,e))
      !				 PRINT*, 'Bounday Information'
      !				 PRINT*, bc_point_1(:, IN1), bc_point_2(:, IN1)
      !				 PRINT*
      !		      ENDIF
      ! 		      IF(IN1 == -10 .AND. IN2 /= -10) THEN
      !			     PRINT*, 'It is single Node on the Bounday', IN2, 'Please Check it carefully'
      !				 PRINT*, 'Nodes Coordinate'
      !				 PRINT*, p_basic(:,t_basic(eg,e)), p_basic(:,t_basic(1,e))
      !				 PRINT*, 'Bounday Information'
      !				 PRINT*, bc_point_1(:, IN2), bc_point_2(:, IN2)
      !				 PRINT*
      !		      ENDIF
    ENDIF

  ENDDO
ENDDO

ALLOCATE(e_basic(5, edge_count))

e_basic = e_basic_tmp(:, 1:edge_count)

DEALLOCATE(e_basic_tmp)



END