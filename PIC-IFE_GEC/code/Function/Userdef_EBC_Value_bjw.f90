SUBROUTINE Userdef_EBC_Value_bjw(delta, p_basic, t_basic_int, node_type, EBC_Value, EBC_Value_xt)


USE IFE_MAIN_PARAM
USE IFE_Data , ONLY:Global_Beta!$ ab.ZWZ for checking convergence

IMPLICIT NONE

INTEGER, INTENT(IN) :: delta        !$ ab.ZWZ 2021/7/7
REAL(8), DIMENSION(:,:), INTENT(IN)			::	p_basic
INTEGER, DIMENSION(:,:), INTENT(IN)			::	t_basic_int, node_type
REAL(8), DIMENSION(:,:), INTENT(OUT)		::	EBC_Value
REAL(8), DIMENSION(:), INTENT(OUT)		    ::	EBC_Value_xt

INTEGER		i, j, k ,e ,num_point_elem, node_index, num_elem
REAL(8)		x, y, xmin, ymin, xmax, ymax, r1, r2, beta, true_value, index

index = 2.2 ! useless
num_point_elem = 4
num_elem = SIZE(t_basic_int,2)

DO i=1,4
	Do j=1,num_elem
		EBC_Value(i,j) = Zero
	ENDDO
ENDDO


DO e=1, num_elem
	DO i=1, num_point_elem
		node_index = t_basic_int(i,e)
    IF(node_type(1,node_index) == 1)THEN
      xmin = MINVAL(p_basic(1,:))
      xmax = MAXVAL(p_basic(1,:))
      ymin = MINVAL(p_basic(2,:))
      ymax = MAXVAL(p_basic(2,:))

      IF( p_basic(1, node_index) == xmin .OR. p_basic(1, node_index) == xmax .OR. &
        p_basic(2, node_index) == ymin .OR. p_basic(2, node_index) == ymax ) THEN

        x = p_basic(1, node_index)
        y = p_basic(2, node_index)
        !EBC_Value(i,e) = x*x+y*y
        !EBC_Value(i,e) = EXP(x*x+y*y)/10.0
        !
        !EBC_Value_xt(node_index) = EXP(x*x+y*y)/10.0
        !$ ======================= mb.ZWZ 2021/7/7 ================= \\
        IF (delta == 0) THEN
          !EBC_Value(i,e) = EXP(x*x+y*y)/Global_Beta(1)
          !EBC_Value_xt(node_index) = EXP(x*x+y*y)/Global_Beta(1)
          CALL Function_True(delta, index, x, y, 0, 0, true_value)
          EBC_Value(i,e) = true_value
          EBC_Value_xt(node_index) = true_value
        ELSEIF (delta == 1) THEN
          CALL Function_True(delta, index, x, y, 0, 0, true_value)
          EBC_Value(i,e) = true_value
          EBC_Value_xt(node_index) = true_value
        ENDIF
        !$ ======================= mb.ZWZ 2021/7/7 ================= //

      END IF

    ENDIF
	ENDDO
ENDDO



END SUBROUTINE