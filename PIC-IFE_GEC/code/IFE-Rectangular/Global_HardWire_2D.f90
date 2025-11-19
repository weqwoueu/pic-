SUBROUTINE Global_HardWire_2D(	Coeff_FUN, p_basic, t_basic_int, element_index	,			&
							    G_Stiff_HW, G_RHS_HW, G_Mass_HW, G_Interpol_HW, delta)

USE IFE_MAIN_PARAM

IMPLICIT NONE

EXTERNAL								Coeff_FUN
REAL(8), DIMENSION(:,:), POINTER	::	p_basic
INTEGER, DIMENSION(:,:), POINTER	::	t_basic_int

REAL(8)							::	G_Stiff_HW(2,4,4), G_RHS_HW(2,4,4), &
									G_Mass_HW(2,4,4), G_Interpol_HW(4,4)

REAL(8)							::	vert(2,4), El_Stiff_HW(4,4), El_RHS_HW(4,4), &
									El_Mass_HW(4,4), El_Interpol_HW(4,4)
INTEGER							::	e, i, j, k, n_nodes_in_elem, n_elem, n_elem_hw
INTEGER                         ::  delta
!===========================NEW ADD===================================================
INTEGER, DIMENSION(:), POINTER				::	  element_index
!======================================================================================

n_nodes_in_elem = 4	! Retangle

n_elem = SIZE(t_basic_int,2)
k = 0 

DO e=1, n_elem

   IF ((element_index(e) < -1 .OR. element_index(e) == -1).AND.(k==0)) THEN

      n_elem_hw =e          !> find the first non-interface element No.
	  k = k+1

  ENDIF
ENDDO

DO e=n_elem_hw,n_elem_hw
	El_Interpol_HW = Zero

    IF (element_index(e) < -1 .OR. element_index(e) == -1) THEN     ! Non-interface element

		vert = p_basic(:, t_basic_int(1:4,e));

		CALL FE_Interpol_HW_2D( vert, El_Interpol_HW )

		DO i=1,4
			DO j=1,4
				G_Interpol_HW(i,j) = El_Interpol_HW(i,j)
			END DO
		END DO
	END IF
END DO

END