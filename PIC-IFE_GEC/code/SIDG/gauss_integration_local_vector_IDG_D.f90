SUBROUTINE gauss_integration_local_vector_IDG_D(delta, Global_Beta, information1_vector, information2_vector, &
                                                Gauss_coefficient_local_1D_linear, Gauss_point_local_1D_linear, &
                                                this_flag1, vertices, test_basis_type, test_basis_index, location, &
                                                con_penalty, boundary_value, int_value, element_Gauss)
  
! -------------------LY add vector of IDG edge for IDG solver 2021-9-28-------------------
!=========LY modification, 2022-7-25=========
Use IMPIC_Data_2D
!=========LY modification, 2022-7-25=========
USE Gauss_Data
USE IFE_MAIN_PARAM !LY REVISE, 2022-1-6, Float number error
USE IFE_INTERFACE, ONLY: Generate_Gauss_local_triangle, Retangular_local_basis_IFE

IMPLICIT NONE

INTEGER                        :: delta
REAL(8), DIMENSION(:), POINTER :: Global_Beta
INTEGER                        :: information1_vector(18)
REAL(8)                        :: information2_vector(8)
REAL(8), DIMENSION(:)          :: Gauss_coefficient_local_1D_linear
REAL(8), DIMENSION(:, :)       :: Gauss_point_local_1D_linear
INTEGER, DIMENSION(4, 2)       :: this_flag1
REAL(8), DIMENSION(2, 4)       :: vertices
INTEGER                        :: test_basis_type, test_basis_index, location
Real(8) :: boundary_value
REAL(8)                        :: con_penalty, int_value

INTEGER                        :: Gpn, out_index, in_index, piece_flag, nx, ny, i, &
                                    derivative_degree_x, derivative_degree_y
REAL(8)                        :: x, y, temp1, temp2, temp3, g_value, index, test_part_S, test_part_P
REAL(8)                        :: beta1, beta2   !LY REVISE, 2022-1-6, Float number error
!=========LY modification, 2022-4-25=========
Integer :: element_Gauss
Real(8) :: function_coefficient_impic1, function_coefficient_impic2
Real(8) :: function_coefficient_impic3, function_coefficient_impic4
Real(8) :: Beta
!=========LY modification, 2022-4-25=========

piece_flag = 0
Beta = 0
nx = 0
ny = 0
i = 0
x = 0.0
y = 0.0
temp1 = 0.0
temp2 = 0.0
temp3 = 0.0
g_value = 0.0
test_part_S = 0.0
test_part_P = 0.0


Gpn = SIZE(Gauss_coefficient_local_1D_linear)
beta1 = information2_vector(1)
beta2 = information2_vector(2)
out_index = -1
in_index = -2

derivative_degree_x = 0
derivative_degree_y = 0
g_value = 0.0
index = 0
!index should be equal Global_Beta(1) or Global_Beta(2).

int_value = 0.0

!LY REVISE, 2022-1-6, Float number error
!initial code is 'beta2 == Global_Beta(1)'
IF  (ABS(beta2-Global_Beta(1)) < SmallValue) THEN

  IF (location == out_index) THEN
    piece_flag = 2
    Beta = beta2
  ELSE IF (location == in_index) THEN
    piece_flag = 1
    Beta = beta1
  END IF
  
!initial code is 'beta2 == Global_Beta(2)'
ELSEIF (ABS(beta2-Global_Beta(2)) < SmallValue) THEN

  IF (location == out_index) THEN
    piece_flag = 1
    Beta = beta1
  ELSE IF (location == in_index) THEN
    piece_flag = 2
    Beta = beta2
  END IF
  
ELSE
  WRITE(6,*) 'error in gauss_integration_local_vector_IDG_D'
END IF

DO i = 1, Gpn
  
  x = Gauss_point_local_1D_linear(1, i)
  y = Gauss_point_local_1D_linear(2, i)
  nx = this_flag1(2, 1)
  ny = this_flag1(3, 1)
  
  !=========LY modification, 2022-7-25=========
  CALL Retangular_local_basis_IFE(x, y, vertices, information1_vector, information2_vector, &
                                  piece_flag, test_basis_type, test_basis_index, 1, 0, temp1, element_Gauss)  !test_basis_x
  
  CALL Retangular_local_basis_IFE(x, y, vertices, information1_vector, information2_vector, &
                                  piece_flag, test_basis_type, test_basis_index, 0, 1, temp2, element_Gauss)  !test_basis_y
  
  CALL Retangular_local_basis_IFE(x, y, vertices, information1_vector, information2_vector, &
                                  piece_flag, test_basis_type, test_basis_index, 0, 0, temp3, element_Gauss)  !test_basis_0
  
  Call Function_coefficient_2D(function_coefficient_impic1, Beta, &
                              x, y, 1, 0, 1, 0, element_Gauss)    !Chi11_Bfield
    
  Call Function_coefficient_2D(function_coefficient_impic2, Beta, &
                              x, y, 0, 1, 0, 1, element_Gauss)    !Chi22_Bfield
                                  
  Call Function_coefficient_2D(function_coefficient_impic3, Beta, &
                              x, y, 0, 1, 1, 0, element_Gauss)    !Chi12_Bfield
    
  Call Function_coefficient_2D(function_coefficient_impic4, Beta, &
                              x, y, 1, 0, 0, 1, element_Gauss)    !Chi21_Bfield
  
  If (Bfiled_index .AND. IMPIC_index) Then  !Bfield and Impic

    test_part_S = (function_coefficient_impic1 * temp1 + function_coefficient_impic3 * temp2) * nx + &
                  (function_coefficient_impic4 * temp1 + function_coefficient_impic2 * temp2) * ny 

  Else  !No-Bfield and Impic OR Expic
    test_part_S = function_coefficient_impic1 * (temp1 * nx + temp2 * ny)
  End If
  
  !=========IFE=========
  !CALL Function_True(delta, index, x, y, derivative_degree_x, derivative_degree_y, g_value)
  !=========IFE=========
  !=========PIC=========
  g_value = boundary_value
  !=========PIC=========
  
  test_part_P = temp3
  
  IF (delta == 0) THEN

    int_value = int_value - (test_part_S * g_value * Gauss_coefficient_local_1D_linear(i)) + &
                            (test_part_P * g_value * Gauss_coefficient_local_1D_linear(i) * con_penalty)

  Else If (delta == 1) Then
    
    int_value = int_value - (y * test_part_S * g_value * Gauss_coefficient_local_1D_linear(i)) + &
                            (y * test_part_P * g_value * Gauss_coefficient_local_1D_linear(i) * con_penalty)

  END IF
  !=========LY modification, 2022-7-25=========
END DO

END SUBROUTINE