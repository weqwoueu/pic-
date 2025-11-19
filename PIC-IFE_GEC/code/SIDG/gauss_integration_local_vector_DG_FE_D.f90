SUBROUTINE gauss_integration_local_vector_DG_FE_D(delta, function_coefficient, Gauss_coefficient_local_1D_linear, &
                                                Gauss_point_local_1D_linear, &
                    this_flag1, vertices, test_basis_type, test_basis_index, con_penalty, boundary_value, int_value, element_Gauss)

! -------------------LY add vector of IDG edge for IDG solver 2021-9-28-------------------
!=========LY modification, 2022-7-25==========
Use IMPIC_Data_2D
!=========LY modification, 2022-7-25==========
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Retangular_local_basis

IMPLICIT NONE

INTEGER                 :: delta
REAL(8)                 :: function_coefficient
REAL(8), DIMENSION(:)   :: Gauss_coefficient_local_1D_linear
REAL(8), DIMENSION(:,:) :: Gauss_point_local_1D_linear
INTEGER, DIMENSION(4,2) :: this_flag1
REAL(8), DIMENSION(2,4) :: vertices
INTEGER                 :: test_basis_type, test_basis_index
!=========LY modification, 2022-7-25=========
Real(8) :: boundary_value
!=========LY modification, 2022-7-25=========
REAL(8)                 :: con_penalty, int_value

REAL(8)                 :: x, y, temp1, temp2, temp3, g_value, index, test_part_S, test_part_P
REAL(8), DIMENSION(2)   :: h_partition
INTEGER                 :: Gpn, i, nx, ny, derivative_degree_x, derivative_degree_y
!=========LY modification, 2022-7-25=========
Integer :: element_Gauss
Real(8) :: function_coefficient_impic1, function_coefficient_impic2
Real(8) :: function_coefficient_impic3, function_coefficient_impic4
!=========LY modification, 2022-7-25=========

!x = 0.0
!y = 0.0
!temp1 = 0.0
!temp2 = 0.0
!temp3 = 0.0
!g_value = 0.0
!test_part_S = 0.0
!test_part_P = 0.0
derivative_degree_x = 0
derivative_degree_y = 0

!Gpn = 0
!nx = 0
!ny = 0
int_value = 0.0

Gpn = SIZE(Gauss_coefficient_local_1D_linear)
nx = this_flag1(2, 1)
ny = this_flag1(3, 1)
h_partition(1) = vertices(1, 2) - vertices(1, 1)
h_partition(2) = vertices(2, 4) - vertices(2, 1)

index = 0
!index should be equal Global_Beta(1) or Global_Beta(2).

DO i = 1, Gpn

  x = Gauss_point_local_1D_linear(1, i)
  y = Gauss_point_local_1D_linear(2, i)
  
  CALL Retangular_local_basis(x, y, vertices(:, 1), h_partition, &
                              test_basis_type,test_basis_index, 1, 0, temp1)
      
  CALL Retangular_local_basis(x, y, vertices(:, 1), h_partition, &
                              test_basis_type,test_basis_index, 0, 1, temp2)
      
  CALL Retangular_local_basis(x, y, vertices(:, 1), h_partition, &
                              test_basis_type,test_basis_index, 0, 0, temp3)
  
  !=========LY modification, 2022-7-25=========
  Call Function_coefficient_2D(function_coefficient_impic1, function_coefficient, &
                              x, y, 1, 0, 1, 0, element_Gauss)    !Chi11_Bfield

  Call Function_coefficient_2D(function_coefficient_impic2, function_coefficient, &
                              x, y, 0, 1, 0, 1, element_Gauss)    !Chi22_Bfield
                              
  Call Function_coefficient_2D(function_coefficient_impic3, function_coefficient, &
                              x, y, 0, 1, 1, 0, element_Gauss)    !Chi12_Bfield

  Call Function_coefficient_2D(function_coefficient_impic4, function_coefficient, &
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
  
  IF (delta==0) THEN

    int_value = int_value - (test_part_S * g_value * Gauss_coefficient_local_1D_linear(i)) + &
                            (test_part_P * g_value * Gauss_coefficient_local_1D_linear(i) * con_penalty)

  Else If (delta==1) Then
  
    int_value = int_value - (y * test_part_S * g_value * Gauss_coefficient_local_1D_linear(i)) + &
                            (y * test_part_P * g_value * Gauss_coefficient_local_1D_linear(i) * con_penalty)

  END IF
  !=========LY modification, 2022-7-25=========
END DO
  
END SUBROUTINE