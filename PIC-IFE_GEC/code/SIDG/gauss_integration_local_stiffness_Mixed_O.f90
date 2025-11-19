SUBROUTINE gauss_integration_local_stiffness_Mixed_O(delta, Global_Beta, information1_vector, information2_vector, &
                                                Gauss_coefficient_local_1D_linear,Gauss_point_local_1D_linear, &
	                                            this_flag1, vertices1, vertices2, trial_basis_type, trial_basis_index, &
											    test_basis_type, test_basis_index, location, this_basic, &
                                                next_basic, con_penalty, int_value, element_Gauss)
   
!=========LY modification, 2022-7-25=========
Use IMPIC_Data_2D
!=========LY modification, 2022-7-25=========
USE Gauss_Data
USE IFE_MAIN_PARAM !LY REVISE, 2022-1-6, Float number error
USE IFE_INTERFACE, ONLY: Generate_Gauss_local_triangle, Retangular_local_basis_IFE, Retangular_local_basis
IMPLICIT NONE

! consider the brother element

INTEGER                     :: delta
REAL(8), DIMENSION(:), POINTER	        	::	  Global_Beta
INTEGER                     :: information1_vector(18)
REAL(8)                     :: information2_vector(8)
REAL(8), DIMENSION(:)       :: Gauss_coefficient_local_1D_linear
REAL(8), DIMENSION(:, :)    :: Gauss_point_local_1D_linear
INTEGER, DIMENSION(4, 2)    :: this_flag1
REAL(8), DIMENSION(2, 4)    :: vertices1, vertices2
INTEGER                     :: trial_basis_type, trial_basis_index, test_basis_index, test_basis_type, location
REAL(8)                     :: con_penalty, int_value
REAL(8),DIMENSION(2)        :: h_partition1, h_partition2
INTEGER                     :: Gpn, out_index, in_index, piece_flag, nx, ny, i, this_basic, next_basic
REAL(8)                     :: temp1, temp2, temp3, temp4, temp5, temp6, x, y, trial_part, test_part, &
                                trial_part_S, test_part_S, trial_part_P, test_part_P
REAL(8)                     :: beta1, beta2, Beta   !LY REVISE, 2022-1-6, Float number error
!=========LY modification, 2022-7-25=========
Integer :: element_Gauss
Real(8) :: function_coefficient_impic1, function_coefficient_impic2
Real(8) :: function_coefficient_impic3, function_coefficient_impic4
!=========LY modification, 2022-7-25=========

Gpn = size(Gauss_coefficient_local_1D_linear)
out_index = -1
in_index = -2
int_value = 0
beta1 = information2_vector(1)
beta2 = information2_vector(2)
h_partition1(1) = vertices1(1, 2) - vertices1(1, 1)
h_partition1(2) = vertices1(2, 4) - vertices1(2, 1)
h_partition2(1) = vertices2(1, 2) - vertices2(1, 1)
h_partition2(2) = vertices2(2, 4) - vertices2(2, 1)

!LY REVISE, 2022-1-6, Float number error
!initial code is 'beta2 == Global_Beta(1)'
if  (ABS(beta2-Global_Beta(1)) < SmallValue) then
    
    if (location == out_index) then
        piece_flag = 2
        Beta = beta2
    elseif (location == in_index) then
        piece_flag = 1
        Beta = beta1
    end if
    
!initial code is 'beta2 == Global_Beta(2)'
elseif (ABS(beta2-Global_Beta(2)) < SmallValue) then
    
    if (location == out_index) then
        piece_flag = 1
        Beta = beta1
    elseif (location == in_index) then
        piece_flag = 2
        Beta = beta2
    end if
    
else
    WRITE(6,*) 'error in gauss_integration_local_stiffness_IDG_D'
    WRITE(6,*) 'beta1 =', beta1, 'beta2 =', beta2
    STOP
end if

do i = 1, Gpn
    
    x = Gauss_point_local_1D_linear(1, i)
    y = Gauss_point_local_1D_linear(2, i)
    nx = this_flag1(2, 1)
    ny = this_flag1(3, 1)
    
    !=========LY modification, 2022-7-25=========
    if (this_basic == 1) then
        CALL Retangular_local_basis_IFE(x, y, vertices1, information1_vector, information2_vector, &
                                    piece_flag, trial_basis_type, trial_basis_index, 1, 0, temp1, element_Gauss)    !trial_basis_x 
    
        CALL Retangular_local_basis_IFE(x, y, vertices1, information1_vector, information2_vector, &
                                    piece_flag, trial_basis_type, trial_basis_index, 0, 1, temp2, element_Gauss)    !trial_basis_y
    
        CALL Retangular_local_basis_IFE(x, y, vertices1, information1_vector, information2_vector, &
                                    piece_flag, trial_basis_type, trial_basis_index, 0, 0, temp3, element_Gauss)    !trial_basis_0
    
    elseif (this_basic == 0) then
        
        CALL Retangular_local_basis(x, y, vertices1(:, 1), h_partition1, &
                            trial_basis_type,trial_basis_index, 1, 0, temp1)    !trial_basis_x
    
        CALL Retangular_local_basis(x, y, vertices1(:, 1), h_partition1, &
                            trial_basis_type,trial_basis_index, 0, 1, temp2)    !trial_basis_y
    
        CALL Retangular_local_basis(x, y, vertices1(:, 1), h_partition1, &
                            trial_basis_type,trial_basis_index, 0, 0, temp3)    !trial_basis_0
        
    endif
    
    
    if (next_basic == 1) then
        
        CALL Retangular_local_basis_IFE(x, y, vertices2, information1_vector, information2_vector, &
                                    piece_flag, test_basis_type, test_basis_index, 1, 0, temp4, element_Gauss)      !test_basis_x
    
        CALL Retangular_local_basis_IFE(x, y, vertices2, information1_vector, information2_vector, &
                                    piece_flag, test_basis_type, test_basis_index, 0, 1, temp5, element_Gauss)      !test_basis_y
    
        CALL Retangular_local_basis_IFE(x, y, vertices2, information1_vector, information2_vector, &
                                    piece_flag, test_basis_type, test_basis_index, 0, 0, temp6, element_Gauss)      !test_basis_0
        
        
    elseif (next_basic == 0) then
        
        CALL Retangular_local_basis(x, y, vertices2(:, 1), h_partition2, &
                            test_basis_type,test_basis_index, 1, 0, temp4)      !test_basis_x
    
        CALL Retangular_local_basis(x, y, vertices2(:, 1), h_partition2, &
                            test_basis_type,test_basis_index, 0, 1, temp5)      !test_basis_y
    
        CALL Retangular_local_basis(x, y, vertices2(:, 1), h_partition2, &
                            test_basis_type,test_basis_index, 0, 0, temp6)      !test_basis_0
        
    endif
    
    Call Function_coefficient_2D(function_coefficient_impic1, Beta, &
                                x, y, 1, 0, 1, 0, element_Gauss)    !Chi11_Bfield
        
    Call Function_coefficient_2D(function_coefficient_impic2, Beta, &
                                x, y, 0, 1, 0, 1, element_Gauss)    !Chi22_Bfield
                                        
    Call Function_coefficient_2D(function_coefficient_impic3, Beta, &
                                x, y, 0, 1, 1, 0, element_Gauss)    !Chi12_Bfield
        
    Call Function_coefficient_2D(function_coefficient_impic4, Beta, &
                                x, y, 1, 0, 0, 1, element_Gauss)    !Chi21_Bfield
    
    If (Bfiled_index .AND. IMPIC_index) Then  !Bfield and Impic
        trial_part = ((function_coefficient_impic1 * temp1 + function_coefficient_impic3 * temp2) * nx + &
                      (function_coefficient_impic4 * temp1 + function_coefficient_impic2 * temp2) * ny) / 2.0

        test_part = temp6 * this_flag1(4, 2)

        trial_part_S = temp3 * this_flag1(4, 1)

        test_part_S = ((function_coefficient_impic1 * temp4 + function_coefficient_impic3 * temp5) * nx + &
                       (function_coefficient_impic4 * temp4 + function_coefficient_impic2 * temp5) * ny) / 2.0

    Else    !No-Bfield and Impic OR Expic
        trial_part = function_coefficient_impic1 * (temp1 * nx + temp2 * ny) / 2.0
        test_part = temp6 * this_flag1(4, 2)

        trial_part_S = temp3 * this_flag1(4, 1)
        test_part_S = function_coefficient_impic1 * (temp4 * nx + temp5 * ny) / 2.0
    End If
    
    trial_part_P = temp3 * this_flag1(4, 1)
    test_part_P = temp6 * this_flag1(4, 2)
    
    if (delta == 0) then

        int_value = int_value - trial_part * test_part * Gauss_coefficient_local_1D_linear(i) - &
                                trial_part_S * test_part_S * Gauss_coefficient_local_1D_linear(i) + &
                                trial_part_P * test_part_P * Gauss_coefficient_local_1D_linear(i) * con_penalty

    Else If (delta == 1) Then

        int_value = int_value - y * trial_part * test_part * Gauss_coefficient_local_1D_linear(i) - &
                                y * trial_part_S * test_part_S * Gauss_coefficient_local_1D_linear(i) + &
                                y * trial_part_P * test_part_P * Gauss_coefficient_local_1D_linear(i) * con_penalty

    endif
    !=========LY modification, 2022-7-25=========
    
enddo

END