SUBROUTINE gauss_integration_local_stiffness_DG_FE_O(delta, function_coefficient,Gauss_coefficient_local_1D_linear, &
                                                Gauss_point_local_1D_linear, &
	                                            this_flag1, vertices, brother_vertices, trial_basis_type, trial_basis_index, &
											    test_basis_type, test_basis_index, con_penalty, int_value, element_Gauss)

USE IFE_INTERFACE, ONLY: Retangular_local_basis
!=========LY modification, 2022-7-25=========
Use IMPIC_Data_2D
!=========LY modification, 2022-7-25=========
IMPLICIT NONE
    
REAL(8)                                           function_coefficient
INTEGER                                           trial_basis_type, test_basis_type

INTEGER                                           i, j ,k
REAL(8)                                           int_value

REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_local_1D_linear
REAL(8),DIMENSION(:,:)			            ::    Gauss_point_local_1D_linear
REAL(8),DIMENSION(2,4)                      ::    vertices, brother_vertices
INTEGER,DIMENSION(4,2)                      ::    this_flag1
REAL(8)                                     ::    x, y, temp1, temp2, temp3, temp4, temp5, temp6, con_penalty
REAL(8),DIMENSION(2)                    	::    h_partition, h_partition_br
REAL(8)                               ::    trial_part, test_part, test_part_S, trial_part_S, trial_part_P, test_part_P

INTEGER                                     ::    trial_basis_index
INTEGER                                     ::    test_basis_index
INTEGER                                     ::    Gpn
INTEGER                                     ::	  delta, nx, ny 

!=========LY modification, 2022-7-25=========
Integer :: element_Gauss
Real(8) :: function_coefficient_impic1, function_coefficient_impic2
Real(8) :: function_coefficient_impic3, function_coefficient_impic4
!=========LY modification, 2022-7-25=========

Gpn = SIZE(Gauss_coefficient_local_1D_linear)
int_value = 0
nx = this_flag1(2, 1)
ny = this_flag1(3, 1)
h_partition(1) = vertices(1, 2) - vertices(1, 1)
h_partition(2) = vertices(2, 4) - vertices(2, 1)
h_partition_br(1) = brother_vertices(1, 2) - brother_vertices(1, 1)
h_partition_br(2) = brother_vertices(2, 4) - brother_vertices(2, 1)

do i = 1, Gpn
     
    x = Gauss_point_local_1D_linear(1, i)
    y = Gauss_point_local_1D_linear(2, i)
    
    CALL Retangular_local_basis(x, y, vertices(:, 1), h_partition, &
                            trial_basis_type,trial_basis_index, 1, 0, temp1)
    
    CALL Retangular_local_basis(x, y, vertices(:, 1), h_partition, &
                            trial_basis_type,trial_basis_index, 0, 1, temp2)
    
    CALL Retangular_local_basis(x, y, vertices(:, 1), h_partition, &
                            trial_basis_type,trial_basis_index, 0, 0, temp3)
    
    CALL Retangular_local_basis(x, y, brother_vertices(:, 1), h_partition_br, &
                            test_basis_type,test_basis_index, 1, 0, temp4)
    
    CALL Retangular_local_basis(x, y, brother_vertices(:, 1), h_partition_br, &
                            test_basis_type,test_basis_index, 0, 1, temp5)
    
    CALL Retangular_local_basis(x, y, brother_vertices(:, 1), h_partition_br, &
                            test_basis_type,test_basis_index, 0, 0, temp6)
    
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

    Else

        Write(*,*) 'wrong delta pleas check'
        stop
        
    end if
    !=========LY modification, 2022-7-25=========
    
enddo

END
