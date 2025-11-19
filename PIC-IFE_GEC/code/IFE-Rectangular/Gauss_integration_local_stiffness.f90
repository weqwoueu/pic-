SUBROUTINE Gauss_integration_local_stiffness(delta, coefficient_function_name,Gauss_coefficient_local,Gauss_point_local, &
	                                         left_lower_point,h_partition,trial_basis_type, trial_basis_index, &
											 trial_derivative_degree_x,trial_derivative_degree_y, test_basis_type, &
											 test_basis_index,test_derivative_degree_x,test_derivative_degree_y, r, element_Gauss)

USE IFE_INTERFACE, ONLY: Retangular_local_basis
IMPLICIT NONE

INTEGER, DIMENSION(:,:), POINTER			::    t_c
REAL(8)                                           coefficient_function_name
REAL(8)                                           coefficient_function_name_impic !!!bjw add for impic 2019-6-3
INTEGER                                           trial_basis_type, trial_derivative_degree_x,trial_derivative_degree_y, &
												  test_basis_type,test_derivative_degree_x,test_derivative_degree_y

INTEGER                                           i, j ,k
REAL(8)								        	  r, temp1, temp2
REAL(8)	                                          h_partition(2)

REAL(8),DIMENSION(:)		                ::    Gauss_coefficient_local
REAL(8),DIMENSION(:,:)			            ::    Gauss_point_local
REAL(8)                           			::    left_lower_point(2)
REAL(8), DIMENSION(:,:), POINTER	        ::	  p_basic

INTEGER                                           trial_basis_index
INTEGER                                           test_basis_index
INTEGER                                           Gpn
INTEGER                                     ::	  delta

!==========LY modification, 2022-7-25==========
Integer :: element_Gauss
!==========LY modification, 2022-7-25==========

Gpn=SIZE(Gauss_coefficient_local)
r=0
DO i = 1, Gpn

 !  IF (delta == 0) THEN
 !
	!CALL Retangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),left_lower_point,h_partition, &
	!                            trial_basis_type,trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y, temp1)
	!
 !   CALL Retangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),left_lower_point,h_partition, &
 !                             	test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y, temp2)
 !
 !   !!! ************************ bjw add for impic 2019-6-3 **********************************************
 !   CALL Function_coefficient_2D(coefficient_function_name_impic,coefficient_function_name, &
 !                                Gauss_point_local(i,1),Gauss_point_local(i,2), &
 !                                trial_derivative_degree_x,trial_derivative_degree_y, &
 !                                test_derivative_degree_x,test_derivative_degree_y)
 !
 !   r=r+Gauss_coefficient_local(i) * coefficient_function_name_impic * temp1 *temp2
 !   !!! ************************ bjw add for impic 2019-6-3 **********************************************
 !   
 !  	 !r=r+Gauss_coefficient_local(i) * coefficient_function_name * temp1 *temp2 
 !  ELSEIF (delta ==1) THEN
 !
	!CALL Retangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),left_lower_point,h_partition, &
	!                            trial_basis_type,trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y, temp1)
 !
 !   CALL Retangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),left_lower_point,h_partition, &
 !                             	test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y, temp2)
 !
 !
 !  	 r=r+Gauss_coefficient_local(i) * coefficient_function_name * temp1 *temp2* Gauss_point_local(i,2)
 !
 !
 !  ELSE 
	!	  PRINT*, ' delta value is wrong, check again, stop'
	!	  STOP
 !  ENDIF
   
   !$ ====================================== mb.ZWZ for impic and code simplification 2021/4/19 ======================================= \\
    
    !$ trial evalution on local basis function 
    CALL Retangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),left_lower_point,h_partition, &
	                            trial_basis_type,trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y, temp1)  
    !$ test evalution on local basis function
    CALL Retangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),left_lower_point,h_partition, &
                                    test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y, temp2)       
    !$ coefficient function evaluation on gauss nodes
    CALL Function_coefficient_2D(coefficient_function_name_impic,coefficient_function_name, &
                                 Gauss_point_local(i,1),Gauss_point_local(i,2), &
                                 trial_derivative_degree_x,trial_derivative_degree_y, &
                                 test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
    
    IF (delta == 0) THEN
   	    r=r+Gauss_coefficient_local(i) * coefficient_function_name_impic * temp1 *temp2     
    ELSEIF (delta ==1) THEN      
   	    r=r+Gauss_coefficient_local(i) * coefficient_function_name_impic * temp1 *temp2* Gauss_point_local(i,2)
    ELSE 
		  PRINT*, ' delta value is wrong, check again, stop'
		  STOP
    ENDIF
    !$ ===================================== mb.ZWZ for impic and code simplification 2021/4/19 ======================================= //

ENDDO


END