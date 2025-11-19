SUBROUTINE  Gauss_integration_local_error_IFE(delta, Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
			                                      vertices,information_vector_1,information_vector_2,elem_index, &
												  test_basis_type,test_derivative_degree_x,test_derivative_degree_y,U_fe,error, element_Gauss)
    
USE IFE_Data      
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Generate_Gauss_local_triangle, Retangular_local_basis_IFE		
IMPLICIT NONE

REAL(8),DIMENSION(:)				        ::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:)		                ::    Gauss_point_reference_triangle

INTEGER                                           information_vector_1(18)
REAL(8)                                           information_vector_2(8)
REAL(8)                                           vertices(2,4)
INTEGER                                           test_basis_type,test_derivative_degree_x,test_derivative_degree_y
INTEGER                                           elem_index
INTEGER                                           test_basis_index
REAL(8)								        	  error, r1, r2, r
INTEGER                                           Gpn,i
INTEGER                                           interface_element_type
INTEGER                                           pointer_reference_to_local(4)
REAL(8)	                                          beta1, beta2, Dx, Dy, Ex, Ey
REAL(8)                                           vertices_triangle(2,3)
INTEGER                                     ::	  delta
REAL(8)                                           coefficient_function_name
REAL(8)                                           coefficient_function_name_impic !!!bjw add for impic 2019-6-3

REAL(8),DIMENSION(:),ALLOCATABLE    :: U_fe
REAL(8)                                           uh(9), u_true(9), A, temp2
INTEGER                                           derivative_degree_x, derivative_degree_y, j
INTEGER                                           piece_flag
!=========LY modification, 2022-7-25=========
Integer :: element_Gauss
!=========LY modification, 2022-7-25=========

derivative_degree_x = test_derivative_degree_x
derivative_degree_y = test_derivative_degree_y

Gpn = SIZE(Gauss_coefficient_reference_triangle) 
interface_element_type=information_vector_1(6)
pointer_reference_to_local=information_vector_1(11:14)
beta1=information_vector_2(1)
beta2=information_vector_2(2)
Dx=information_vector_2(3)
Dy=information_vector_2(4)
Ex=information_vector_2(5)
Ey=information_vector_2(6)

error = 0.0

IF (interface_element_type==1) THEN
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(1))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(1))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey

    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

    uh = 0.0
    u_true = 0.0
    piece_flag = 1
    DO i = 1, Gpn

        DO j = 1, 4

		    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                    information_vector_1,information_vector_2,piece_flag,test_basis_type,j, &
										    test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 
                
            uh(i) = uh(i) + U_fe(HT(j,elem_index)) * r2

        END DO
            
        IF (beta1 == Global_Beta(1)) THEN
            coefficient_function_name = Global_Beta(1)
        ELSE
            coefficient_function_name = Global_Beta(2)
        ENDIF
            
        CALL Function_True(delta,coefficient_function_name, Gauss_point_local_triangle_Nine(i,1),&
                            Gauss_point_local_triangle_Nine(i,2), derivative_degree_x, derivative_degree_y,temp2)
        u_true(i) = temp2

        !error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= \\
        IF (delta == 0) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        ELSEIF (delta == 1) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))* Gauss_point_local_triangle_Nine(i,2)    
        ENDIF 
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= //
        
    END DO

    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(4))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = vertices(1,pointer_reference_to_local(3))
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(4))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = vertices(2,pointer_reference_to_local(3))

    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
    uh = 0.0
    u_true = 0.0
    piece_flag = 2
    DO i = 1, Gpn

        DO j = 1, 4
                
		    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                    information_vector_1,information_vector_2,piece_flag,test_basis_type,j, &
										    test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 
            uh(i) = uh(i) + U_fe(HT(j,elem_index)) * r2
            
        END DO
            
        IF (beta2 == Global_Beta(1)) THEN
            coefficient_function_name = Global_Beta(1)
        ELSE
            coefficient_function_name = Global_Beta(2)
        ENDIF
            
        CALL Function_True(delta,coefficient_function_name, Gauss_point_local_triangle_Nine(i,1),&
                            Gauss_point_local_triangle_Nine(i,2), derivative_degree_x, derivative_degree_y,temp2)
        u_true(i) = temp2
            
        !error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= \\
        IF (delta == 0) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        ELSEIF (delta == 1) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))* Gauss_point_local_triangle_Nine(i,2)    
        ENDIF 
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= //
            
    END DO

    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(3))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
    uh = 0.0
    u_true = 0.0
    piece_flag = 2
    DO i = 1, Gpn

        DO j = 1, 4
                
		    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                    information_vector_1,information_vector_2,piece_flag,test_basis_type,j, &
										    test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 
            uh(i) = uh(i) + U_fe(HT(j,elem_index)) * r2
            
        END DO
            
        IF (beta2 == Global_Beta(1)) THEN
            coefficient_function_name = Global_Beta(1)
        ELSE
            coefficient_function_name = Global_Beta(2)
        ENDIF
            
        CALL Function_True(delta,coefficient_function_name, Gauss_point_local_triangle_Nine(i,1),&
                            Gauss_point_local_triangle_Nine(i,2), derivative_degree_x, derivative_degree_y,temp2)
        u_true(i) = temp2
            
        !error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= \\
        IF (delta == 0) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        ELSEIF (delta == 1) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))* Gauss_point_local_triangle_Nine(i,2)    
        ENDIF 
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= //

    END DO

    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(2))
	vertices_triangle(1,2) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(2))
	vertices_triangle(2,2) = vertices(2,pointer_reference_to_local(3))
	vertices_triangle(2,3) = Ey
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine)  
    uh = 0.0
    u_true = 0.0
    piece_flag = 2
    DO i = 1, Gpn

        DO j = 1, 4
                
		    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                    information_vector_1,information_vector_2,piece_flag,test_basis_type,j, &
										    test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 
            uh(i) = uh(i) + U_fe(HT(j,elem_index)) * r2
            
        END DO
            
        IF (beta2 == Global_Beta(1)) THEN
            coefficient_function_name = Global_Beta(1)
        ELSE
            coefficient_function_name = Global_Beta(2)
        ENDIF
            
        CALL Function_True(delta,coefficient_function_name, Gauss_point_local_triangle_Nine(i,1),&
                            Gauss_point_local_triangle_Nine(i,2), derivative_degree_x, derivative_degree_y,temp2)
        u_true(i) = temp2
            
        !error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= \\
        IF (delta == 0) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        ELSEIF (delta == 1) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))* Gauss_point_local_triangle_Nine(i,2)    
        ENDIF 
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= //

    END DO

ELSEIF (interface_element_type==2) THEN
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(4))
	vertices_triangle(1,2) = vertices(1,pointer_reference_to_local(1))
	vertices_triangle(1,3) = Dx
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(4))
	vertices_triangle(2,2) = vertices(2,pointer_reference_to_local(1))
	vertices_triangle(2,3) = Dy
 
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
    uh = 0.0
    u_true = 0.0
    piece_flag = 1
    DO i = 1, Gpn
 
        DO j = 1, 4
                
		    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                    information_vector_1,information_vector_2,piece_flag,test_basis_type,j, &
										    test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 
            uh(i) = uh(i) + U_fe(HT(j,elem_index)) * r2
            
        END DO
            
        IF (beta1 == Global_Beta(1)) THEN
            coefficient_function_name = Global_Beta(1)
        ELSE
            coefficient_function_name = Global_Beta(2)
        ENDIF
            
        CALL Function_True(delta,coefficient_function_name, Gauss_point_local_triangle_Nine(i,1),&
                            Gauss_point_local_triangle_Nine(i,2), derivative_degree_x, derivative_degree_y,temp2)
        u_true(i) = temp2
            
        !error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= \\
        IF (delta == 0) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        ELSEIF (delta == 1) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))* Gauss_point_local_triangle_Nine(i,2)    
        ENDIF 
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= //
 
    ENDDO
    
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(1))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(1))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
    uh = 0.0
    u_true = 0.0
    piece_flag = 1
    DO i = 1, Gpn
 
        DO j = 1, 4
                
		    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                    information_vector_1,information_vector_2,piece_flag,test_basis_type,j, &
										    test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 
            uh(i) = uh(i) + U_fe(HT(j,elem_index)) * r2
            
        END DO
            
        IF (beta1 == Global_Beta(1)) THEN
            coefficient_function_name = Global_Beta(1)
        ELSE
            coefficient_function_name = Global_Beta(2)
        ENDIF
            
        CALL Function_True(delta,coefficient_function_name, Gauss_point_local_triangle_Nine(i,1),&
                            Gauss_point_local_triangle_Nine(i,2), derivative_degree_x, derivative_degree_y,temp2)
        u_true(i) = temp2
            
        !error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= \\
        IF (delta == 0) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        ELSEIF (delta == 1) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))* Gauss_point_local_triangle_Nine(i,2)    
        ENDIF 
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= //
 
    ENDDO
    
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(2))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(2))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
    uh = 0.0
    u_true = 0.0
    piece_flag = 2
    DO i = 1, Gpn
 
        DO j = 1, 4
                
		    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                    information_vector_1,information_vector_2,piece_flag,test_basis_type,j, &
										    test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 
            uh(i) = uh(i) + U_fe(HT(j,elem_index)) * r2
            
        END DO
            
        IF (beta2 == Global_Beta(1)) THEN
            coefficient_function_name = Global_Beta(1)
        ELSE
            coefficient_function_name = Global_Beta(2)
        ENDIF
            
        CALL Function_True(delta,coefficient_function_name, Gauss_point_local_triangle_Nine(i,1),&
                            Gauss_point_local_triangle_Nine(i,2), derivative_degree_x, derivative_degree_y,temp2)
        u_true(i) = temp2
            
        !error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= \\
        IF (delta == 0) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        ELSEIF (delta == 1) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))* Gauss_point_local_triangle_Nine(i,2)    
        ENDIF 
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= //         
 
    ENDDO
    
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(1,2) = vertices(1,pointer_reference_to_local(2))
	vertices_triangle(1,3) = Dx
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(3))
	vertices_triangle(2,2) = vertices(2,pointer_reference_to_local(2))
	vertices_triangle(2,3) = Dy
 
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
 
    uh = 0.0
    u_true = 0.0
    piece_flag = 2
    DO i = 1, Gpn
 
        DO j = 1, 4
                
		    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                    information_vector_1,information_vector_2,piece_flag,test_basis_type,j, &
										    test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 
            uh(i) = uh(i) + U_fe(HT(j,elem_index)) * r2
            
        END DO
            
        IF (beta2 == Global_Beta(1)) THEN
            coefficient_function_name = Global_Beta(1)
        ELSE
            coefficient_function_name = Global_Beta(2)
        ENDIF
            
        CALL Function_True(delta,coefficient_function_name, Gauss_point_local_triangle_Nine(i,1),&
                            Gauss_point_local_triangle_Nine(i,2), derivative_degree_x, derivative_degree_y,temp2)
        u_true(i) = temp2
            
        !error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= \\
        IF (delta == 0) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))
        ELSEIF (delta == 1) THEN
            error = error + Gauss_coefficient_local_triangle_Nine(i) * (uh(i)-u_true(i)) * (uh(i)-u_true(i))* Gauss_point_local_triangle_Nine(i,2)    
        ENDIF 
        !$ ============== ab.ZWZ for adding cylindrical coordinate ============= //
    
    ENDDO
    
ENDIF


END