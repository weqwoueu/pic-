SUBROUTINE  Gauss_integration_local_stiffness_IFE(delta, Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
			                                      vertices,information_vector_1,information_vector_2,trial_basis_type, &
												  trial_basis_index,trial_derivative_degree_x,trial_derivative_degree_y,test_basis_type, &
												  test_basis_index,test_derivative_degree_x,test_derivative_degree_y,r, element_Gauss)
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Generate_Gauss_local_triangle, Retangular_local_basis_IFE		
IMPLICIT NONE

!$ 2021/4/19 zwz: 改了判断delta的位置，减少代码量, 并且为了impic加了coefficient函数

REAL(8),DIMENSION(:)				        ::    Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:)		                ::    Gauss_point_reference_triangle

INTEGER                                           information_vector_1(18)
REAL(8)                                           information_vector_2(8)
REAL(8)                                           vertices(2,4)
INTEGER                                           trial_basis_type, trial_derivative_degree_x,trial_derivative_degree_y, &
												  test_basis_type,test_derivative_degree_x,test_derivative_degree_y
INTEGER                                           trial_basis_index
INTEGER                                           test_basis_index
REAL(8)								        	  r, r1, r2
INTEGER                                           Gpn,i
INTEGER                                           interface_element_type
INTEGER                                           pointer_reference_to_local(4)
REAL(8)	                                          beta1, beta2, Dx, Dy, Ex, Ey
REAL(8)                                           vertices_triangle(2,3)
INTEGER                                     ::	  delta
REAL(8)                                           coefficient_function_name_impic !!!bjw add for impic 2019-6-3

!==========LY modification, 2022-7-25==========
Integer :: element_Gauss
!==========LY modification, 2022-7-25==========

Gpn = SIZE(Gauss_coefficient_reference_triangle) 
interface_element_type=information_vector_1(6)
pointer_reference_to_local=information_vector_1(11:14)
beta1=information_vector_2(1)
beta2=information_vector_2(2)
Dx=information_vector_2(3)
Dy=information_vector_2(4)
Ex=information_vector_2(5)
Ey=information_vector_2(6)


r=0;

IF (interface_element_type==1) THEN
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(1))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(1))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey
        
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
        
    DO i = 1, Gpn
        
	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,1,trial_basis_type,trial_basis_index,&
										trial_derivative_degree_x,trial_derivative_degree_y,r1, element_Gauss) 

		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,1,test_basis_type,test_basis_index, &
										test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 

        !!! ************************ bjw add for impic 2019-6-3 **********************************************
        CALL Function_coefficient_2D(coefficient_function_name_impic,beta1, &
                                        Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2), &
                                        trial_derivative_degree_x,trial_derivative_degree_y, &
                                        test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
        !!! ************************ bjw add for impic 2019-6-3 **********************************************
            
        IF (delta == 0) THEN

            r=r+Gauss_coefficient_local_triangle_Nine(i) * coefficient_function_name_impic * r1*r2
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta1*r1*r2

        ELSEIF (delta ==1) THEN

            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta1*r1*r2*Gauss_point_local_triangle_Nine(i,2)
            r=r+Gauss_coefficient_local_triangle_Nine(i)*coefficient_function_name_impic*r1*r2*Gauss_point_local_triangle_Nine(i,2)

        ELSE 
		      PRINT*, ' delta value is wrong, check again, stop'
		      STOP
        ENDIF
	
    END DO
    
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(4))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = vertices(1,pointer_reference_to_local(3))
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(4))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = vertices(2,pointer_reference_to_local(3))

    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
    DO i = 1, Gpn
 
	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,trial_basis_type,trial_basis_index,&
										trial_derivative_degree_x,trial_derivative_degree_y,r1, element_Gauss) 

		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
										test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 

        !!! ************************ bjw add for impic 2019-6-3 **********************************************
        CALL Function_coefficient_2D(coefficient_function_name_impic,beta2, &
                                        Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2), &
                                        trial_derivative_degree_x,trial_derivative_degree_y, &
                                        test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
        !!! ************************ bjw add for impic 2019-6-3 **********************************************
            
        IF (delta == 0) THEN

            r=r+Gauss_coefficient_local_triangle_Nine(i) * coefficient_function_name_impic * r1*r2
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2

        ELSEIF (delta ==1) THEN

            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2*Gauss_point_local_triangle_Nine(i,2)
            r=r+Gauss_coefficient_local_triangle_Nine(i)*coefficient_function_name_impic*r1*r2*Gauss_point_local_triangle_Nine(i,2)

        ELSE 
		      PRINT*, ' delta value is wrong, check again, stop'
		      STOP
        END IF
 
    END DO
    
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(3))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
    DO i = 1, Gpn

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,trial_basis_type,trial_basis_index,&
										trial_derivative_degree_x,trial_derivative_degree_y,r1, element_Gauss) 

		CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
										test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 

        !!! ************************ bjw add for impic 2019-6-3 **********************************************
        CALL Function_coefficient_2D(coefficient_function_name_impic,beta2, &
                                        Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2), &
                                        trial_derivative_degree_x,trial_derivative_degree_y, &
                                        test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
        !!! ************************ bjw add for impic 2019-6-3 **********************************************
            
        IF (delta == 0) THEN
            r=r+Gauss_coefficient_local_triangle_Nine(i) * coefficient_function_name_impic * r1*r2
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2

        ELSEIF (delta ==1) THEN

            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2*Gauss_point_local_triangle_Nine(i,2)
            r=r+Gauss_coefficient_local_triangle_Nine(i)*coefficient_function_name_impic*r1*r2*Gauss_point_local_triangle_Nine(i,2)

        ELSE 
		      PRINT*, ' delta value is wrong, check again, stop'
		      STOP
        END IF
    END DO
    
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(2))
	vertices_triangle(1,2) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(2))
	vertices_triangle(2,2) = vertices(2,pointer_reference_to_local(3))
	vertices_triangle(2,3) = Ey
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine)  
    DO i = 1, Gpn

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,trial_basis_type,trial_basis_index,&
									    trial_derivative_degree_x,trial_derivative_degree_y,r1, element_Gauss) 

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
									    test_derivative_degree_x, test_derivative_degree_y,r2, element_Gauss) 

        !!! ************************ bjw add for impic 2019-6-3 **********************************************
        CALL Function_coefficient_2D(coefficient_function_name_impic,beta2, &
                                        Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2), &
                                        trial_derivative_degree_x,trial_derivative_degree_y, &
                                        test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
        !!! ************************ bjw add for impic 2019-6-3 **********************************************
            
        IF (delta == 0) THEN
            r=r+Gauss_coefficient_local_triangle_Nine(i) * coefficient_function_name_impic * r1*r2
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2

        ELSEIF (delta ==1) THEN
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2*Gauss_point_local_triangle_Nine(i,2)
            r=r+Gauss_coefficient_local_triangle_Nine(i)*coefficient_function_name_impic*r1*r2*Gauss_point_local_triangle_Nine(i,2)
        ELSE 
		        PRINT*, ' delta value is wrong, check again, stop'
		        STOP
        END IF

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
    DO i = 1, Gpn

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,1,trial_basis_type,trial_basis_index, &
									    trial_derivative_degree_x,trial_derivative_degree_y,r1, element_Gauss)

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,1,test_basis_type,test_basis_index, &
									    test_derivative_degree_x,test_derivative_degree_y,r2, element_Gauss) 

        !!! ************************ bjw add for impic 2019-6-3 **********************************************
        CALL Function_coefficient_2D(coefficient_function_name_impic,beta1, &
                                        Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2), &
                                        trial_derivative_degree_x,trial_derivative_degree_y, &
                                        test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
        !!! ************************ bjw add for impic 2019-6-3 **********************************************
    
        IF (delta == 0) THEN
            r=r+Gauss_coefficient_local_triangle_Nine(i) * coefficient_function_name_impic * r1*r2
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta1*r1*r2
        ELSEIF (delta ==1) THEN
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta1*r1*r2*Gauss_point_local_triangle_Nine(i,2)
            r=r+Gauss_coefficient_local_triangle_Nine(i)*coefficient_function_name_impic*r1*r2*Gauss_point_local_triangle_Nine(i,2)
        ELSE 
		    PRINT*, ' delta value is wrong, check again, stop'
		    STOP
        ENDIF  
    ENDDO
    
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(1))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(1))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
    DO i = 1, Gpn

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,1,trial_basis_type,trial_basis_index, &
									    trial_derivative_degree_x,trial_derivative_degree_y,r1, element_Gauss)

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,1,test_basis_type,test_basis_index, &
									    test_derivative_degree_x,test_derivative_degree_y,r2, element_Gauss) 

        !!! ************************ bjw add for impic 2019-6-3 **********************************************
        CALL Function_coefficient_2D(coefficient_function_name_impic,beta1, &
                                        Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2), &
                                        trial_derivative_degree_x,trial_derivative_degree_y, &
                                        test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
        !!! ************************ bjw add for impic 2019-6-3 **********************************************
        
        IF (delta == 0) THEN
            r=r+Gauss_coefficient_local_triangle_Nine(i) * coefficient_function_name_impic * r1*r2
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta1*r1*r2
        ELSEIF (delta ==1) THEN
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta1*r1*r2*Gauss_point_local_triangle_Nine(i,2)
            r=r+Gauss_coefficient_local_triangle_Nine(i)*coefficient_function_name_impic*r1*r2*Gauss_point_local_triangle_Nine(i,2)
        ELSE 
		    PRINT*, ' delta value is wrong, check again, stop'
		    STOP
        ENDIF 

    ENDDO
    
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(2))
	vertices_triangle(1,2) = Dx
	vertices_triangle(1,3) = Ex
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(2))
	vertices_triangle(2,2) = Dy
	vertices_triangle(2,3) = Ey
    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 
    DO i = 1, Gpn

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,trial_basis_type,trial_basis_index, &
									    trial_derivative_degree_x,trial_derivative_degree_y,r1, element_Gauss)

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
									    test_derivative_degree_x,test_derivative_degree_y,r2, element_Gauss) 

        !!! ************************ bjw add for impic 2019-6-3 **********************************************
        CALL Function_coefficient_2D(coefficient_function_name_impic,beta2, &
                                        Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2), &
                                        trial_derivative_degree_x,trial_derivative_degree_y, &
                                        test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
        !!! ************************ bjw add for impic 2019-6-3 **********************************************

        IF (delta == 0) THEN
            r=r+Gauss_coefficient_local_triangle_Nine(i) * coefficient_function_name_impic * r1*r2
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2
        ELSEIF (delta ==1) THEN
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2*Gauss_point_local_triangle_Nine(i,2)
            r=r+Gauss_coefficient_local_triangle_Nine(i)*coefficient_function_name_impic*r1*r2*Gauss_point_local_triangle_Nine(i,2)
        ELSE 
		    PRINT*, ' delta value is wrong, check again, stop'
		    STOP
        ENDIF 

    ENDDO
    
    vertices_triangle(1,1) = vertices(1,pointer_reference_to_local(3))
	vertices_triangle(1,2) = vertices(1,pointer_reference_to_local(2))
	vertices_triangle(1,3) = Dx
    vertices_triangle(2,1) = vertices(2,pointer_reference_to_local(3))
	vertices_triangle(2,2) = vertices(2,pointer_reference_to_local(2))
	vertices_triangle(2,3) = Dy

    CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
	                                   vertices_triangle,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine) 

    DO i = 1, Gpn

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,trial_basis_type,trial_basis_index, &
									    trial_derivative_degree_x,trial_derivative_degree_y,r1, element_Gauss)

	    CALL Retangular_local_basis_IFE(Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2),vertices, &
		                                information_vector_1,information_vector_2,2,test_basis_type,test_basis_index, &
									    test_derivative_degree_x,test_derivative_degree_y,r2, element_Gauss) 

        !!! ************************ bjw add for impic 2019-6-3 **********************************************
        CALL Function_coefficient_2D(coefficient_function_name_impic,beta2, &
                                        Gauss_point_local_triangle_Nine(i,1),Gauss_point_local_triangle_Nine(i,2), &
                                        trial_derivative_degree_x,trial_derivative_degree_y, &
                                        test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
        !!! ************************ bjw add for impic 2019-6-3 **********************************************

        IF (delta == 0) THEN
            r=r+Gauss_coefficient_local_triangle_Nine(i) * coefficient_function_name_impic * r1*r2
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2
        ELSEIF (delta ==1) THEN
            !r=r+Gauss_coefficient_local_triangle_Nine(i)*beta2*r1*r2*Gauss_point_local_triangle_Nine(i,2)
            r=r+Gauss_coefficient_local_triangle_Nine(i)*coefficient_function_name_impic*r1*r2*Gauss_point_local_triangle_Nine(i,2)
        ELSE 
		    PRINT*, ' delta value is wrong, check again, stop'
		    STOP
        ENDIF 

    ENDDO
    
ENDIF

!print*,'r=',r

END