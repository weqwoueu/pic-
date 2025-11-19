SUBROUTINE  Retangular_local_basis_IFE(x,y,vertices,information_vector_1,information_vector_2,piece_flag,basis_type, &
		                               basis_index,derivative_degree_x,derivative_degree_y,r,element_Gauss) 

IMPLICIT NONE

REAL(8)                                           x, y
INTEGER                                           information_vector_1(18)
REAL(8)                                           information_vector_2(8)
REAL(8)                                           vertices(2,4)
INTEGER                                           basis_type, basis_index, derivative_degree_x, derivative_degree_y  
REAL(8)								        	  r, r1, r2
INTEGER                                           Gpn
INTEGER                                           interface_element_type
INTEGER                                           pointer_reference_to_local(4), pointer_local_to_reference(4)
REAL(8)	                                          beta1, beta2, Dx, Dy, Ex, Ey, x1, x2, x4, y1, y2, y4, bottom, kesai, eita, a, b
REAL(8)                                           vertices_triangle(2,3)
INTEGER                                           piece_flag

!=========LY modification, 2022-7-25=========
Integer :: element_Gauss
!=========LY modification, 2022-7-25=========

IF (basis_type == 1) THEN
    
    interface_element_type = information_vector_1(6) 
    pointer_local_to_reference = information_vector_1(7:10) 
    pointer_reference_to_local = information_vector_1(11:14) 
    beta1 = information_vector_2(1) 
    beta2 = information_vector_2(2) 
    Dx = information_vector_2(3) 
    Dy = information_vector_2(4) 
    Ex = information_vector_2(5) 
    Ey = information_vector_2(6)   
    
    
    x1 = vertices(1,pointer_reference_to_local(1)) 
    y1 = vertices(2,pointer_reference_to_local(1)) 
    x2 = vertices(1,pointer_reference_to_local(2)) 
    y2 = vertices(2,pointer_reference_to_local(2))     
    x4 = vertices(1,pointer_reference_to_local(4)) 
    y4 = vertices(2,pointer_reference_to_local(4)) 

    bottom = (x2-x1)*(y4-y1)-(x4-x1)*(y2-y1) 
    kesai = ((y4-y1)*(x-x1)-(x4-x1)*(y-y1))/bottom 
    eita = ((x2-x1)*(y-y1)-(y2-y1)*(x-x1))/bottom 


!in the real coordinate position


    IF (interface_element_type == 1) THEN

        a = ((y4-y1)*(Ex-x1)-(x4-x1)*(Ey-y1))/bottom 
        b = ((x2-x1)*(Dy-y1)-(y2-y1)*(Dx-x1))/bottom 
		

    ELSEIF (interface_element_type == 2) THEN

        a = ((y4-y1)*(Ex-x1)-(x4-x1)*(Ey-y1))/bottom 
        b = ((y4-y1)*(Dx-x1)-(x4-x1)*(Dy-y1))/bottom 
    
    ENDIF
    
    !$ =-------------------要传给function_coefficient affine mapping 之前的x,y值 -------------------------------------
    
    IF (derivative_degree_x == 0 .AND. derivative_degree_y == 0) THEN

         CALL Retangular_reference_basis_IFE(x,y,kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
		                                     pointer_local_to_reference(basis_index),0,0,r,element_Gauss) 
   
    ELSEIF (derivative_degree_x == 1 .AND. derivative_degree_y == 0) THEN
    
	     CALL Retangular_reference_basis_IFE(x,y,kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
		                                     pointer_local_to_reference(basis_index),1,0,r1,element_Gauss) 

	   	 CALL Retangular_reference_basis_IFE(x,y,kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
		                                     pointer_local_to_reference(basis_index),0,1,r2,element_Gauss)
			
		 r = r1*(y4-y1)/bottom+r2*(y1-y2)/bottom 
    
	ELSEIF (derivative_degree_x == 0 .AND. derivative_degree_y == 1) THEN
   
         CALL Retangular_reference_basis_IFE(x,y,kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
		                                     pointer_local_to_reference(basis_index),1,0,r1,element_Gauss)
!!!										 *(x1-x4)/bottom+
		 CALL Retangular_reference_basis_IFE(x,y,kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
		                                     pointer_local_to_reference(basis_index),0,1,r2,element_Gauss)
!!!  									 *(x2-x1)/bottom 
         
		 r = r1*(x1-x4)/bottom+r2*(x2-x1)/bottom 

     
	ELSEIF (derivative_degree_x == 1 .AND. derivative_degree_y == 1) THEN
   
         CALL Retangular_reference_basis_IFE(x,y,kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
		                                     pointer_local_to_reference(basis_index),1,1,r1,element_Gauss)
	 										 
         CALL Retangular_reference_basis_IFE(x,y,kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
		                                     pointer_local_to_reference(basis_index),1,1,r2,element_Gauss)
											 
		r = r1*(x2-x1)*(y4-y1)/(bottom**2)+r2*(x1-x4)*(y1-y2)/(bottom**2) 
   
    ENDIF

!    IF (derivative_degree_x == 0 .AND. derivative_degree_y == 0) THEN
!
!         CALL Retangular_reference_basis_IFE(kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
!		                                     pointer_local_to_reference(basis_index),0,0,r) 
!   
!    ELSEIF (derivative_degree_x == 1 .AND. derivative_degree_y == 0) THEN
!    
!	     CALL Retangular_reference_basis_IFE(kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
!		                                     pointer_local_to_reference(basis_index),1,0,r1) 
!
!	   	 CALL Retangular_reference_basis_IFE(kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
!		                                     pointer_local_to_reference(basis_index),0,1,r2)
!			
!		 r = r1*(y4-y1)/bottom+r2*(y1-y2)/bottom 
!    
!	ELSEIF (derivative_degree_x == 0 .AND. derivative_degree_y == 1) THEN
!   
!         CALL Retangular_reference_basis_IFE(kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
!		                                     pointer_local_to_reference(basis_index),1,0,r1)
!!!!										 *(x1-x4)/bottom+
!		 CALL Retangular_reference_basis_IFE(kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
!		                                     pointer_local_to_reference(basis_index),0,1,r2)
!!!!  									 *(x2-x1)/bottom 
!         
!		 r = r1*(x1-x4)/bottom+r2*(x2-x1)/bottom 
!
!     
!	ELSEIF (derivative_degree_x == 1 .AND. derivative_degree_y == 1) THEN
!   
!         CALL Retangular_reference_basis_IFE(kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
!		                                     pointer_local_to_reference(basis_index),1,1,r1)
!	 										 
!         CALL Retangular_reference_basis_IFE(kesai,eita,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type, &
!		                                     pointer_local_to_reference(basis_index),1,1,r2)
!											 
!		r = r1*(x2-x1)*(y4-y1)/(bottom**2)+r2*(x1-x4)*(y1-y2)/(bottom**2) 
!   
!    ENDIF
    
ENDIF

END