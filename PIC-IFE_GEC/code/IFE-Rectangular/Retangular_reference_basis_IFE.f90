!SUBROUTINE  Retangular_reference_basis_IFE(x,y,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type,basis_index,  &
!                                           derivative_degree_x,derivative_degree_y, r1)
SUBROUTINE  Retangular_reference_basis_IFE(x_beforeaffine,y_beforeaffine, &
                                           x,y,interface_element_type,piece_flag,a,b,beta1,beta2,basis_type,basis_index,  &
                                           derivative_degree_x,derivative_degree_y, r1, element_Gauss)


IMPLICIT NONE
REAL(8)     x_beforeaffine, y_beforeaffine

REAL(8)                                           x, y

INTEGER                                           basis_type, basis_index, derivative_degree_x, derivative_degree_y  
REAL(8)								        	  r1
INTEGER                                           Gpn
INTEGER                                           interface_element_type

REAL(8)	                                          beta1, beta2, bottom, kesai, eita, a, b, R, C1, C2, C3, C4 !, Dx, Dy, Ex, Ey, x1, x2, x4, y1, y2, y4, 

INTEGER                                           piece_flag

REAL(8)    k11, k12, k21, k22

!=========LY modification, 2022-7-25=========
Integer :: element_Gauss
!=========LY modification, 2022-7-25=========

!CALL Function_coefficient_2D(k11,0,0,0, 1,0,1,0)
!CALL Function_coefficient_2D(k12,0,0,0, 1,0,0,1)
!CALL Function_coefficient_2D(k21,0,0,0, 0,1,1,0)
!CALL Function_coefficient_2D(k22,0,0,0, 0,1,0,1)
!CALL Function_coefficient_2D(k11,0._8,0._8,0._8, 1,0,1,0)
!CALL Function_coefficient_2D(k12,0._8,0._8,0._8, 1,0,0,1)
!CALL Function_coefficient_2D(k21,0._8,0._8,0._8, 0,1,1,0)
!CALL Function_coefficient_2D(k22,0._8,0._8,0._8, 0,1,0,1)

!=========LY modification, 2022-7-25=========
!=========PIC=========
CALL Function_coefficient_2D(k11,1._8, x_beforeaffine, y_beforeaffine, 1,0,1,0, element_Gauss)
CALL Function_coefficient_2D(k12,0._8, x_beforeaffine, y_beforeaffine, 0,1,1,0, element_Gauss)
CALL Function_coefficient_2D(k21,0._8, x_beforeaffine, y_beforeaffine, 1,0,0,1, element_Gauss)
CALL Function_coefficient_2D(k22,1._8, x_beforeaffine, y_beforeaffine, 0,1,0,1, element_Gauss)
!=========PIC=========
!=========LY modification for checking convergence=========
!=========IFE=========
!CALL Function_coefficient_2D(k11,1._8, x_beforeaffine, y_beforeaffine, 1,0,1,0, element_Gauss)
!CALL Function_coefficient_2D(k12,1._8, x_beforeaffine, y_beforeaffine, 0,1,1,0, element_Gauss)
!CALL Function_coefficient_2D(k21,1._8, x_beforeaffine, y_beforeaffine, 1,0,0,1, element_Gauss)
!CALL Function_coefficient_2D(k22,1._8, x_beforeaffine, y_beforeaffine, 0,1,0,1, element_Gauss)
!=========IFE=========
!=========LY modification for checking convergence=========
!=========LY modification, 2022-7-25=========

!$ ============= mb.ZWZ 2021/5/10 ======================== \\
!CALL Function_coefficient_2D(k11,1._8, x_beforeaffine, y_beforeaffine, 1,0,1,0)
!CALL Function_coefficient_2D(k12,0._8, x_beforeaffine, y_beforeaffine, 1,0,0,1)
!CALL Function_coefficient_2D(k21,0._8, x_beforeaffine, y_beforeaffine, 0,1,1,0)
!CALL Function_coefficient_2D(k22,1._8, x_beforeaffine, y_beforeaffine, 0,1,0,1)
!$ ====== ab.ZWZ when checking convergence ============== \\
!CALL Function_coefficient_2D(k11,1._8,0._8,0._8, 1,0,1,0)
!CALL Function_coefficient_2D(k12,1._8,0._8,0._8, 1,0,0,1)
!CALL Function_coefficient_2D(k21,1._8,0._8,0._8, 0,1,1,0)
!CALL Function_coefficient_2D(k22,1._8,0._8,0._8, 0,1,0,1)
!$ ====== ab.ZWZ when checking convergence ============== //

IF (basis_type == 1) THEN

    R = beta1/beta2 

    IF (interface_element_type ==1) THEN  !!邻边
        !!! 参考何老师博士论文P21-P27... bottom = W, 
        !!! bjw 各项异性基函数推导
                
        bottom = 2*R*(a**2)*k22 - (a**2)*(b**2)*k21 - (a**2)*(b**2)*k12 + 2*R*(b**2)*k11 + 2*a*(b**2)*k11 + 2*a*(b**2)*k12  &
                 - a*(b**3)*k11 + 2*(a**2)*b*k21 + 2*(a**2)*b*k22 - (a**3)*b*k22 - 2*R*a*(b**2)*k11 - 2*R*a*(b**2)*k12  &
                 + R*a*(b**3)*k11 - 2*R*(a**2)*b*k21 - 2*R*(a**2)*b*k22 + R*(a**3)*b*k22 + R*(a**2)*(b**2)*k12  &
                 + R*(a**2)*(b**2)*k21 + 2*R*a*b*k12 + 2*R*a*b*k21
                         
        IF (basis_index ==1) THEN
            
            IF (piece_flag ==1) THEN
                C1 = 1 
                C2 =  ((b**3)*k11 - 2*(b**2)*k12 - 2*(b**2)*k11 - 2*a*b*k21 - 2*a*b*k22 - 2*R*(a**2)*k22 + 2*R*(b**2)*k12  &
                     - R*(b**3)*k11 + a*(b**2)*k12 + a*(b**2)*k21 + (a**2)*b*k22 - R*a*(b**2)*k12 - R*a*(b**2)*k21 -  &
                     R*(a**2)*b*k22 - 2*R*a*b*k12 + 2*R*a*b*k22)/bottom 
                C3 = ((a**3)*k22 - 2*(a**2)*k22 - 2*(a**2)*k21 - 2*a*b*k11 - 2*a*b*k12 + 2*R*(a**2)*k21 - R*(a**3)*k22 - &
                     2*R*(b**2)*k11 + a*(b**2)*k11 + (a**2)*b*k12 + (a**2)*b*k21 - R*a*(b**2)*k11 - R*(a**2)*b*k12 - &
                     R*(a**2)*b*k21 + 2*R*a*b*k11 - 2*R*a*b*k21)/bottom 
                C4 = 2*R*((a**2)*k22 + (b**2)*k11 + a*b*k12 + a*b*k21)/bottom
             
            ELSEIF (piece_flag ==2) THEN
                C1 = 2*R*((a**2)*k22 + (b**2)*k11 + a*b*k12 + a*b*k21)/bottom
                C2 = -C1 
                C3 = -C1 
                C4 = C1
              
            ENDIF
            
        ELSEIF (basis_index ==2) THEN
                             
            IF (piece_flag ==1) THEN
                C1 = 0 
                C2 = (2*(b**2)*k11 - (b**3)*k11 + 2*a*b*k21 + 2*R*(a**2)*k22 + R*(b**3)*k11 + a*(b**2)*k12 - a*(b**2)*k21 + &
                    (a**2)*b*k22 - R*a*(b**2)*k12 + R*a*(b**2)*k21 - R*(a**2)*b*k22 + 2*R*a*b*k12)/bottom 
                C3 = (2*(a**2)*k21 - (a**3)*k22 + (a**2)*(b**2)*k12 + (a**2)*(b**2)*k21 + 2*a*b*k11 - 2*R*(a**2)*k21 + &
                     R*(a**3)*k22 - 3*a*(b**2)*k11 + a*(b**3)*k11 - (a**2)*b*k12 - 3*(a**2)*b*k21 + (a**3)*b*k22 + &
                     3*R*a*(b**2)*k11 - R*a*(b**3)*k11 + R*(a**2)*b*k12 + 3*R*(a**2)*b*k21 - R*(a**3)*b*k22 - &
                     R*(a**2)*(b**2)*k12 - R*(a**2)*(b**2)*k21 - 2*R*a*b*k11)/bottom 
                C4 = (2*R*a*(b**2)*k12 - 2*R*(b**2)*k11 - 2*a*(b**2)*k12 - 2*(a**2)*b*k22 - 2*R*(a**2)*k22 + &
                     2*R*(a**2)*b*k22 - 2*R*a*b*k12 - 2*R*a*b*k21)/bottom
        
            ELSEIF (piece_flag ==2) THEN
                C1 = (2*a*(b**2)*k11 - (a**2)*(b**2)*k21 - (a**2)*(b**2)*k12 - a*(b**3)*k11 + 2*(a**2)*b*k21 - (a**3)*b*k22 - &
                     2*R*a*(b**2)*k11 + R*a*(b**3)*k11 - 2*R*(a**2)*b*k21 + R*(a**3)*b*k22 + R*(a**2)*(b**2)*k12 + &
                     R*(a**2)*(b**2)*k21)/bottom 
                C2 = 1-C1 
                C3 = -C1
                C4 = C1-1
               
            ENDIF
            
        ELSEIF (basis_index ==3) THEN
 !           piece_flag=2
            IF (piece_flag ==1) THEN
                C1 = 0 
                C2 = ((b**3)*k11 - (a**2)*(b**2)*k12 - (a**2)*(b**2)*k21 - R*(b**3)*k11 + a*(b**2)*k12 - a*(b**3)*k11 + &
                     a*(b**2)*k21 + (a**2)*b*k22 - (a**3)*b*k22 - R*a*(b**2)*k12 + R*a*(b**3)*k11 - R*a*(b**2)*k21 - &
                     R*(a**2)*b*k22 + R*(a**3)*b*k22 + R*(a**2)*(b**2)*k12 + R*(a**2)*(b**2)*k21)/bottom 
                C3 = ((a**3)*k22 - (a**2)*(b**2)*k12 - (a**2)*(b**2)*k21 - R*(a**3)*k22 + a*(b**2)*k11 - a*(b**3)*k11 + &
                     (a**2)*b*k12 + (a**2)*b*k21 - (a**3)*b*k22 - R*a*(b**2)*k11 + R*a*(b**3)*k11 - R*(a**2)*b*k12 - &
                     R*(a**2)*b*k21 + R*(a**3)*b*k22 + R*(a**2)*(b**2)*k12 + R*(a**2)*(b**2)*k21)/bottom 
                C4 = (2*R*(a**2)*k22 + 2*R*(b**2)*k11 + 2*a*(b**2)*k11 + 2*a*(b**2)*k12 + 2*(a**2)*b*k21 + 2*(a**2)*b*k22 - &
                     2*R*a*(b**2)*k11 - 2*R*a*(b**2)*k12 - 2*R*(a**2)*b*k21 - 2*R*(a**2)*b*k22 + 2*R*a*b*k12 + 2*R*a*b*k21)/bottom
   
            ELSEIF (piece_flag ==2) THEN
                C1 = ((a**2)*(b**2)*k12 + (a**2)*(b**2)*k21 + a*(b**3)*k11 + (a**3)*b*k22 - R*a*(b**3)*k11 - R*(a**3)*b*k22 - &
                     R*(a**2)*(b**2)*k12 - R*(a**2)*(b**2)*k21)/bottom 
                C2 = -C1 
                C3 = -C1 
                C4 = 1+C1
               
            ENDIF
            
        ELSEIF (basis_index ==4) THEN
            
            IF (piece_flag ==1) THEN
                C1 = 0 
                C2 = (2*(b**2)*k12 - (b**3)*k11 + (a**2)*(b**2)*k12 + (a**2)*(b**2)*k21 + 2*a*b*k22 - 2*R*(b**2)*k12 + &
                     R*(b**3)*k11 - 3*a*(b**2)*k12 + a*(b**3)*k11 - a*(b**2)*k21 - 3*(a**2)*b*k22 + (a**3)*b*k22 + &
                     3*R*a*(b**2)*k12 - R*a*(b**3)*k11 + R*a*(b**2)*k21 + 3*R*(a**2)*b*k22 - R*(a**3)*b*k22 - &
                     R*(a**2)*(b**2)*k12 - R*(a**2)*(b**2)*k21 - 2*R*a*b*k22)/bottom 
                C3 = (2*(a**2)*k22 - (a**3)*k22 + 2*a*b*k12 + R*(a**3)*k22 + 2*R*(b**2)*k11 + a*(b**2)*k11 - &
                     (a**2)*b*k12 + (a**2)*b*k21 - R*a*(b**2)*k11 + R*(a**2)*b*k12 - R*(a**2)*b*k21 + 2*R*a*b*k21)/bottom 
                C4 = (2*R*a*(b**2)*k11 - 2*R*(b**2)*k11 - 2*a*(b**2)*k11 - 2*(a**2)*b*k21 - 2*R*(a**2)*k22 + &
                      2*R*(a**2)*b*k21 - 2*R*a*b*k12 - 2*R*a*b*k21)/bottom

            ELSEIF (piece_flag ==2) THEN
                C1 = (2*a*(b**2)*k12 - (a**2)*(b**2)*k21 - (a**2)*(b**2)*k12 - a*(b**3)*k11 + 2*(a**2)*b*k22 - &
                     (a**3)*b*k22 - 2*R*a*(b**2)*k12 + R*a*(b**3)*k11 - 2*R*(a**2)*b*k22 + R*(a**3)*b*k22 + &
                      R*(a**2)*(b**2)*k12 + R*(a**2)*(b**2)*k21)/bottom 
                C2 = -C1 
                C3 = 1-C1 
                C4 = C1-1
              
            ENDIF
            
        ENDIF     
        
    ELSEIF (interface_element_type ==2) THEN !!对边
        
        bottom = 2*R*k11 + a*k11 + 2*a*k12 + b*k11 - 2*b*k12 - (a**2)*k12 + (a**2)*k21 + 2*(a**2)*k22 - (a**3)*k22 + &
                (b**2)*k12 - (b**2)*k21 + 2*(b**2)*k22 - (b**3)*k22 - R*a*k11 + 2*R*a*k21 - R*b*k11 - 2*R*b*k21 - &
                 4*a*b*k22 + R*(a**2)*k12 - R*(a**2)*k21 + R*(a**3)*k22 - R*(b**2)*k12 + R*(b**2)*k21 + R*(b**3)*k22 + &
                 a*(b**2)*k22 + (a**2)*b*k22 - R*a*(b**2)*k22 - R*(a**2)*b*k22
        IF (basis_index ==1) THEN
            
            IF (piece_flag ==1) THEN
                C1 = 1 
                C2 = (2*R*k12 - 2*k12 - R*k11 - k11 + a*k12 - a*k21 - 2*a*k22 - b*k11 + 3*b*k12 + b*k21 + 2*b*k22 + &
                     (a**2)*k22 - (b**2)*k12 + (b**2)*k21 - 3*(b**2)*k22 + (b**3)*k22 - 3*R*a*k12 - R*a*k21 + &
                      2*R*a*k22 + R*b*k11 - R*b*k12 + R*b*k21 - 2*R*b*k22 - a*b*k12 - a*b*k21 + 2*a*b*k22 - &
                      3*R*(a**2)*k22 + R*(b**2)*k12 - R*(b**2)*k21 + R*(b**2)*k22 - R*(b**3)*k22 - (a**2)*b*k22 + &
                      R*(a**2)*b*k22 + R*a*b*k12 + R*a*b*k21 + 2*R*a*b*k22)/bottom 
                C3 = -1 
                C4 = (2*R*k11 + 2*b*k11 - 2*(b**2)*k21 + 2*R*a*k12 + 2*R*a*k21 - 2*R*b*k11 - 2*R*b*k12 - 2*R*b*k21 + &
                      2*a*b*k21 + 2*R*(a**2)*k22 + 2*R*(b**2)*k21 + 2*R*(b**2)*k22 - 2*R*a*b*k21 - 4*R*a*b*k22)/bottom

            ELSEIF (piece_flag ==2) THEN
                C1 = (2*R*k11 + b*k11 - 2*b*k12 + (b**2)*k12 - (b**2)*k21 + 2*(b**2)*k22 - (b**3)*k22 + 2*R*a*k12 + &
                      2*R*a*k21 - R*b*k11 - 2*R*b*k21 + a*b*k12 + a*b*k21 - 2*a*b*k22 + 2*R*(a**2)*k22 - R*(b**2)*k12 + &
                      R*(b**2)*k21 + R*(b**3)*k22 + (a**2)*b*k22 - R*(a**2)*b*k22 - R*a*b*k12 - R*a*b*k21 - 2*R*a*b*k22)/bottom 
                C2 = -C1 
                C3 = (2*(b**2)*k21 - 2*b*k11 - 2*R*k11 - 2*R*a*k12 - 2*R*a*k21 + 2*R*b*k11 + 2*R*b*k12 + 2*R*b*k21 - &
                      2*a*b*k21 - 2*R*(a**2)*k22 - 2*R*(b**2)*k21 - 2*R*(b**2)*k22 + 2*R*a*b*k21 + 4*R*a*b*k22)/bottom 
                C4 = -C3
                
            ENDIF
            
        ELSEIF (basis_index ==2) THEN
            
            IF (piece_flag ==1) THEN
                C1 = 0 
                C2 = (k11 + R*k11 + a*k12 + a*k21 + b*k11 - 3*b*k12 - b*k21 + (a**2)*k22 + (b**2)*k12 - (b**2)*k21 + &
                      3*(b**2)*k22 - (b**3)*k22 + R*a*k12 + R*a*k21 - R*b*k11 + R*b*k12 - R*b*k21 + a*b*k12 + a*b*k21 - &
                      4*a*b*k22 + R*(a**2)*k22 - R*(b**2)*k12 + R*(b**2)*k21 - R*(b**2)*k22 + R*(b**3)*k22 + (a**2)*b*k22 - &
                      R*(a**2)*b*k22 - R*a*b*k12 - R*a*b*k21)/bottom 
                C3 = 0 
                C4 = (2*b*k12 - 2*a*k12 - 2*b*k11 - 2*R*k11 - 2*(a**2)*k22 + 2*(b**2)*k21 - 2*(b**2)*k22 - 2*R*a*k21 + &
                      2*R*b*k11 + 2*R*b*k21 - 2*a*b*k21 + 4*a*b*k22 - 2*R*(b**2)*k21 + 2*R*a*b*k21)/bottom

            ELSEIF (piece_flag ==2) THEN
                C1 = (a*k11 - (a**2)*k12 + (a**2)*k21 - (a**3)*k22 - R*a*k11 - a*b*k12 - a*b*k21 + R*(a**2)*k12 - &
                      R*(a**2)*k21 + R*(a**3)*k22 + a*(b**2)*k22 - R*a*(b**2)*k22 + R*a*b*k12 + R*a*b*k21)/bottom 
                C2 = 1-C1 
                C3 = (b*k11 - a*k11 + (a**2)*k12 - (a**2)*k21 + (a**3)*k22 - (b**2)*k12 - (b**2)*k21 + (b**3)*k22 + &
                      R*a*k11 - R*b*k11 + 2*a*b*k21 - R*(a**2)*k12 + R*(a**2)*k21 - R*(a**3)*k22 + R*(b**2)*k12 + &
                      R*(b**2)*k21 - R*(b**3)*k22 - a*(b**2)*k22 - (a**2)*b*k22 + R*a*(b**2)*k22 + R*(a**2)*b*k22 - &
                      2*R*a*b*k21)/bottom 
                C4 = -C3-1
                
            ENDIF
            
        ELSEIF (basis_index ==3) THEN
            
            IF (piece_flag ==1) THEN
                C1 = 0 
                C2 = (k11 - R*k11 - a*k11 + a*k12 + a*k21 + b*k12 - b*k21 - (a**2)*k12 - (a**2)*k21 + (a**2)*k22 - &
                     (a**3)*k22 - (b**2)*k22 + R*a*k11 - R*a*k12 - R*a*k21 - R*b*k12 + R*b*k21 - a*b*k12 + a*b*k21 + &
                      R*(a**2)*k12 + R*(a**2)*k21 - R*(a**2)*k22 + R*(a**3)*k22 + R*(b**2)*k22 + a*(b**2)*k22 - &
                      R*a*(b**2)*k22 + R*a*b*k12 - R*a*b*k21)/bottom 
                C3 = 0 
                C4 = (2*R*k11 + 2*a*k11 + 2*a*k12 - 2*b*k12 + 2*(a**2)*k21 + 2*(a**2)*k22 + 2*(b**2)*k22 - 2*R*a*k11 + &
                      2*R*a*k21 - 2*R*b*k21 - 2*a*b*k21 - 4*a*b*k22 - 2*R*(a**2)*k21 + 2*R*a*b*k21)/bottom
               
            ELSEIF (piece_flag ==2) THEN
                C1 = a*(k11 - R*k11 + a*k12 + a*k21 + b*k12 - b*k21 + (a**2)*k22 - (b**2)*k22 - R*a*k12 - R*a*k21 - &
                     R*b*k12 + R*b*k21 - R*(a**2)*k22 + R*(b**2)*k22)/bottom 
                C2 = -C1 
                C3 = (b*k11 - a*k11 - (a**2)*k12 - (a**2)*k21 - (a**3)*k22 + (b**2)*k12 - (b**2)*k21 - (b**3)*k22 + &
                      R*a*k11 - R*b*k11 + 2*a*b*k21 + R*(a**2)*k12 + R*(a**2)*k21 + R*(a**3)*k22 - R*(b**2)*k12 + &
                      R*(b**2)*k21 + R*(b**3)*k22 + a*(b**2)*k22 + (a**2)*b*k22 - R*a*(b**2)*k22 - R*(a**2)*b*k22 - &
                      2*R*a*b*k21)/bottom 
                C4 = 1-C3
                
            ENDIF
            
        ELSEIF (basis_index ==4) THEN
            
            IF (piece_flag ==1) THEN
                C1 = 0 
                C2 = (2*k12 - k11 + R*k11 - 2*R*k12 + a*k11 - 3*a*k12 - a*k21 + 2*a*k22 - b*k12 + b*k21 - 2*b*k22 + &
                     (a**2)*k12 + (a**2)*k21 - 3*(a**2)*k22 + (a**3)*k22 + (b**2)*k22 - R*a*k11 + 3*R*a*k12 + R*a*k21 - &
                      2*R*a*k22 + R*b*k12 - R*b*k21 + 2*R*b*k22 + a*b*k12 - a*b*k21 + 2*a*b*k22 - R*(a**2)*k12 - &
                      R*(a**2)*k21 + 3*R*(a**2)*k22 - R*(a**3)*k22 - R*(b**2)*k22 - a*(b**2)*k22 + R*a*(b**2)*k22 - &
                      R*a*b*k12 + R*a*b*k21 - 2*R*a*b*k22)/bottom 
                C3 = 1 
                C4 = (2*R*a*k11 - 2*a*k11 - 2*(a**2)*k21 - 2*R*k11 - 2*R*a*k12 - 2*R*a*k21 + 2*R*b*k12 + 2*R*b*k21 + &
                      2*a*b*k21 + 2*R*(a**2)*k21 - 2*R*(a**2)*k22 - 2*R*(b**2)*k22 - 2*R*a*b*k21 + 4*R*a*b*k22)/bottom

            ELSEIF (piece_flag ==2) THEN
                C1 = -a*(k11 - 2*k12 - R*k11 + 2*R*k12 + a*k12 + a*k21 - 2*a*k22 + b*k12 - b*k21 + 2*b*k22 + (a**2)*k22 - &
                     (b**2)*k22 - R*a*k12 - R*a*k21 + 2*R*a*k22 - R*b*k12 + R*b*k21 - 2*R*b*k22 - R*(a**2)*k22 + &
                      R*(b**2)*k22)/bottom 
                C2 = -C1 
                C3 = (2*R*k11 + 2*a*k11 + 2*(a**2)*k21 - 2*R*a*k11 + 2*R*a*k12 + 2*R*a*k21 - 2*R*b*k12 - 2*R*b*k21 - &
                      2*a*b*k21 - 2*R*(a**2)*k21 + 2*R*(a**2)*k22 + 2*R*(b**2)*k22 + 2*R*a*b*k21 - 4*R*a*b*k22)/bottom 
                C4 = -C3
               
            ENDIF
            
        ENDIF 
        
    ENDIF
   
    IF (derivative_degree_x ==0 .AND. derivative_degree_y ==0) THEN       
        r1 = C1+C2*x+C3*y+C4*x*y                  
    ELSEIF (derivative_degree_x ==1 .AND. derivative_degree_y ==0) THEN
        r1 = C2+C4*y 
    ELSEIF (derivative_degree_x ==0 .AND. derivative_degree_y ==1) THEN
        r1 = C3+C4*x 
    ELSEIF (derivative_degree_x ==1 .AND. derivative_degree_y ==1) THEN
        !r1 = C4+x-x 
        r1 = C4 
    ENDIF
                
ENDIF

!print*,'r1=',r1
END
!IF (basis_type == 1) THEN
!
!    R = beta1/beta2 
!
!    IF (interface_element_type ==1) THEN
!        
!        bottom = 2*R*(b**2)-2*R*(b**2)*a+2*R*(a**2)-2*R*(a**2)*b+R*(a**3)*b+R*(b**3)*a+2*(b**2)*a+2*(a**2)*b-(a**3)*b-(b**3)*a 
!    
!        IF (basis_index ==1) THEN
!            
!            IF (piece_flag ==1) THEN
!                C1 = 1 
!                C2 = -(-2*R*b*a+R*(a**2)*b+R*(b**3)+2*(b**2)+2*a*b-(a**2)*b-(b**3)+2*R*(a**2))/bottom 
!                C3 = -(-2*R*b*a+R*(a**3)+R*(b**2)*a+2*a*b+2*(a**2)-(a**3)-(b**2)*a+2*R*(b**2))/bottom 
!                C4 = 2*R*((a**2)+(b**2))/bottom
!             
!            ELSEIF (piece_flag ==2) THEN
!                C1 = 2*R*((a**2)+(b**2))/bottom 
!                C2 = -C1 
!                C3 = -C1 
!                C4 = C1
!              
!            ENDIF
!            
!        ELSEIF (basis_index ==2) THEN
!                             
!            IF (piece_flag ==1) THEN
!                C1 = 0 
!                C2 = -(R*(a**2)*b-R*(b**3)-2*(b**2)-(a**2)*b+(b**3)-2*R*(a**2))/bottom 
!                C3 = -a*(2*R*b-R*(a**2)-3*R*(b**2)-2*b+(a**2)+3*(b**2)+R*(a**2)*b+R*(b**3)-(a**2)*b-(b**3))/bottom 
!                C4 = 2*(-R*(a**2)-R*(b**2)+R*(a**2)*b-(a**2)*b)/bottom
!        
!            ELSEIF (piece_flag ==2) THEN
!                C1 = a*b*(-2*R*b+R*(a**2)+R*(b**2)+2*b-(a**2)-(b**2))/bottom 
!                C2 = 1-C1 
!                C3 = -C1 
!                C4 = C1-1
!               
!            ENDIF
!            
!        ELSEIF (basis_index ==3) THEN
! !           piece_flag=2
!            IF (piece_flag ==1) THEN
!                C1 = 0 
!                C2 = b*(-(b**2)+R*(a**2)+R*(b**2)-(a**2))*(-1+a)/bottom 
!                C3 = a*(-(b**2)+R*(a**2)+R*(b**2)-(a**2))*(b-1)/bottom 
!                C4 = -2*(-R*(b**2)+R*(b**2)*a-R*(a**2)+R*(a**2)*b-(b**2)*a-(a**2)*b)/bottom
!   
!            ELSEIF (piece_flag ==2) THEN
!                C1 = -a*b*(-(b**2)+R*(a**2)+R*(b**2)-(a**2))/bottom 
!                C2 = -C1 
!                C3 = -C1 
!                C4 = 1+C1
!               
!            ENDIF
!            
!        ELSEIF (basis_index ==4) THEN
!            
!            IF (piece_flag ==1) THEN
!                C1 = 0 
!                C2 = -b*(2*R*a-3*R*(a**2)-R*(b**2)-2*a+3*(a**2)+(b**2)+R*(a**3)+R*(b**2)*a-(a**3)-(b**2)*a)/bottom 
!                C3 = (R*(a**3)-R*(b**2)*a+2*(a**2)-(a**3)+(b**2)*a+2*R*(b**2))/bottom 
!                C4 = -2*(R*(b**2)+(b**2)*a-R*(b**2)*a+R*(a**2))/bottom
!
!            ELSEIF (piece_flag ==2) THEN
!                C1 = a*b*(-2*R*a+R*(a**2)+R*(b**2)+2*a-(a**2)-(b**2))/bottom 
!                C2 = -C1 
!                C3 = 1-C1 
!                C4 = C1-1
!              
!            ENDIF
!            
!        ENDIF     
!        
!    ELSEIF (interface_element_type ==2) THEN
!        
!        bottom = R*(a**3)-(a**3)-R*(a**2)*b+(a**2)*b+2*(a**2)-R*(b**2)*a-R*a-4*a*b+(b**2)*a+a+2*R+R*(b**3)-R*b+b-(b**3)+2*(b**2) 
!!         bottom = -a-2*(a**2)+a**3-b+4*a*b-(a**2)*b-2*(b**2)-a*(b**2)+b**3 
!        IF (basis_index ==1) THEN
!            
!            IF (piece_flag ==1) THEN
!                C1 = 1 
!                C2 = (R*(a**2)*b-3*R*(a**2)-(a**2)*b+(a**2)+2*b*R*a+2*R*a-2*a+ &
!                        2*a*b+R*(b**2)-R*(b**3)-R-R*b+(b**3)+b-3*(b**2)-1)/bottom 
!                C3 = -1 
!                C4 = 2*(R*(a**2)-2*b*R*a+R*(b**2)+R-R*b+b)/bottom
!
!            ELSEIF (piece_flag ==2) THEN
!                C1 = -(-(a**2)*b+R*b-b+2*b*R*a-2*R*(a**2)-2*R+2*a*b-2*(b**2)+R*(a**2)*b-R*(b**3)+(b**3))/bottom 
!                C2 = -C1 
!                C3 = -2*(R*(a**2)-2*b*R*a+R*(b**2)+R-R*b+b)/bottom 
!                C4 = -C3
!                
!            ENDIF
!            
!        ELSEIF (basis_index ==2) THEN
!            
!            IF (piece_flag ==1) THEN
!                C1 = 0 
!                C2 = -(-R*(a**2)+R*(a**2)*b-(a**2)-(a**2)*b+4*a*b-R*(b**3)+R*(b**2)+R*b-R-b-3*(b**2)+(b**3)-1)/bottom 
!                C3 = 0 
!                C4 = -2*((a**2)-2*a*b+R-R*b+(b**2)+b)/bottom
!
!            ELSEIF (piece_flag ==2) THEN
!                C1 = a*(1-R+(b**2)+R*(a**2)-(a**2)-R*(b**2))/bottom 
!                C2 = 1-C1 
!                C3 = -(R*(a**3)-(a**3)-R*(a**2)*b+(a**2)*b-R*(b**2)*a-R*a+a+(b**2)*a+R*(b**3)+R*b-(b**3)-b)/bottom 
!                C4 = -C3-1
!                
!            ENDIF
!            
!        ELSEIF (basis_index ==3) THEN
!            
!            IF (piece_flag ==1) THEN
!                C1 = 0 
!                C2 = (-(a**3)+R*(a**3)-R*(a**2)+(a**2)+R*a-R*(b**2)*a+(b**2)*a-a-R+R*(b**2)-(b**2)+1)/bottom 
!                C3 = 0 
!                C4 = 2*((a**2)-R*a-2*a*b+a+R+(b**2))/bottom
!               
!            ELSEIF (piece_flag ==2) THEN
!                C1 = -a*(-1+R+(b**2)+R*(a**2)-(a**2)-R*(b**2))/bottom 
!                C2 = -C1 
!                C3 = (R*a-a+R*(a**3)-(a**3)-R*(a**2)*b+(a**2)*b-R*(b**2)*a+(b**2)*a+R*(b**3)-R*b+b-(b**3))/bottom 
!                C4 = 1-C3
!                
!            ENDIF
!            
!        ELSEIF (basis_index ==4) THEN
!            
!            IF (piece_flag ==1) THEN
!                C1 = 0 
!                C2 = -(-(a**3)+R*(a**3)-3*R*(a**2)+3*(a**2)+2*b*R*a+3*R*a- &
!                    R*(b**2)*a-2*a*b-3*a+(b**2)*a+R*(b**2)-R-2*R*b-(b**2)+1+2*b)/bottom 
!                C3 = 1 
!                C4 = -2*(R*(a**2)-R*a-2*b*R*a+a+R*(b**2)+R)/bottom
!
!            ELSEIF (piece_flag ==2) THEN
!                C1 = a*(-1+R+2*R*b+(b**2)-2*R*a-2*b+2*a+R*(a**2)-(a**2)-R*(b**2))/bottom 
!                C2 = -C1 
!                C3 = 2*(R*(a**2)-R*a-2*b*R*a+a+R*(b**2)+R)/bottom 
!                C4 = -C3
!               
!            ENDIF
!            
!        ENDIF 
!        
!    ENDIF
!   
!    IF (derivative_degree_x ==0 .AND. derivative_degree_y ==0) THEN       
!        r1 = C1+C2*x+C3*y+C4*x*y                  
!    ELSEIF (derivative_degree_x ==1 .AND. derivative_degree_y ==0) THEN
!        r1 = C2+C4*y 
!    ELSEIF (derivative_degree_x ==0 .AND. derivative_degree_y ==1) THEN
!        r1 = C3+C4*x 
!    ELSEIF (derivative_degree_x ==1 .AND. derivative_degree_y ==1) THEN
!        r1 = C4+x-x 
!    ENDIF
!                
!ENDIF
!
!
!END