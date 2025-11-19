SUBROUTINE Function_coefficient_2D(coefficient_function_name_impic,coefficient_function_name,x,y, &
                                 trial_derivative_degree_x,trial_derivative_degree_y,&
                                 test_derivative_degree_x,test_derivative_degree_y, element_Gauss)
    
!!! ************************ bjw add for impic 2019-6-3 **********************************************
!!! ***** 2019-6-21 bjw 修改，认为在 Function_coefficient_2D.f90 里加 Epsilon更合理 ******
!!! ***** 有磁场时只在对角线上加 Epsilon ，此处要注意 ******
!=========LY modification, 2022-7-25=========
Use IFE_Data
!=========LY modification, 2022-7-25=========
USE Domain_2D  
USE IMPIC_Data_2D
IMPLICIT NONE
REAL(8)   :: coefficient_function_name
REAL(8)   :: coefficient_function_name_impic
REAL(8)   :: x,y
INTEGER   :: trial_derivative_degree_x,trial_derivative_degree_y, &
		     test_derivative_degree_x,test_derivative_degree_y
REAL(8)   :: rxp, ryp, dx, dy, f, g
INTEGER   :: i, j
!$ ==================== ab.ZWZ 2021/4/12=======================\\
REAL(8) :: xcellmdx, ycellmdy
REAL(8) :: R1, R2, den
REAL(8) :: P1, P2, P3, P4
REAL(8) :: BETA = 1.
!INTEGER :: delta 
REAL(8) :: Chi_p
!$ ==================== ab.ZWZ 2021/4/12=======================//
!=========LY modification, 2022-7-25=========
Integer :: element_Gauss
Real(8) :: hx_partition, hy_partition
Real(8), Parameter :: SmallValue_LY = 1.0D-12
!=========LY modification, 2022-7-25=========

!=========LY modification, 2022-7-25=========
coefficient_function_name_impic = 0.0

hx_partition = HP(1, HT(2, element_Gauss)) - HP(1, HT(1, element_Gauss))
hy_partition = HP(2, HT(4, element_Gauss)) - HP(2, HT(1, element_Gauss))

dx = (x - HP(1, HT(1, element_Gauss))) / hx_partition
dy = (y - HP(2, HT(1, element_Gauss))) / hy_partition
xcellmdx = 1.0 - dx
ycellmdy = 1.0 - dy

!$ ==================== mb.ZWZ 2021/4/12=======================\\
IF( delta_global == 0 ) THEN
  P1 = xcellmdx * ycellmdy    !local 1
  P2 = dx       * ycellmdy    !local 2
  P3 = dx       * dy          !local 3
  P4 = xcellmdx * dy          !local 4
ELSEIF( delta_global == 1 ) THEN
  If (ABS(HP(2, HT(1, element_Gauss))-dymin) < SmallValue_LY) Then
    BETA = 0.75
  Else
    BETA =1.0
  End If
  R1 = f_left_wall(2) + HP(2, HT(1, element_Gauss))
  R2 = f_left_wall(2) + HP(2, HT(4, element_Gauss))
  den=R2*R2-R1*R1
    
  P1 = xcellmdx * (R2*R2-y*y) / den * BETA    !local 1
  P2 = dx       * (R2*R2-y*y) / den * BETA    !local 2
  P3 = dx       * (y*y-R1*R1) / den           !local 3
  P4 = xcellmdx * (y*y-R1*R1) / den           !local 4
ENDIF

IF (Bfiled_index) THEN   !!!! 有磁场
    
    IF (trial_derivative_degree_x==1 .AND. trial_derivative_degree_y==0 .AND. & 
        test_derivative_degree_x==1 .AND. test_derivative_degree_y==0) THEN !!! xx, C11

        Chi_p = TransChi(1,HT(1,element_Gauss),1)*P1 + TransChi(1,HT(2,element_Gauss),1)*P2 + &
                TransChi(1,HT(3,element_Gauss),1)*P3 + TransChi(1,HT(4,element_Gauss),1)*P4
        
        !coefficient_function_name_impic = 2.0
        !=========IFE=========
        !coefficient_function_name_impic = Chi_p*coefficient_function_name + coefficient_function_name   !$ when Chi is given specific value in checking convergence
        !=========IFE=========
        !=========PIC=========
        coefficient_function_name_impic = Chi_p + coefficient_function_name
        !=========PIC=========
    
    ELSEIF (trial_derivative_degree_x==1 .AND. trial_derivative_degree_y==0 .AND. & 
            test_derivative_degree_x==0 .AND. test_derivative_degree_y==1) THEN !!! xy, C21
    
        Chi_p = TransChi(4,HT(1,element_Gauss),1)*P1 + TransChi(4,HT(2,element_Gauss),1)*P2 + &
                TransChi(4,HT(3,element_Gauss),1)*P3 + TransChi(4,HT(4,element_Gauss),1)*P4
        
        !coefficient_function_name_impic = 2.0
        !=========IFE=========
        !coefficient_function_name_impic = Chi_p*coefficient_function_name + 0.0
        !=========IFE=========
        !=========PIC=========
        coefficient_function_name_impic = Chi_p + 0.0
        !=========PIC=========
    
    ELSEIF (trial_derivative_degree_x==0 .AND. trial_derivative_degree_y==1 .AND. & 
            test_derivative_degree_x==1 .AND. test_derivative_degree_y==0) THEN !!! yx, C12

        Chi_p = TransChi(2,HT(1,element_Gauss),1)*P1 + TransChi(2,HT(2,element_Gauss),1)*P2 + &
                TransChi(2,HT(3,element_Gauss),1)*P3 + TransChi(2,HT(4,element_Gauss),1)*P4
        
        !coefficient_function_name_impic = 2.0
        !=========IFE=========
        !coefficient_function_name_impic = Chi_p*coefficient_function_name + 0.0
        !=========IFE=========
        !=========PIC=========
        coefficient_function_name_impic = Chi_p + 0.0
        !=========PIC=========
    
    ELSEIF (trial_derivative_degree_x==0 .AND. trial_derivative_degree_y==1 .AND. & 
            test_derivative_degree_x==0 .AND. test_derivative_degree_y==1) THEN !!! yy, C22
    
        Chi_p = TransChi(5,HT(1,element_Gauss),1)*P1 + TransChi(5,HT(2,element_Gauss),1)*P2 + &
                TransChi(5,HT(3,element_Gauss),1)*P3 + TransChi(5,HT(4,element_Gauss),1)*P4
    
        !coefficient_function_name_impic = 5.0
        !=========IFE=========
        !coefficient_function_name_impic = Chi_p*coefficient_function_name + coefficient_function_name
        !=========IFE=========
        !=========PIC=========
        coefficient_function_name_impic = Chi_p + coefficient_function_name
        !=========PIC=========
        
    ENDIF
        
ELSE !!!! 无磁场

    
    !$ =============== mb.ZWZ 2021/7/8 ================ \\
    Chi_p = Chi(HT(1,element_Gauss),1)*P1 + Chi(HT(2,element_Gauss),1)*P2 + &
            Chi(HT(3,element_Gauss),1)*P3 + Chi(HT(4,element_Gauss),1)*P4
    
    !coefficient_function_name_impic = (f + dy * (g-f)) + coefficient_function_name
    !Chi_p = OneAndChi(i,j)*P1 + OneAndChi(i,j)*P2 + OneAndChi(i,j)*P3 + OneAndChi(i,j)*P4
    
    !=========IFE=========
    !coefficient_function_name_impic = Chi_p*coefficient_function_name + coefficient_function_name
    !=========IFE=========
        
    !=========PIC=========
    coefficient_function_name_impic = Chi_p + coefficient_function_name
    !=========PIC=========
    !$ =============== mb.ZWZ 2021/7/8 ================ //
    
END IF
!$ ==================== mb.ZWZ 2021/4/12=======================//
!=========LY modification, 2022-7-25=========

!IF (Bfiled_index) THEN   !!!! 有磁场
!    
!    IF (trial_derivative_degree_x==1 .AND. trial_derivative_degree_y == 0 .AND. & 
!        test_derivative_degree_x==1 .AND. test_derivative_degree_y==0) THEN !!! xx, C11
!
!        f = TransChi(1,i,j)+dx*(TransChi(1,i+1,j)-TransChi(1,i,j))
!        g = TransChi(1,i,j+1)+dx*(TransChi(1,i+1,j+1)-TransChi(1,i,j+1))
!        !coefficient_function_name_impic = (f + dy * (g-f)) * coefficient_function_name
!        coefficient_function_name_impic = (f + dy * (g-f)) + coefficient_function_name
!    
!    ELSEIF (trial_derivative_degree_x==1 .AND. trial_derivative_degree_y == 0 .AND. & 
!        test_derivative_degree_x==0 .AND. test_derivative_degree_y==1) THEN !!! xy, C12
!    
!        f = TransChi(2,i,j)+dx*(TransChi(2,i+1,j)-TransChi(2,i,j))
!        g = TransChi(2,i,j+1)+dx*(TransChi(2,i+1,j+1)-TransChi(2,i,j+1))
!        !coefficient_function_name_impic = (f + dy * (g-f)) * coefficient_function_name
!        coefficient_function_name_impic = (f + dy * (g-f)) + 0.0
!    
!    ELSEIF (trial_derivative_degree_x==0 .AND. trial_derivative_degree_y == 1 .AND. & 
!        test_derivative_degree_x==1 .AND. test_derivative_degree_y==0) THEN !!! yx, C21
!    
!        f = TransChi(4,i,j)+dx*(TransChi(4,i+1,j)-TransChi(4,i,j))
!        g = TransChi(4,i,j+1)+dx*(TransChi(4,i+1,j+1)-TransChi(4,i,j+1))
!        !coefficient_function_name_impic = (f + dy * (g-f)) * coefficient_function_name
!        coefficient_function_name_impic = (f + dy * (g-f)) + 0.0
!    
!    ELSEIF (trial_derivative_degree_x==0 .AND. trial_derivative_degree_y == 1 .AND. & 
!        test_derivative_degree_x==0 .AND. test_derivative_degree_y==1) THEN !!! yy, C22
!    
!        f = TransChi(5,i,j)+dx*(TransChi(5,i+1,j)-TransChi(5,i,j))
!        g = TransChi(5,i,j+1)+dx*(TransChi(5,i+1,j+1)-TransChi(5,i,j+1))
!        !coefficient_function_name_impic = (f + dy * (g-f)) * coefficient_function_name
!        coefficient_function_name_impic = (f + dy * (g-f)) + coefficient_function_name
!    
!    ENDIF
!        
!ELSE !!!! 无磁场
!    
!    f = OneAndChi(i,j)+dx*(OneAndChi(i+1,j)-OneAndChi(i,j))
!    g = OneAndChi(i,j+1)+dx*(OneAndChi(i+1,j+1)-OneAndChi(i,j+1))
!    !coefficient_function_name_impic = (f + dy * (g-f)) * coefficient_function_name
!    coefficient_function_name_impic = (f + dy * (g-f)) + coefficient_function_name
!    
!END IF

END
    