SUBROUTINE Compute_IFE_Norm_error(delta, Gauss_coefficient_reference, Gauss_point_reference, Gauss_point_reference_triangle, &
                                 Gauss_coefficient_reference_triangle, derivative_degree_x, derivative_degree_y, error_norm)
    
  !!! bjw add 
USE Domain_2D
USE Field_2D
USE IFE_Data    
USE IFE_MAIN_PARAM
USE Object_Data_2D
USE Gauss_Data
USE IFE_INTERFACE, ONLY: Retangular_local_basis, Generate_Gauss_local, Gauss_integration_local_error_IFE
IMPLICIT NONE    
INTEGER, INTENT(IN) :: delta !$ ab.ZWZ 2021/7/7
INTEGER     :: i,j,k, derivative_degree_x, derivative_degree_y, nindex, i_Object, N_Objects
REAL(8)     :: A,temp1,temp2,error,error_norm,error1,error2

REAL(8),DIMENSION(:),POINTER                ::  Gauss_coefficient_reference
REAL(8),DIMENSION(:,:),POINTER              ::  Gauss_point_reference
REAL(8),DIMENSION(:),POINTER                ::  Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:),POINTER              ::  Gauss_point_reference_triangle
REAL(8),DIMENSION(:),POINTER                ::  Gauss_coefficient_local
REAL(8),DIMENSION(:,:),POINTER              ::  Gauss_point_local
INTEGER                                           information_vector_1(18)
!REAL(8)                                           information_vector_2(12)
REAL(8)                                           information_vector_2(8)       !$ mb.ZWZ 2021/7/7
REAL(8)                                           vertices(2,4)

REAL(8),DIMENSION(:),ALLOCATABLE    :: uh, u_true, error_ele, U_fe
TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
REAL(8)                   plus_coefficient_function_name, temp



ALLOCATE(Gauss_coefficient_local(9), Gauss_point_local(9,2))


ALLOCATE(uh(9), u_true(9))
ALLOCATE(error_ele(SIZE(HT,2)))
ALLOCATE(U_fe(SIZE(HP,2)))

! ---- wsy revise for SIDG----
DO nindex = 1, SIZE(HP,2)
    U_fe(nindex) = Phi(nindex,1)
END DO
!------------------------

error = 0.
DO i=1,SIZE(HT,2)
    
    h_partition(1) = HP(1,HT(2,i)) - HP(1,HT(1,i))
    h_partition(2) = HP(2,HT(4,i)) - HP(2,HT(1,i))
    
   CALL Generate_Gauss_local(Gauss_coefficient_reference,Gauss_point_reference, &
               HP(1:2,HT(1,i)), h_partition, Gauss_coefficient_local, Gauss_point_local)
   
   uh = 0.0
   error = 0.0
   plus_coefficient_function_name = Global_Beta(1)
   IF (element_index(i) == -1) THEN  !OUTSIDE THE OBJECTS
       DO k = 1, SIZE(Gauss_coefficient_local)    
           DO j = 1, 4
       
                CALL Retangular_local_basis( Gauss_point_local(k,1),Gauss_point_local(k,2), HP(1:2,HT(1,i)), &
                                                h_partition, 1, j,&
                                                derivative_degree_x, derivative_degree_y, temp1)

               uh(k) = uh(k) + U_fe(HT(j,i)) * temp1
    
           END DO
           
           CALL Function_True(delta,plus_coefficient_function_name, Gauss_point_local(k,1),Gauss_point_local(k,2),&
                             derivative_degree_x, derivative_degree_y, temp2)
           u_true(k) = temp2
           
           !error = error + Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k))
           !$ ============== ab.ZWZ for adding cylindrical coordinate ===== \\
           IF (delta == 0) THEN
                error = error + Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k))
           ELSEIF (delta == 1) THEN
                error = error + Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k)) *Gauss_point_local(k,2)
           ENDIF
           !$ ============== ab.ZWZ for adding cylindrical coordinate ===== //
           
           error2 =  error2 + temp
           
       END DO    
   ELSEIF (element_index(i) < -1) THEN  !INSIDE THE OBJECTS

       plus_coefficient_function_name = Global_Beta(2)
 
       DO k = 1, SIZE(Gauss_coefficient_local)    
           DO j = 1, 4
       
               CALL Retangular_local_basis( Gauss_point_local(k,1),Gauss_point_local(k,2), HP(1:2,HT(1,i)), &
                                                h_partition, 1, j, &
                                                derivative_degree_x, derivative_degree_y, temp1)

               uh(k) = uh(k) + U_fe(HT(j,i)) * temp1
       
           END DO
       
           CALL Function_True(delta,plus_coefficient_function_name, Gauss_point_local(k,1),Gauss_point_local(k,2), &
                                derivative_degree_x, derivative_degree_y, temp2)
           u_true(k) = temp2
           
           !error = error + Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k))
           !$ ============== ab.ZWZ for adding cylindrical coordinate ===== \\
           IF (delta == 0) THEN
                error = error + Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k))
           ELSEIF (delta == 1) THEN
                error = error + Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k)) *Gauss_point_local(k,2)
           ENDIF
           !$ ============== ab.ZWZ for adding cylindrical coordinate ===== //
           
           error2 =  error2 + temp
       END DO   
   ELSEIF (element_index(i) > 0) THEN  !THE INTERSECT ELEMENTS
       vertices = HP(1:2,HT(1:4,i))
       information_vector_1 = information_1(:,element_index(i))
       information_vector_2 = information_2(:,element_index(i))    
       
       !CALL Gauss_integration_local_error_IFE(0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
       CALL Gauss_integration_local_error_IFE(delta,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
			                                    vertices,information_vector_1,information_vector_2,i, &
												1, derivative_degree_x, derivative_degree_y, U_fe, temp, i)
      
       error = temp
       
       error1 =  error1 + temp

   END IF
   
 
   
   error_ele(i) = error
   
END DO   

error_norm = DSQRT(SUM(error_ele(:)))
WRITE(*,*) '++++',MAXVAL(error_ele(:)),  error_norm

   

END
   
    
    
    
! ALLOCATE(Gauss_coefficient_local(9), Gauss_point_local(9,2))
!
!
!ALLOCATE(uh(9), u_true(9))
!ALLOCATE(error_ele(SIZE(t_c,2)))
!ALLOCATE(U_fe(nx*ny))
!
!DO i=1,nx
!	DO j=1,ny
!			nindex = j+(i-1)*ny
!			U_fe(nindex) = Phi(i,j)
!	END DO
!END DO
!WRITE(*,*) '**************'
!WRITE(*,*) U_fe
!WRITE(*,*) '**************'
!A = hx(1) * hx(2) !!! element area
!error = 0.
!DO i=1,SIZE(t_c,2)
!
!   CALL Generate_Gauss_local(Gauss_coefficient_reference,Gauss_point_reference, &
!               p_basic(1:2,t_c(1,i)), hx, Gauss_coefficient_local, Gauss_point_local)
!   
!   uh = 0.0
!   error = 0.0
!   plus_coefficient_function_name = Global_Beta(1)
!   IF (element_index(i) == -1) THEN  !OUTSIDE THE OBJECTS
!       DO k = 1, SIZE(Gauss_coefficient_local)    
!           DO j = 1, 4
!       
!               CALL Retangular_local_basis( Gauss_point_local(k,1),Gauss_point_local(k,2), p_basic(1:2,t_c(1,i)), &
!                                                h_partition, 1, j,&
!                                                derivative_degree_x, derivative_degree_y, temp1)
!       
!               !r=r+Gauss_coefficient_local(i) * temp
!               uh(k) = uh(k) + U_fe(t_c(j,i)) * temp1
!       
!               CALL Function_True(plus_coefficient_function_name, Gauss_point_local(k,1),Gauss_point_local(k,2),&
!                                  derivative_degree_x, derivative_degree_y, temp2)
!               u_true(k) = temp2
!       
!           END DO
!       
!           error = error + Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k))
!           !error = error + A * Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k))
!           error2 =  error2 + temp
!       END DO    
!   ELSEIF (element_index(i) < -1) THEN  !INSIDE THE OBJECTS
!
!       plus_coefficient_function_name = Global_Beta(2)
! 
!       DO k = 1, SIZE(Gauss_coefficient_local)    
!           DO j = 1, 4
!       
!               CALL Retangular_local_basis( Gauss_point_local(k,1),Gauss_point_local(k,2), p_basic(1:2,t_c(1,i)), &
!                                                h_partition, 1, j, &
!                                                derivative_degree_x, derivative_degree_y, temp1)
!       
!               !r=r+Gauss_coefficient_local(i) * temp
!               uh(k) = uh(k) + U_fe(t_c(j,i)) * temp1
!       
!               CALL Function_True(plus_coefficient_function_name, Gauss_point_local(k,1),Gauss_point_local(k,2), &
!                                  derivative_degree_x, derivative_degree_y, temp2)
!               u_true(k) = temp2
!       
!           END DO
!       
!           error = error + Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k))
!           !error = error + A * Gauss_coefficient_local(k) * (uh(k)-u_true(k)) * (uh(k)-u_true(k))
!           error2 =  error2 + temp
!       END DO   
!   ELSEIF (element_index(i) > 0) THEN  !THE INTERSECT ELEMENTS
!       vertices = p_basic(1:2,t_c(1:4,i))
!       information_vector_1 = information_1(:,element_index(i))
!       information_vector_2 = information_2(:,element_index(i))    
!       
!       CALL Gauss_integration_local_error_IFE(0,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle, &
!			                                    vertices,information_vector_1,information_vector_2,i, &
!												1, derivative_degree_x, derivative_degree_y, U_fe, temp)
!      
!       error = temp
!       
!       error1 =  error1 + temp
!       !WRITE(*,*) 'iii',i,error
!       !STOP
!   END IF
!   
!   !IF (i==2197) THEN
!   !     WRITE(*,*) '++++iiiii',i,error
!   !     STOP
!   !ENDIF
!   
!   
!   error_ele(i) = error
!
!   !IF (error>1.0D-3)THEN
!   !    WRITE(*,*) i, element_index(i), error
!   !ENDIF
!   
!END DO   
!
!WRITE(*,*) '++++', error1,error2
!PAUSE
!error_norm = DSQRT(SUM(error_ele(:)))
!WRITE(*,*) '++++',MAXVAL(error_ele(:)),  error_norm
!
!   
!
!END
      
    
    
    
    