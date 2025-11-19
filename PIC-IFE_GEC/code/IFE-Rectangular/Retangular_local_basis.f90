SUBROUTINE  Retangular_local_basis( x, y, left_lower_point, h_partition, basis_type, &
                                basis_index, derivative_degree_x, derivative_degree_y, r)


IMPLICIT NONE

REAL(8)                                               x, y, kesai, eita
REAL(8)	                                              h_partition(2) 
INTEGER                                               basis_type, basis_index, derivative_degree_x, derivative_degree_y
REAL(8)                           	        	      left_lower_point(2)  
REAL(8)                                               r

kesai=(2*x-2*left_lower_point(1)-h_partition(1))/h_partition(1)
eita=(2*y-2*left_lower_point(2)-h_partition(2))/h_partition(2)

IF (derivative_degree_x==0 .AND. derivative_degree_y==0) THEN

 CALL Retangular_reference_basis(kesai,eita,basis_type,basis_index,0,0,r)

ELSEIF (derivative_degree_x==1 .AND. derivative_degree_y==0) THEN

 CALL Retangular_reference_basis(kesai,eita,basis_type,basis_index,1,0,r)

   r = r*2/h_partition(1)

ELSEIF (derivative_degree_x==0 .AND. derivative_degree_y==1) THEN

 CALL Retangular_reference_basis(kesai,eita,basis_type,basis_index,0,1,r)
   r = r*2/h_partition(2)

ELSEIF (derivative_degree_x==1 .AND. derivative_degree_y==1) THEN

 CALL Retangular_reference_basis(kesai,eita,basis_type,basis_index,1,1,r)
   r = r*4/(h_partition(1)*h_partition(2))


ENDIF

END


