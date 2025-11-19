SUBROUTINE	retangular_local_basis_IFE_nonhomogeneous(x,y,vertices,information_vector_1,	&
														information_vector_2,piece_flag,basis_type,derivative_degree_x,derivative_degree_y,r)

USE IFE_INTERFACE, ONLY: Left_Divide
IMPLICIT NONE

!INCLUDE	'Left_Divide.inc'			

REAL(8)										::	x, y
REAL(8)										::	vertices(2,4)
INTEGER										::	information_vector_1(18)
REAL(8)										::	information_vector_2(8)
INTEGER										::	piece_flag
INTEGER                                     ::  basis_type, derivative_degree_x, derivative_degree_y
REAL(8)										::	r

INTEGER										::	interface_element_type
INTEGER										::	pointer_reference_to_local(4)
INTEGER										::	plus_piece_flag
REAL(8)										::	beta1, beta2
REAL(8)										::	Dx, Dy, Ex, Ey
REAL(8)										::	nx, ny, magnitude
REAL(8)										::	x1, x2, x3, x4
REAL(8)										::	y1, y2, y3, y4
REAL(8)										::	R_beta
INTEGER										::	S1
REAL(8)										::	temp, rhs(8,1)
REAL(8)										::	coe_matrix(8,8), coe(8,1)
INTEGER										::	i,j
REAL(8)										::	C1, C2, C3, C4	

interface_element_type = information_vector_1(6)
pointer_reference_to_local(:) = information_vector_1(11:14)
plus_piece_flag = information_vector_1(15)

beta1 = information_vector_2(1)
beta2 = information_vector_2(2)
Dx = information_vector_2(3)
Dy = information_vector_2(4)
Ex = information_vector_2(5)
Ey = information_vector_2(6)
nx = information_vector_2(7)
ny = information_vector_2(8)

x1 = vertices(1,pointer_reference_to_local(1))
y1 = vertices(2,pointer_reference_to_local(1))
x2 = vertices(1,pointer_reference_to_local(2))
y2 = vertices(2,pointer_reference_to_local(2)) 
x3 = vertices(1,pointer_reference_to_local(3))
y3 = vertices(2,pointer_reference_to_local(3))    
x4 = vertices(1,pointer_reference_to_local(4))
y4 = vertices(2,pointer_reference_to_local(4))

R_beta = beta1 / beta2

IF(plus_piece_flag == 2)THEN
	S1 = 1
ELSEIF(plus_piece_flag == 1)THEN
	s1 = -1
endIF
 
 
temp = beta2 * S1 * SQRT((Ex - Dx)**2 + (Ey - Dy)**2)   



rhs(1,1) = 0
rhs(2,1) = 0
rhs(3,1) = 0
rhs(4,1) = 0
rhs(5,1) = 0
rhs(6,1) = 0
rhs(7,1) = 0
rhs(8,1) = 1 / (beta2 * S1 * SQRT((Ex - Dx)**2 + (Ey - Dy)**2))


DO i=1,8
	DO j=1,8
		coe_matrix(i,j) = 0.
	ENDDO
ENDDO

coe_matrix(1,1) = 1
coe_matrix(1,2) = x1
coe_matrix(1,3) = y1
coe_matrix(1,4) = x1 * y1

coe_matrix(2,5) = 1
coe_matrix(2,6) = x2
coe_matrix(2,7) = y2
coe_matrix(2,8) = x2 * y2

coe_matrix(3,5) = 1
coe_matrix(3,6) = x3
coe_matrix(3,7) = y3
coe_matrix(3,8) = x3 * y3

IF(interface_element_type == 1)THEN
	coe_matrix(4,5) = 1
	coe_matrix(4,6) = x4
	coe_matrix(4,7) = y4
	coe_matrix(4,8) = x4 * y4
ELSEIF(interface_element_type == 2)THEN
	coe_matrix(4,1) = 1
	coe_matrix(4,2) = x4
	coe_matrix(4,3) = y4
	coe_matrix(4,4) = x4 * y4
ENDIF

coe_matrix(5,1) = 1
coe_matrix(5,2) = Dx
coe_matrix(5,3) = Dy
coe_matrix(5,4) = Dx * Dy
coe_matrix(5,5) = -1
coe_matrix(5,6) = -Dx
coe_matrix(5,7) = -Dy
coe_matrix(5,8) = -Dx * Dy

coe_matrix(6,1) = 1
coe_matrix(6,2) = Ex
coe_matrix(6,3) = Ey
coe_matrix(6,4) = Ex * Ey
coe_matrix(6,5) = -1
coe_matrix(6,6) = -Ex
coe_matrix(6,7) = -Ey
coe_matrix(6,8) = -Ex * Ey

IF( ABS( (Dx-Ex)/MIN(Dx,Ex) )<1.0E-10 .OR. ABS( (Dy-Ey)/MIN(Dy,Ey) )<1.0E-10 )THEN
	coe_matrix(7,1) = 0
	coe_matrix(7,2) = 0
	coe_matrix(7,3) = 0
	coe_matrix(7,4) = 1
	coe_matrix(7,5) = 0
	coe_matrix(7,6) = 0
	coe_matrix(7,7) = 0
	coe_matrix(7,8) = -1
ELSE
	coe_matrix(7,1) = 1
	coe_matrix(7,2) = (Dx + Ex) / 2
	coe_matrix(7,3) = (Dy + Ey) / 2
	coe_matrix(7,4) = (Dx + Ex) * (Dy + Ey) / 4
	coe_matrix(7,5) = -1
	coe_matrix(7,6) = -(Dx + Ex) / 2
	coe_matrix(7,7) = -(Dy + Ey) / 2
	coe_matrix(7,8) = -(Dx + Ex) * (Dy + Ey) / 4
ENDIF

coe_matrix(8,1) = 0
coe_matrix(8,2) = -R_beta * nx
coe_matrix(8,3) = -R_beta * ny
coe_matrix(8,4) = -0.5 * R_beta * ((Dx + Ex) * ny + (Dy + Ey) * nx)
coe_matrix(8,5) = 0
coe_matrix(8,6) = nx
coe_matrix(8,7) = ny
coe_matrix(8,8) = 0.5 * ((Dx + Ex) * ny + (Dy + Ey) * nx)

CALL Left_Divide(coe, coe_matrix, rhs)

IF(piece_flag == 1)THEN
	C1 = coe(1,1)
	C2 = coe(2,1)
	C3 = coe(3,1)
	C4 = coe(4,1)
ELSEIF(piece_flag == 2)THEN
	C1 = coe(5,1)
	C2 = coe(6,1)
	C3 = coe(7,1)
	C4 = coe(8,1)
ENDIF


IF(derivative_degree_x == 0 .AND. derivative_degree_y == 0)THEN
	r = C1 + C2 * x + C3 * y + C4 * x * y
ELSEIF(derivative_degree_x == 1 .AND. derivative_degree_y == 0)THEN
	r = C2 + C4 * y
ELSEIF(derivative_degree_x == 0 .AND. derivative_degree_y == 1)THEN
	r = C3 + C4 * x
ELSEIF(derivative_degree_x == 1 .AND. derivative_degree_y == 1)THEN
	r = C4
ENDIF

END SUBROUTINE
	
