SUBROUTINE Userdef_EBC_Value(p_basic, t_basic_int, node_type, EBC_Value)


USE IFE_MAIN_PARAM

IMPLICIT NONE

REAL(8), DIMENSION(:,:), INTENT(IN)			::	p_basic
INTEGER, DIMENSION(:,:), INTENT(IN)			::	t_basic_int, node_type
REAL(8), DIMENSION(:,:), INTENT(OUT)		::	EBC_Value

INTEGER		i, j, k ,e ,num_point_elem, node_index, num_elem
REAL(8)		x, y, xmin, ymin, xmax, ymax, r1, r2, beta

num_point_elem = 4
num_elem = SIZE(t_basic_int,2)

DO i=1,4
	Do j=1,num_elem
		EBC_Value(i,j) = Zero
	ENDDO
ENDDO

r1=19.9
r2=89.9
beta = 1

DO e=1, num_elem
	DO i=1, num_point_elem
		node_index = t_basic_int(i,e)
		IF(node_type(1,node_index) == 1)THEN
				xmin = MINVAL(p_basic(1,:))
				xmax = MAXVAL(p_basic(1,:))
				ymin = MINVAL(p_basic(2,:))
				ymax = MAXVAL(p_basic(2,:))
                x = p_basic(1, node_index)
                y = p_basic(2, node_index)
                IF(p_basic(1, node_index) == xmin) THEN
                    IF (y<r1) THEN
                        EBC_Value(i,e) = -5!(x*x+y*y)/beta+(1-1/beta)*r1*r1 
                    ELSEIF (y>r1 .AND. y<r2) THEN
                        EBC_Value(i,e) = x*x+y*y
                    ELSEIF (y>r2) THEN
                        EBC_Value(i,e) = -5!(x*x+y*y)/beta+(1-1/beta)*r2*r2 
                    ENDIF
                    
                ELSEIF(p_basic(1, node_index) == xmax) THEN
                    EBC_Value(i,e) = -5!(x*x+y*y)/beta+(1-1/beta)*r2*r2 
                ELSEIF((p_basic(2, node_index) == ymin)) THEN
				    IF (x<r1) THEN
                        EBC_Value(i,e) = -5!(x*x+y*y)/beta+(1-1/beta)*r1*r1 
                    ELSEIF (x>r1 .AND. x<r2) THEN
                        EBC_Value(i,e) = x*x+y*y
                    ELSEIF (x>r2) THEN
                        EBC_Value(i,e) = -5!(x*x+y*y)/beta+(1-1/beta)*r2*r2 
                    ENDIF	
                    
				ELSEIF((p_basic(2, node_index) == ymax) ) THEN
					EBC_Value(i,e) = -5!(x*x+y*y)/beta+(1-1/beta)*r2*r2			
				ENDIF
                
                
!				IF( (p_basic(1, node_index) == xmin) .OR. (p_basic(1, node_index) == xmax) ) THEN
!					y = p_basic(2, node_index)
!!					EBC_Value(i,e) = (1 + y**2)**2.5
!					EBC_Value(i,e) = y
!				ELSEIF((p_basic(2, node_index) == ymin)) THEN
!				    EBC_Value(i,e) = 0.0487805			
!				ELSEIF((p_basic(2, node_index) == ymax) ) THEN
!					y = p_basic(2, node_index)
!!					EBC_Value(i,e) = (x**2 + 1)**2.5
!					EBC_Value(i,e) = y				
!				ENDIF
		ENDIF
	ENDDO
ENDDO



END SUBROUTINE