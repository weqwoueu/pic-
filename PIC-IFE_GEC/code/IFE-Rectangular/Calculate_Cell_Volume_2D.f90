SUBROUTINE Calculate_Cell_Volume_2D(t_basic, p_basic, element_index, information_1, information_2, node_index, cell_volume_temp)

IMPLICIT NONE

INTEGER, DIMENSION(:,:), INTENT(IN)		::	t_basic
REAL(8), DIMENSION(:,:), INTENT(IN)		::	p_basic
INTEGER, DIMENSION(:), INTENT(IN)		::	element_index
INTEGER, DIMENSION(:,:),INTENT(IN)		::  information_1
REAL(8), DIMENSION(:,:),INTENT(IN)		::  information_2
INTEGER, DIMENSION(:), INTENT(IN)		::	node_index
REAL(8), DIMENSION(:), POINTER          ::  cell_volume_temp

INTEGER                                    ::  num_of_element
INTEGER                                    ::  i,j
REAL(8), DIMENSION(2,4)                  ::  vert
REAL(8)                                    ::  area
INTEGER, DIMENSION(4)                    ::  node_index_el
INTEGER                                    ::  n

num_of_element = SIZE(t_basic,2)
n = 0

ALLOCATE(cell_volume_temp(num_of_element))

DO i=1, num_of_element
    DO j=1,4
        vert(:,j) = p_basic(:,t_basic(j,i))
        node_index_el(j) = node_index(t_basic(j,i))
    ENDDO
    IF(element_index(i) <  0) THEN             ! non-interface element
        CALL Cell_Volume_FE_2D(vert, area)
    ELSEIF(element_index(i) > 0) THEN         ! interface element
        n = n+1
        CALL Cell_Volume_IFE_2D(vert, information_1(:,n), information_2(:,n), node_index_el, area)
    ENDIF
    
    cell_volume_temp(i) = area
ENDDO


END SUBROUTINE