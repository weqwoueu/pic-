SUBROUTINE Gauss_Nodes_2D(p_basic, t_basic, g_x, g_y)

IMPLICIT NONE

REAL(8), DIMENSION(:,:), INTENT(IN)		::	p_basic
INTEGER, DIMENSION(:,:), INTENT(IN)		::	t_basic
REAL(8), DIMENSION(:,:), POINTER		::	g_x, g_y

INTEGER									::	n_gnodes, n_t, k
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	vert, gnodes

n_gnodes = 3

n_t		 = SIZE(t_basic,2)

ALLOCATE(g_x(n_gnodes, n_t), g_y(n_gnodes, n_t))

ALLOCATE(vert(2,3), gnodes(2,3))

DO k = 1,n_t
	
	vert = p_basic(:,t_basic(:,k))
	
	CALL Gauss_Nodes_Elem_2D(vert, gnodes)
		
	g_x(:,k) = gnodes(1,:)
	g_y(:,k) = gnodes(2,:)
	
END DO

DEALLOCATE(vert, gnodes)

END
