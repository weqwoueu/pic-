SUBROUTINE Partition_Cube_2D( p_basic, e_basic, dimensions, 	&
							node_type, num_of_unknowns)
USE IFE_Boundary
IMPLICIT NONE

REAL(8), DIMENSION(:,:), POINTER		::	p_basic
INTEGER, DIMENSION(:,:), POINTER		::	e_basic
REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions
INTEGER, DIMENSION(:,:), POINTER		::	node_type
INTEGER, INTENT(OUT)					::	num_of_unknowns

INTEGER									::	n_nodes, k, unknown_count, n_edges, i
REAL(8)									::	x_min, x_max, y_min, y_max
REAL(8)									::	vert(2), tol

x_min = dimensions(1,1);
x_max = dimensions(1,2);
y_min = dimensions(2,1);
y_max = dimensions(2,2);

n_nodes		= SIZE(p_basic,2);
n_edges		= SIZE(e_basic,2);

tol = EPSILON(1.);

unknown_count = 0;

ALLOCATE (node_type(2,n_nodes))

! node_type(1,:) = 1 Drichilet, 0 Neumann, -1 Interion #1, -2 Interior #2, ... etc
! node_type(2,:) = 0 Known node, +ve Unknown node

node_type = -200
node_type(1,68)=1   !四边诺伊曼，内部一点为迪利克雷（自己选定）
node_type(2,68)=-1

DO k = 1,n_nodes
	
		IF ( node_type(2,k) == -200 ) THEN ! this node has not been indexed
!			vert = p_basic(:,k);

!			IF (ABS(vert(1) - x_min) < tol) THEN ! this is node is on the left surface
!				IF (ABS(vert(2) - y_min) < tol) THEN  ! this node is also on right edge of the left surface
!					node_type(1,k) = MAX(bc_index(1), bc_index(3));
!				ELSEIF (ABS(vert(2) - y_max) < tol) THEN !this node is also on left edge of the left surface	
!					node_type(1,k) = MAX(bc_index(1), bc_index(4));
!				ELSE !this node is in the interior of the left surface
!					node_type(1,k) = bc_index(1);
!				END	IF			
!
!			ELSEIF (ABS(vert(1) - x_max) < tol) THEN ! this is node is on the right surface
!				IF (ABS(vert(2) - y_min) < tol) THEN  ! this node is also on left edge of the right surface
!					node_type(1,k) = MAX(bc_index(2), bc_index(3));
!				ELSEIF (ABS(vert(2) - y_max) < tol) THEN !this node is also on right edge of the right surface	
!					node_type(1,k) = MAX(bc_index(2), bc_index(4));
!				ELSE !this node is in the interior of the right surface
!					node_type(1,k) = bc_index(2);
!				END	IF		
!
!			ELSEIF (ABS(vert(2) - y_min) < tol) THEN ! this is node is on the front surface
!				IF (ABS(vert(1) - x_min) < tol) THEN  ! this node is also on left edge of the front surface
!					node_type(1,k) = MAX(bc_index(3), bc_index(1));
!				ELSEIF (ABS(vert(1) - x_max) < tol) THEN !this node is also on right edge of the front surface	
!					node_type(1,k) = MAX(bc_index(3), bc_index(2));
!				ELSE !this node is in the interior of the front surface
!					node_type(1,k) = bc_index(3);
!				END	IF
!			
!			ELSEIF (ABS(vert(2) - y_max) < tol) THEN ! this is node is on the back surface
!				IF (ABS(vert(1) - x_min) < tol) THEN  ! this node is also on right edge of the back surface
!					node_type(1,k) = MAX(bc_index(4), bc_index(1));
!				ELSEIF (ABS(vert(1) - x_max) < tol) THEN !this node is also on left edge of the back surface	
!					node_type(1,k) = MAX(bc_index(4), bc_index(2));
!				ELSE !this node is in the interior of the back surface
!					node_type(1,k) = bc_index(4);
!				END	IF	
!		
!			END IF

			DO i = 1, n_edges
			   IF (k == e_basic(1, i)) THEN
				 node_type(1,k) = MAX(bc_index(e_basic(4, i)), node_type(1,k))
			   ENDIF

			   IF (k == e_basic(2, i)) THEN
				 node_type(1,k) = MAX(bc_index(e_basic(5, i)), node_type(1,k))
			   ENDIF

			ENDDO 

			IF (node_type(1,k) < 1) THEN			
			! This is an INTERIOR node or a NEUMANN BOUNDARY node
            	unknown_count = unknown_count + 1;
            	node_type(2,k) = unknown_count;
            ELSE										
			! This is a DRICHILET BOUNDARY node
            	node_type(2,k) = -node_type(1,k); ! = -1
            END IF

		END IF
END	DO

num_of_unknowns = unknown_count


END