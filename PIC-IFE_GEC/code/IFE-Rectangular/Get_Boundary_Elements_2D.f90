SUBROUTINE Get_Boundary_Elements_2D(n_elements, t_basic, node_type, node_loc, bnd_elem_index)


! Purpose:		Get boundary elements; either on outer or inner(immersed) boundaries
! Last Update:	10/6/2003 04:15 AM

USE Object_Data_2D

IMPLICIT NONE

INTEGER, DIMENSION(:,:), POINTER		::	t_basic, node_type
INTEGER, DIMENSION(:), POINTER			::	node_loc
INTEGER, DIMENSION(:), POINTER			::	bnd_elem_index

INTEGER								::	m_t, n_t, bc_elem_count, e, i, n_elements
INTEGER, DIMENSION(:), ALLOCATABLE	::	bnd_elem_index_tmp

m_t			= 4	  ! node number for a elemnet
n_t			= n_elements ! element number
	
ALLOCATE(bnd_elem_index_tmp(n_t))

bc_elem_count = 0


DO e=1,n_t	
 
	DO i=1,m_t

		IF ( node_type(1,t_basic(i,e))>=0 ) THEN	! Boundary Drichilet or Neumann node
			bc_elem_count = bc_elem_count + 1
			bnd_elem_index_tmp(bc_elem_count) = e

			EXIT
!!Y.C delete node_loc<-1 is not a fixed potential node, delete here, NON-vacuum point is treated as float point 
		ELSEIF(size(node_loc) /= 0 ) THEN

			IF ( node_loc(t_basic(i,e))<vacuum ) THEN	! Fixed potential node ( Interior point)
				bc_elem_count = bc_elem_count + 1
				bnd_elem_index_tmp(bc_elem_count) = e

				EXIT
			ENDIF
!!Y.C delete node_loc<-1 is not a fixed potential node, delete here, NON-vacuum point is treated as float point 
		ENDIF
	END DO

END DO

ALLOCATE(bnd_elem_index(bc_elem_count))
bnd_elem_index = bnd_elem_index_tmp(1:bc_elem_count)
DEALLOCATE(bnd_elem_index_tmp)


END