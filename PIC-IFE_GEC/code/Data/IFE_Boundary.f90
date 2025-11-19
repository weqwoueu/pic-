!					==================================
!					||      IFE_Boundary_Interface      ||
!					==================================


MODULE IFE_Boundary

USE IFE_MAIN_PARAM

!INTEGER, DIMENSION(:,:), POINTER	::	e_basic
REAL(8), DIMENSION(:,:), POINTER	::	bc_point_1, bc_point_2
INTEGER, DIMENSION(:), POINTER	    ::	bc_index
!INTEGER, DIMENSION(:), POINTER		::	bc_value
REAL(8), DIMENSION(:), POINTER		::	bc_value

END MODULE