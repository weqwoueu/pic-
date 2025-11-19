MODULE IFE_Interface_cell_volume

IMPLICIT NONE

INTERFACE

    SUBROUTINE	Mesh_Objects_Intersection_info_2D_BL_cell_volume(t_basic, p_basic, 	&
													objects, Int_El_Frac,	&
													element_index, information_1, information_2, &
													node_index_el)

		USE Object_Data_2D
		IMPLICIT NONE

		REAL(8),				 INTENT(IN)						::	Int_El_Frac
		INTEGER, DIMENSION(:,:), INTENT(IN)						::	t_basic
		REAL(8), DIMENSION(:,:), INTENT(IN)						::	p_basic
		TYPE(ObjectType), DIMENSION(:), INTENT(IN)				::	objects
		INTEGER, DIMENSION(:), POINTER							::	element_index
		INTEGER, DIMENSION(:,:), POINTER			                ::  information_1
		REAL(8), DIMENSION(:,:), POINTER			                ::  information_2
		INTEGER, DIMENSION(:), INTENT(IN)						::	node_index_el
	END SUBROUTINE


	SUBROUTINE Setup_Cell_Volume_2D(delta, xmin, xmax, ymin, ymax, nnx, nny, N_Objects, N_Boundary, objects)
		USE Object_Data_2D
		INTEGER		delta 
		REAL(8)		xmin, xmax, ymin, ymax
		INTEGER		nnx, nny, N_Objects, N_Boundary
		TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
	END SUBROUTINE
	
    SUBROUTINE Cubic_Partition_Cell_Volume_2D(dimensions, n_nodes, P, T)
		REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions
		INTEGER, DIMENSION(2), INTENT(IN)		::	n_nodes
		REAL(8), DIMENSION(:,:), POINTER 	::	P
		INTEGER, DIMENSION(:,:), POINTER 	::	T
	END SUBROUTINE
	
	SUBROUTINE Calculate_Cell_Volume_2D(t_basic, p_basic, element_index, information_1, information_2, node_index, cell_volume_temp)
        INTEGER, DIMENSION(:,:), INTENT(IN)						::	t_basic
        REAL(8), DIMENSION(:,:), INTENT(IN)						::	p_basic
        INTEGER, DIMENSION(:), INTENT(IN)						::	element_index
        INTEGER, DIMENSION(:,:),INTENT(IN)						::  information_1
        REAL(8), DIMENSION(:,:),INTENT(IN)						::  information_2
        INTEGER, DIMENSION(:), INTENT(IN)						::	node_index
        REAL(8), DIMENSION(:), POINTER                          ::  cell_volume_temp
	END SUBROUTINE
	
END INTERFACE

END MODULE