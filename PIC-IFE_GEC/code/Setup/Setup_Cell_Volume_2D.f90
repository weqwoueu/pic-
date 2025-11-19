!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: Setup_Cell_Volume.f90                              C
!
!  Purpose: Get cell volume for every grid point, with cell volume to get rho density on every point
!                                                                      C
!  Author: Yuchuan Chu                                Date: 13-Dec-12  C
!  Comments:                                                           C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
SUBROUTINE Setup_Cell_Volume_2D(delta, xmin, xmax, ymin, ymax, nnx, nny, N_Objects, N_Boundary, objects)

USE IFE_MAIN_PARAM
USE IFE_Data
USE IFE_Boundary
USE Object_Data_2D
USE IFE_INTERFACE
USE Cell_Volume_Data
USE IFE_Interface_cell_volume, ONLY:Cubic_Partition_Cell_Volume_2D, Mesh_Objects_Intersection_info_2D_BL_cell_volume, &
                                    Calculate_Cell_Volume_2D

IMPLICIT NONE


INTEGER		                                ::  delta 
REAL(8)		                                ::  xmin, xmax, ymin, ymax
INTEGER		                                ::  nnx, nny, N_Objects, N_Boundary
TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
INTEGER                                        ::  num_of_nodes, n_int_elements, i, &
                                                    j, FID, New_IFE_Mesh_Option, num_of_elements
REAL(8)		                                ::  dimensions(2,2)
INTEGER		                                ::  nodes(2)
INTEGER, DIMENSION(:), POINTER              ::  node_index_cell_volume
INTEGER, DIMENSION(:,:), POINTER            ::  t_c_cell_volume
REAL(8), DIMENSION(:,:), POINTER            ::  p_basic_cell_volume
INTEGER, DIMENSION(:), POINTER              ::  element_index_cell_volume
INTEGER,DIMENSION(:,:),POINTER              ::  information_1_cell_volume
REAL(8),DIMENSION(:,:),POINTER              ::  information_2_cell_volume


FID = 1
WRITE(6,*)
WRITE(6,*)'READ "ife.inp"'
WRITE(6,*)
OPEN(FID, ACTION = 'READ', FILE = 'ife.inp')
	READ(FID,*) !New_IFE_Mesh_Option
	READ(FID,*)	!Int_El_Frac
	READ(FID,*)	!NSolver
	READ(FID,*)	!SLSolver
	READ(FID,*)	!Blocking
	READ(FID,*)	!PCG_MaxIT
	READ(FID,*)	!PCG_Tol
	READ(FID,*)	!BGS_MaxIT
	READ(FID,*)	!BGS_Tol
	READ(FID,*)	!Newton_MaxIT			  
	READ(FID,*)	!Newton_Tol
	READ(FID,*)	delta
CLOSE(FID)

IF(delta ==0) THEN
	WRITE(6,*)'delta = 0 : 2-D Cartesian coordinates for cell volume' 
ELSEIF(delta ==1) THEN
	WRITE(6,*)'delta = 1: cylindrical coordinate, axisymmetric situation for cell volume'
ELSE
	WRITE(6,*)'delta != 0 and delta != 1, WRONG delta setup, STOP'
	STOP
ENDIF
WRITE(6,*)

WRITE(6,*)
WRITE(6,*) 'Setup IFE Mesh'
WRITE(6,*) '============= '


dimensions(:,1) = (/xmin, ymin/)
dimensions(:,2) = (/xmax, ymax/)
nodes = (/nnx+1, nny+1/)

!h_partition(1) = ( dimensions(1,2)- dimensions(1,1))/(2*(nnx - 1))
!h_partition(2) = (dimensions(2,2) - dimensions(2,1))/(2*(nny - 1))

CALL Cubic_Partition_Cell_Volume_2D(dimensions, nodes, p_basic_cell_volume, t_c_cell_volume)

num_of_nodes	= SIZE(p_basic_cell_volume,2)

ALLOCATE(node_index_cell_volume(num_of_nodes))
node_index_cell_volume = -1
num_of_elements = SIZE(t_c_cell_volume,2)
ALLOCATE(element_index_cell_volume(num_of_elements))
element_index_cell_volume = -1



IF( N_Objects /= 0) THEN

    DO i = 1, num_of_nodes
	   DO j = 1, N_Objects
	   	  	IF( node_index_cell_volume(i)==-1) THEN
		      CALL Locate_Node_2D (p_basic_cell_volume(:,i), objects(j), node_index_cell_volume(i))
			ENDIF
	   END DO
    END DO 

    ALLOCATE(information_1_cell_volume(18,SIZE(t_c_cell_volume,2)), information_2_cell_volume(8,SIZE(t_c_cell_volume,2)))
	CALL Mesh_Objects_Intersection_info_2D_BL_cell_volume(t_c_cell_volume,	p_basic_cell_volume, 	&
											objects, Int_El_Frac,	&
											element_index_cell_volume, information_1_cell_volume,   &
											information_2_cell_volume, node_index_cell_volume)

ENDIF

CALL Calculate_Cell_Volume_2D(t_c_cell_volume, p_basic_cell_volume, element_index_cell_volume, &
        information_1_cell_volume, information_2_cell_volume, node_index_cell_volume, cell_volume)

!IF( N_Objects /= 0) THEN
!   DEALLOCATE(objects)
!ENDIF

END SUBROUTINE