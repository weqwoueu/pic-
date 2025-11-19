SUBROUTINE Setup_IFE_Mesh_2D(delta, xmin, xmax, ymin, ymax, nnx, nny, N_Objects, N_Boundary, objects)

! Purpose:		Setup IFE Mesh, do objects-mesh intersections and save mesh information to a IFE_Mesh_Filename.
! Last Updated:	4/21/2004 01:21 PM

USE IFE_MAIN_PARAM
USE IFE_Data
USE IFE_Boundary
USE Object_Data_2D
USE IFE_INTERFACE, ONLY: Cubic_Partition_2D, Boundary_Condition_info_2D, &
                        modify_interface_elements, Mesh_Objects_Intersection_info_2D_BL, &
                        Update_Mesh_Objects_Intersection_2D, Update_Locate_Node_2D, &
                        generate_Adaptive_DGT_IFE, generate_HP_HT_2D, generate_HE_2D, &
                        generate_DGP_DGT_2D, generate_DGE_2D, generate_Adaptive_DG_P_T_E, &
                        generate_Adaptive_DG_P_T_E_repeat, generate_Adaptive_DG_P_average, &
                        generate_E_Dirichlet_2D, generate_DG_edge_flag, &
                        Update_Edge_Objects_Intersection_2D, generate_DG_node_flag, &
                        generate_Tier_initial_2D, Periodic_pair_boundary_nodes, Improve_HT, &
                        generate_Elementmap

IMPLICIT NONE



INTEGER		delta 
!===================new add===============================
INTEGER		LOG_modify
!=========================================================

!==================Add for surface jump===================
INTEGER, DIMENSION(:),   POINTER	    ::    modified_element_index
INTEGER, DIMENSION(:,:), POINTER      ::    modified_information_1
REAL(8), DIMENSION(:,:), POINTER      ::    modified_information_2
!=========================================================

REAL(8)		xmin, xmax, ymin, ymax   !, PIC_Phi_BC(6)
INTEGER		nnx, nny, N_Objects, N_Boundary


INTEGER		 num_of_nodes, n_int_elements, i, j, FID, New_IFE_Mesh_Option, num_of_elements
REAL(8)		dimensions(2,2)
INTEGER		nodes(2)

!--------- DG------------------
INTEGER   :: area_flag
Real(8), Allocatable :: node_map_temp(:,:)
Real(8)              :: h_min(2)
!------------------------------

TYPE(ObjectType), DIMENSION(:), POINTER		::	objects

REAL(8)                                         beta_minus, beta_plus, x, y
REAL(8), DIMENSION(:), POINTER				::	U_fe_full

Integer(4) :: nindex, k
!==========LY modification for Speeding Particle Positioning, 2022-7-25==========
INTEGER :: N_left_bottom_child, N_left_top_child, N_right_bottom_child, N_right_top_child
Integer :: MaxArrayNumber, EdgeNumber, count, ChildElementNumber
!==========LY modification for Speeding Particle Positioning, 2022-7-25=========


FID = 1
WRITE(6,*)
WRITE(6,*)'READ "ife.inp"'
WRITE(6,*)
OPEN(FID, ACTION = 'READ', FILE = './INPUT/ife.inp')
!	READ(FID,*) IFE_Input_Filename
!	READ(FID,*) IFE_Mesh_Filename
	READ(FID,*) New_IFE_Mesh_Option
	READ(FID,*)	Int_El_Frac
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
	WRITE(6,*)'delta = 0 : 2-D Cartesian coordinates' 
ELSEIF(delta ==1) THEN
	WRITE(6,*)'delta = 1: cylindrical coordinate, axisymmetric situation'
ELSE
	WRITE(6,*)'delta != 0 and delta != 1, WRONG delta setup, STOP'
	STOP
ENDIF
WRITE(6,*)


IF(New_IFE_Mesh_Option ==1 ) THEN
	WRITE(6,*)'Setup and save a new IFE mesh'
	WRITE(6,*)'Interface elements to total elements fraction (initial guess)='
	WRITE(6,*)Int_El_Frac
ENDIF

IF (New_IFE_Mesh_Option/=1) THEN
	WRITE(6,*)'Load from a saved IFE mesh, RETURN'
	RETURN	
!	STOP
END IF

!	SETUP NEW IFE MESH
WRITE(6,*)
WRITE(6,*) 'Setup IFE Mesh'
WRITE(6,*) '============= '


dimensions(:,1) = (/xmin, ymin/)
dimensions(:,2) = (/xmax, ymax/)
nodes = (/nnx, nny/)

h_partition(1) = ( dimensions(1,2)- dimensions(1,1))/(nodes(1) - 1)
h_partition(2) = (dimensions(2,2) - dimensions(2,1))/(nodes(2) - 1)



WRITE(6,*) 'x(z)min=',dimensions(1,1),' x(z)max=',dimensions(1,2)
WRITE(6,*) 'y(r)min=',dimensions(2,1),' y(r)max=',dimensions(2,2)
WRITE(6,*) 'nodes_x(z)=',nodes(1),' nodes_y(r)=',nodes(2)
WRITE(6,*) 'Number of Objects  =', N_Objects
WRITE(6,*) 'Number of Boudaries=', N_Boundary
WRITE(6,*)

WRITE(6,*) 'Mesh Generation ....'

CALL Cubic_Partition_2D(dimensions, nodes, p_basic, t_c)

OPEN(1,ACTION='READ',FILE='./INPUT/IDG_inf.inp')
  READ(1,*)
  READ(1, *) repeat_refinement
  READ(1, *) area_flag
CLOSE(1)


!-----------------------------------SIDG------------------------------------
CALL generate_Adaptive_DGT_IFE(p_basic, t_c, objects, DGT_IFE_partition, area_flag)
CALL generate_HP_HT_2D(p_basic, t_c, DGT_IFE_partition, HP, HT, P_average)
CALL generate_HE_2D(HP, HT, nny - 1, HE)

!==========LY Add for Multigrid Store, 2022-1-17==========
CALL generate_Tier_initial_2D(HP, HT, CellMesh)
!==========LY Add for Multigrid Store, 2022-1-17==========

IF (repeat_refinement == 0) THEN
  CellMesh(1:SIZE(CellMesh))%isSplitted = 0
END IF
WRITE(6,*) 'SIDG information Generation Done.'
!-----------------------------------SIDG------------------------------------

!----------------------------------REPEAT MESH------------------------------
IF (repeat_refinement > 0) THEN

  !LY REVISE, 2021-12-7
  !If case need to mesh refinement, then we need to save every mesh information about HP, HT, HT_flag after refinement.
  !when repeat_refinement great than 0, then we firstly need to save the 1st mesh information.
  
  !==========First mesh refinement==========
  !CALL generate_Adaptive_DGT_IFE(HP, HT, objects, DGT_IFE_partition, area_flag)
  CALL generate_Adaptive_DG_P_T_E(HP, HT, HE, DGT_IFE_partition, AHP, AHT, AHE, HT_flag, &
                                  CellMesh, CellMesh_Temp, N_Tier_Old)
  
  DEALLOCATE(HP, HT, HE)
  ALLOCATE(HP(2, SIZE(AHP, 2)))
  ALLOCATE(HT(6, SIZE(AHT, 2)))
  ALLOCATE(HE(8, SIZE(AHE, 2)))
  HP(1:2, 1:SIZE(AHP,2)) = 0.0
  HT(1:6, 1:SIZE(AHT,2)) = 0
  HE(1:8, 1:SIZE(AHE,2)) = 0
  
  DO i = 1, SIZE(AHP, 2)
    HP(1:2, i) = AHP(1:2, i)
  END DO
  
  DO i = 1, SIZE(AHT, 2)
    HT(1:6, i) = AHT(1:6, i)
  END DO
  
  DO i = 1, SIZE(AHE, 2)
    HE(1:8, i) = AHE(1:8, i)
  END DO
  !==========First mesh refinement==========
  
  !==========LY Add for Multigrid Store, 2022-1-17==========
  DEALLOCATE(CellMesh)
  ALLOCATE(CellMesh(SIZE(CellMesh_Temp, 1)))
  DO i = 1, SIZE(CellMesh_Temp, 1)
    CellMesh(i) = CellMesh_Temp(i)
  END DO
  !==========LY Add for Multigrid Store, 2022-1-17==========
  
  !LY REVISE, 2021-12-7
  !when finish the first mesh refinement, 
  !we need to construct the connection between 1st mesh information(no refinement) and 2st mesh information(first refinement).

  
  !LY REVISE, 2021-12-7
  !when finish the first mesh refinement, we firstly need to protect the 1st mesh information, then to store the 2st mesh information.
  
  
  DO k = 2, repeat_refinement
  
    CALL generate_Adaptive_DGT_IFE(HP, HT, objects, DGT_IFE_partition, area_flag)
    CALL generate_Adaptive_DG_P_T_E_repeat(HP, HT, HE, DGT_IFE_partition, k, AHP, AHT, AHE, HT_flag, &
                                            CellMesh, CellMesh_Temp, N_Tier_Old)
    
    DEALLOCATE(HP, HT, HE)
    ALLOCATE(HP(2, SIZE(AHP, 2)))
    ALLOCATE(HT(6, SIZE(AHT, 2)))
    ALLOCATE(HE(8, SIZE(AHE, 2)))
    HP(1:2, 1:SIZE(AHP,2)) = 0.0
    HT(1:6, 1:SIZE(AHT,2)) = 0
    HE(1:8, 1:SIZE(AHE,2)) = 0
    
    !==========LY Add for Multigrid Store, 2022-1-18==========
    DEALLOCATE(CellMesh)
    ALLOCATE(CellMesh(SIZE(CellMesh_Temp, 1)))
    DO i = 1, SIZE(CellMesh_Temp, 1)
      CellMesh(i) = CellMesh_Temp(i)
    END DO
    !==========LY Add for Multigrid Store, 2022-1-18==========
    
    DO i = 1, SIZE(AHP, 2)
      HP(1:2, i) = AHP(1:2, i)
    END DO

    DO i = 1, SIZE(AHT, 2)
      HT(1:6, i) = AHT(1:6, i)
    END DO

    DO i = 1, SIZE(AHE, 2)
      HE(1:8, i) = AHE(1:8, i)
    END DO
    
    
  END DO
    
  DEALLOCATE(AHP, AHT, AHE)
  
  CALL generate_Adaptive_DG_P_average(HP, P_average)
  
  Write(*,*) 'The mesh refinement Done.'
  
END IF

CALL generate_DG_node_flag(HP, HT, P_average, P_flag)

IF (repeat_refinement > 0) THEN
    CALL Improve_HT(HT, HP, HE, repeat_refinement) ! the key of ImproveSIDG
END IF

!CALL generate_DG_edge_flag(DGE, dimensions, DGP, DGT, DG_flag)
!----------------------------------REPEAT MESH------------------------------

!---wsy revise at here---
!------------WSY add for diagnosis (x-axis average) ----------------------
count = 0
Allocate(node_map_temp(2, nnx*  (2**(repeat_refinement)) ))
node_map_temp = 0
h_min(:) = h_partition(:)/  (2**(repeat_refinement))
Do i = 1, size(HP,2)
    k = ((HP(1,i)-xmin)/h_min(1)) + 1 
    node_map_temp(1,k) = node_map_temp(1,k) + 1
    node_map_temp(2,k) = HP(1,i)
End DO

Do i = 1, size(node_map_temp, 2)
    If (node_map_temp(1,i) /= 0) Then
        count = count + 1
    End If
End Do

Allocate(node_map(2,count))
count = 0
Do i = 1, size(node_map_temp, 2)
    If (node_map_temp(1,i) /= 0) Then
        count = count + 1
        node_map(1,count) = i
        node_map(2,count) = node_map_temp(2,i)
    End If
End Do
Deallocate(node_map_temp)
!--------------------------------------------------------



!=========LY modification for Speeding Particle Positioning, 2022-7-25=========
Do i = 1, Size(t_c,2)
  If (CellMesh(i)%isSplitted == 1) Then
    N_left_bottom_child   = CellMesh(i)%Child(1)
    N_left_top_child      = CellMesh(i)%Child(2)
    N_right_bottom_child  = CellMesh(i)%Child(3)
    N_right_top_child     = CellMesh(i)%Child(4)

    !Searching the left-bottom child element until finding final element index.
    Do While (CellMesh(N_left_bottom_child)%Finalindex == 0)
      N_left_bottom_child = CellMesh(N_left_bottom_child)%Child(1)
    End Do
    CellMesh(i)%Nodeindex(1) = CellMesh(N_left_bottom_child)%Nodeindex(1)

    !Searching the left-top child element until finding final element index.
    Do While (CellMesh(N_left_top_child)%Finalindex == 0)
      N_left_top_child = CellMesh(N_left_top_child)%Child(2)
    End Do
    CellMesh(i)%Nodeindex(4) = CellMesh(N_left_top_child)%Nodeindex(4)

    !Searching the right-bottom child element until finding final element index.
    Do While (CellMesh(N_right_bottom_child)%Finalindex == 0)
      N_right_bottom_child = CellMesh(N_right_bottom_child)%Child(3)
    End Do
    CellMesh(i)%Nodeindex(2) = CellMesh(N_right_bottom_child)%Nodeindex(2)

    !Searching the right-top child element until finding final element index.
    Do While (CellMesh(N_right_top_child)%Finalindex == 0)
      N_right_top_child = CellMesh(N_right_top_child)%Child(4)
    End Do
    CellMesh(i)%Nodeindex(3) = CellMesh(N_right_top_child)%Nodeindex(3)
  End If
End Do

!Firstly, We need to construct the relationship between initial and no-refined element and it's valid child element
Allocate(ElementParentChildTemp(100*(repeat_refinement+1), Size(t_c,2)))  !Enough big
ElementParentChildTemp(1:100*(repeat_refinement+1), 1:Size(t_c,2)) = 0
Do i = 1, Size(CellMesh)
  If (CellMesh(i)%Finalindex /= 0) Then
    ElementParentChildTemp(1,CellMesh(i)%FinalParent) = ElementParentChildTemp(1,CellMesh(i)%FinalParent) + 1
    ElementParentChildTemp(1+ElementParentChildTemp(1,CellMesh(i)%FinalParent),CellMesh(i)%FinalParent) = CellMesh(i)%Finalindex
  End If
End Do

MaxArrayNumber = 1 + MAXVAL(ElementParentChildTemp(1, 1:Size(t_c,2)))
Allocate(ElementParentChild(MaxArrayNumber, Size(t_c,2)))
ElementParentChild(1:MaxArrayNumber, 1:Size(t_c,2)) = 0
Do i = 1, Size(t_c,2)
  ElementParentChild(1:MaxArrayNumber, i) = ElementParentChildTemp(1:MaxArrayNumber, i)
End Do

!Secondly, We need to construct the relationship between DG edge and it's element
Allocate(EdgeElementTemp(50, Size(HT,2)))  !Enough big
EdgeElementTemp(1:50, 1:Size(HT,2)) = 0
Do i = 1, Size(HE,2)
  EdgeElementTemp(1, HE(5,i)) = EdgeElementTemp(1, HE(5,i)) + 1
  EdgeElementTemp(1+EdgeElementTemp(1,HE(5,i)), HE(5,i)) = i
End Do

MaxArrayNumber = 1 + MAXVAL(EdgeElementTemp(1, 1:Size(HT,2)))
Allocate(EdgeElement(MaxArrayNumber, Size(HT,2)))
EdgeElement(1:MaxArrayNumber, 1:Size(HT,2)) = 0
Do i = 1, Size(HT,2)
  EdgeElement(1:MaxArrayNumber, i) = EdgeElementTemp(1:MaxArrayNumber, i)
End Do

Deallocate(ElementParentChildTemp, EdgeElementTemp)

!Thirdly, We need to construct the relationship between initial and no-refined element and it's DG edge.
Allocate(EdgeParentNumber(Size(t_c,2)))
EdgeParentNumber = 0
Do i = 1, Size(t_c,2)
  EdgeNumber = 0
  ChildElementNumber = ElementParentChild(1,i)
  Do j = 2, 1+ChildElementNumber
    EdgeNumber = EdgeNumber + EdgeElement(1,ElementParentChild(j,i))
  End Do
  EdgeParentNumber(i) = EdgeNumber
End Do
!The EdgeParentNumber is to store the number of DG edge in any initial and no-refined element.

MaxArrayNumber = 1 + MAXVAL(EdgeParentNumber)
Allocate(EdgeParent(MaxArrayNumber, Size(t_c,2)))
EdgeParent = 0
Do i = 1, Size(t_c,2)
  EdgeParent(1,i) = EdgeParentNumber(i)
  count = 2
  Do j = 2, 1+ElementParentChild(1,i)
    Do k = 2, 1+EdgeElement(1,ElementParentChild(j,i))
      EdgeParent(count,i) = EdgeElement(k,ElementParentChild(j,i))
      count = count + 1
    End Do
  End Do
End Do

Deallocate(EdgeElement, ElementParentChild)
!=========LY modification for Speeding Particle Positioning, 2022-7-25=========

!=========LY modification for Periodic Boundary Condition, 2022-6-16=========
If (bc_index(1)==-1 .OR. bc_index(2)==-1 .OR. bc_index(3)==-1 .OR. bc_index(4)==-1) Then
  !There denote that Periodic boundary is exit.
  Call Periodic_pair_boundary_nodes(nnx, nny, bc_index, HP, bc_point_1, bc_point_2, PairNodes)
  !The PairNodes is to store the Periodic boundary nodes index information.
End If
!=========LY modification for Periodic Boundary Condition, 2022-6-16=========

num_of_nodes	= SIZE(HP,2)

!LY Add for SIDG-PIC Method: Phi, Rho, Rho_s size. 2022-1-14
Field_Size = SIZE(HP, 2)
!LY Add for SIDG-PIC Method: Phi, Rho, Rho_s size. 2022-1-14

ALLOCATE(node_index(num_of_nodes))
node_index = -1

num_of_elements = SIZE(HT,2)
ALLOCATE(element_index(num_of_elements))
element_index = -1

ALLOCATE(U_fe_full(num_of_nodes))
U_fe_full = Zero

WRITE(6,*) 'Mesh Generation Done.'

Call generate_Elementmap(Element_Map, repeat_refinement)

IF( N_Objects /= 0) THEN

  DO i = 1, num_of_nodes
    DO j = 1, N_Objects
      IF( node_index(i)==-1) THEN
        CALL Update_Locate_Node_2D (HP(:,i), objects(j), node_index(i),U_fe_full(i))
      ENDIF
    END DO
  END DO

  ALLOCATE(information_1(18,SIZE(HT,2)), information_2(8,SIZE(HT,2)))
  WRITE(6,*) 'Mesh_Object Intersection ....'
  CALL Update_Mesh_Objects_Intersection_2D  (	HT, HP, objects, &
                                              element_index, information_1, information_2, node_index)

  WRITE(6,*) 'Mesh_Object Intersection Done.'
    
!========================CYC ADD for filter======================================================
OPEN(1, ACTION='READ', FILE='./INPUT/ife.inp')
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
	READ(FID,*)	!delta
	READ(FID,*)	LOG_modify
CLOSE(1)

IF(LOG_modify == 1)THEN
	ALLOCATE(modified_element_index(SIZE(element_index,1)))
	ALLOCATE(modified_information_1(SIZE(information_1,1),SIZE(information_1,2)))
	ALLOCATE(modified_information_2(SIZE(information_2,1),SIZE(information_2,2)))

!	ALLOCATE(element_index_for_curve(SIZE(element_index,1)))
!	ALLOCATE(information_1_for_curve(SIZE(information_1,1),SIZE(information_1,2)))
!	ALLOCATE(information_2_for_curve(SIZE(information_2,1),SIZE(information_2,2)))

	CALL modify_interface_elements(h_partition, t_c, p_basic, element_index, information_1, information_2,	&
									modified_element_index, modified_information_1, modified_information_2)
!	element_index_for_curve = element_index
!	information_1_for_curve = information_1
!	information_2_for_curve = information_2
	element_index = 0
	information_1 = 0
	information_2 = 0
	element_index = modified_element_index
	information_1 = modified_information_1
	information_2 = modified_information_2

	DEALLOCATE(modified_element_index, modified_information_1, modified_information_2)
ENDIF


!================================================================================================

    n_int_elements = 0
    DO i = 1, SIZE(element_index,1)
	   IF (element_index(i)>0) THEN	! Interface elemenet: the number refers to the order in t_basic
		   n_int_elements = n_int_elements + 1		
	   END IF
    END DO

    WRITE(6,*) 'No of Mesh Elements      = ',SIZE(HT,2)
    WRITE(6,*) 'No of Mesh Nodes         = ',SIZE(HP,2)
    WRITE(6,*) 'No of Interface Elements = ',n_int_elements
ELSE
    WRITE(6,*) 'No of Mesh Elements      = ',SIZE(HT,2)
    WRITE(6,*) 'No of Mesh Nodes         = ',SIZE(HP,2)
ENDIF

IF (N_Boundary /=0) THEN
   WRITE(6,*)
   WRITE(6,*) 'Boundry Condition Information....'

   CALL Boundary_Condition_info_2D  (SIZE(HT,2), HT, HP, e_basic)
   CALL generate_E_Dirichlet_2D(HP, HT, E_Dirichlet) 
   
   WRITE(6,*) 'Boundry Condition Information Done'
   WRITE(6,*) 'No of Boundary Condition = ',N_Boundary
   WRITE(6,*) 'No of Edge on Boundary   = ',SIZE(e_basic,2)

ENDIF

!-----------------------------LY REVISE, 2021-11-27-------------------------------------
ALLOCATE(edge_index(SIZE(HE,2)))
CALL Update_Edge_Objects_Intersection_2D(HE, HP, HT, information_2, &
                                         node_index, element_index, objects, information_3, edge_index)

CALL Update_Edge_Objects_Intersection_2D(E_Dirichlet, HP, HT, information_2, &
                                         node_index, element_index, objects, information_3_D, edge_index_D)

WRITE(6,*) 'Edge_Object Intersection Done.'
!---------------------------------------------------------------------------------------

WRITE(6,*) 'Saving IFE mesh ...'

OPEN(1, ACTION='WRITE', FILE='ife.msh')
WRITE(1,*) ' Objects material property: num_of_objects'
WRITE(1,*) N_Objects
WRITE(1,*)'! Boundary Condition: num_of_boundary'
WRITE(1,*) N_Boundary

WRITE(1,*) '! Nodes coordinates: num_of_x,num_of_y,num_of_nodes, p_basic'
WRITE(1,*) SIZE(p_basic,2)
DO i=1, SIZE(p_basic,2)
	WRITE(1,*) (p_basic(j,i), j=1,2)
END DO

WRITE(1,*) '! Elements nodal conenctivity: num_of_elements, t_c'
WRITE(1,*) SIZE(t_c,2)
DO i=1,  SIZE(t_c,2)
	WRITE(1,*) (t_c(j,i), j=1,4)
END DO

WRITE(1,*) '! Element location: element_index'
WRITE(1,*) SIZE(element_index)    !LY REVISE, 2021-11-27
DO i=1, SIZE(element_index)
  WRITE(1,*) (element_index(i))
ENDDO

!-------------------LY REVISE, 2021-11-27---------------------------
WRITE(1,*) '! wsy add HP'
WRITE(1,*) SIZE(HP,2)
DO i=1,  SIZE(HP,2)
	WRITE(1,*) (HP(j,i), j=1,2)
END DO

WRITE(1,*) '! wsy add HT'
WRITE(1,*) SIZE(HT,2)
DO i=1,  SIZE(HT,2)
  ! WRITE(*,*) (HT(j,i), j=1,5)
	WRITE(1,*) (HT(j,i), j=1,5)
END DO

WRITE(1,*) '! wsy add HE'
WRITE(1,*) SIZE(HE,2)
DO i=1,  SIZE(HE,2)
	WRITE(1,*) (HE(j,i), j=1,8)
END DO

write(1, *) ' ! wsy add E_Dirichlet'
WRITE(1,*) SIZE(E_Dirichlet,2)
DO i=1,  SIZE(E_Dirichlet,2)
	WRITE(1,*) (E_Dirichlet(j,i), j=1,7)
END DO

write(1, *) ' ! wsy add information_3'
WRITE(1,*) SIZE(information_3,2)
DO i=1,  SIZE(information_3,2)
	WRITE(1,*) (information_3(j,i), j=1,7)
END DO

write(1, *) ' ! wsy add information_3_D'
WRITE(1,*) SIZE(information_3_D,2)
DO i=1,  SIZE(information_3_D,2)
	WRITE(1,*) (information_3_D(j,i), j=1,7)
END DO

write(1, *) ' ! wsy add edge_index'
DO i=1,  SIZE(HE,2)
	WRITE(1,*) (edge_index(i)) 
END DO

write(1, *) ' ! wsy add edge_index_D'
DO i=1,  SIZE(E_Dirichlet,2)
	WRITE(1,*) (edge_index_D(i))
END DO

WRITE(1,*) ' ! LY Add P average,2021-11-24'
WRITE(1,*) SIZE(P_average,2)
DO i=1,  SIZE(P_average,2)
	WRITE(1,*) (P_average(j,i), j=1,6)
END DO

WRITE(1,*) ' ! LY Add P flag,2021-12-3'
WRITE(1,*) SIZE(P_flag,2)
DO i=1,  SIZE(P_flag,2)
	WRITE(1,*) (P_flag(j,i), j=1,5)
END DO
!-------------------------------------------------------------------

IF( N_Objects /= 0) THEN
  WRITE(1,*) '! Interface elements: num_of_interface_elements'
  WRITE(1,*) n_int_elements

  IF (n_int_elements/=0) THEN
    WRITE(1,*) '! Intersection information: '
    DO i=1, n_int_elements

      WRITE(1,*) information_1(:,i)
      WRITE(1,*) information_2(:,i)

    END DO
  ELSEIF (n_int_elements == 0) THEN

  ENDIF

  WRITE(1,*) '! Nodes location index: node_index'
  DO i=1,SIZE(node_index)
    WRITE(1,*) node_index(i)
  END DO

ENDIF

IF( N_Boundary /= 0) THEN
   WRITE(1,*)'! Boundary Condition'
   DO i = 1,  N_Boundary
      WRITE(1, *) bc_index(i), bc_value(i)
      WRITE(1, *) bc_point_1(:, i)
      WRITE(1, *) bc_point_2(:, i)
   ENDDO
   WRITE(1,*) '! Edge Boundary Condition'
   WRITE(1,*) SIZE(e_basic, 2)
   DO i = 1,  SIZE(e_basic, 2)
      WRITE(1, *) e_basic(:,i)
   ENDDO

ENDIF

CLOSE(1)

WRITE(6,*) 'Saving IFE mesh Done.' 

DEALLOCATE(p_basic , t_c)
DEALLOCATE(HE, E_Dirichlet)   !LY, wsy modification: We donot deallcoate HP, HT, P_average and P_flag because SetObject and SetupGrids_2D_QLL. 2022-7-25

IF( N_Objects /= 0) THEN
   !DEALLOCATE(element_index, node_index)  !LY modification, We donot deallcoate node_index because SetObject. 2022-6-13
   DEALLOCATE(element_index)
!   DEALLOCATE(objects)
ENDIF

IF( N_Boundary /= 0) THEN
   DEALLOCATE(bc_index, bc_value, bc_point_1, bc_point_2)
   DEALLOCATE(e_basic)
ENDIF

DEALLOCATE(U_fe_full)

END