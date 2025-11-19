SUBROUTINE generate_Adaptive_DG_P_T_E(HP, HT, HE, HT_IFE, AHP, AHT, AHE, HT_flag, &
                                      CellMesh_Old, CellMesh_New, N_Tier_Old)
  
!------------------LY add AHP, AHT, AHE for Adaptive Mesh 2021-9-25------------------
!Note: HP, HT, DGE, HT_IFE--old, not refinement
!      AHP, AHT, AHE--new, refinement
!LY REVISE, 2021-12-3
!AHT(5,:): number of mesh refinement
!AHT(6,:): CG--0, DG--1

!==========LY Add for Multigrid Store, 2022-1-17==========
USE Cell_Data_2D
!==========LY Add for Multigrid Store, 2022-1-17==========

IMPLICIT NONE

!==========LY Add for Multigrid Store, 2022-1-17==========
INTEGER :: count, local_N, N_Tier_Old
TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh_Old, CellMesh_New
TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh_temp
!==========LY Add for Multigrid Store, 2022-1-17==========

REAL(8), DIMENSION(:,:), POINTER     :: HP, AHP
INTEGER, DIMENSION(:,:), POINTER     :: HT, HE, HT_IFE, AHT, AHE, HT_flag
INTEGER, DIMENSION(:,:), ALLOCATABLE :: element_edge_pointer_old

INTEGER, DIMENSION(:,:), POINTER :: temp_AHE, element_edge_pointer, adaptive_to_uniform_pointer, &
                                    uniform_to_adaptive_pointer
INTEGER :: row_number_of_HT_IFE, number_of_element_old, nn, k, number_of_element_adaptive, &
            n, i, this_ele, DG_node, temp
INTEGER :: neighbor_element_index

row_number_of_HT_IFE = SIZE(HT_IFE, 1)
number_of_element_old = SIZE(HT_IFE, 2)

ALLOCATE(element_edge_pointer_old(4, SIZE(HT, 2)))
element_edge_pointer_old = 0

DO k = 1, SIZE(HE, 2)
   
    this_ele = HE(5, k)
    IF (HE(1, k) == HT(1, this_ele)) THEN
        element_edge_pointer_old(1, this_ele) = k;
    ELSEIF (HE(1, k) == HT(2, this_ele)) THEN
        element_edge_pointer_old(2, this_ele) = k;
    ELSEIF (HE(1, k) == HT(3, this_ele)) THEN
        element_edge_pointer_old(3, this_ele) = k;
    ELSEIF (HE(1, k) == HT(4, this_ele)) THEN
        element_edge_pointer_old(4, this_ele) = k;
    END IF
    
ENDDO




nn = 0
k = 0
number_of_element_adaptive = 0

!----------judge number of element in mesh refinement----------
count = 0
DO k = 1, number_of_element_old
  IF (HT_IFE(row_number_of_HT_IFE, k)<=0) THEN
    nn = nn + 1
  ELSE
    nn = nn + 4
    count = count + 1
  END IF
END DO
!----------judge number of element in mesh refinement----------

!===============LY Add for Multigrid Store, 2022-1-17===============
Local_N = SIZE(CellMesh_Old,1)
ALLOCATE(CellMesh_temp(Local_N))
DO i = 1, Local_N
  CellMesh_temp(i)%Globalindex = 0
  CellMesh_temp(i)%Localindex = 0
  CellMesh_temp(i)%isSplitted = 0
  CellMesh_temp(i)%Boundary(1:4) = 0.0
  CellMesh_temp(i)%Tier = 0
  CellMesh_temp(i)%Parent = 0
  CellMesh_temp(i)%Child(1:4) = 0
  CellMesh_temp(i)%Finalindex = 0
  CellMesh_temp(i)%Nodeindex(1:4) = 0
  CellMesh_temp(i)%FinalParent = 0  !LY modification for Speeding Particle Positioning, 2022-7-25
END DO

!We need to protect the 1st CellMesh information.
DO i = 1, Local_N
  CellMesh_temp(i) = CellMesh_Old(i)
END DO

N_Tier_Old = 4*count
ALLOCATE(CellMesh_New(Local_N+N_Tier_Old))
DO i = 1, SIZE(CellMesh_New,1)
  CellMesh_New(i)%Globalindex = 0
  CellMesh_New(i)%Localindex = 0
  CellMesh_New(i)%isSplitted = 0
  CellMesh_New(i)%Boundary(1:4) = 0.0
  CellMesh_New(i)%Tier = 0
  CellMesh_New(i)%Parent = 0
  CellMesh_New(i)%Child(1:4) = 0
  CellMesh_New(i)%Finalindex = 0
  CellMesh_New(i)%Nodeindex(1:4) = 0
  CellMesh_New(i)%FinalParent = 0 !LY modification for Speeding Particle Positioning, 2022-7-25
END DO

DO i = 1, Local_N
  CellMesh_New(i) = CellMesh_temp(i)
END DO
count = Local_N
!===============LY Add for Multigrid Store, 2022-1-17===============

ALLOCATE(AHT(6, nn))
AHT(1:6,1:nn) = 0

!LY REVISE, 2021-12-7
ALLOCATE(HT_flag(6, SIZE(HT,2)))
HT_flag(1:6, 1:SIZE(HT,2)) = 0
DO k = 1, SIZE(HT,2)
  HT_flag(1,k) = k
END DO
!LY REVISE, 2021-12-7

DG_node = nn * 4
DO k=1, SIZE(HT, 2)
    IF (HT(row_number_of_HT_IFE,k) > 0) THEN
      temp = MINVAL(HT(1:4,k))
      DG_node = MIN(DG_node, temp)
    END IF
END DO


!-----------save old element indices of every new elemenet after refinement----------
ALLOCATE(adaptive_to_uniform_pointer(2, nn))
adaptive_to_uniform_pointer(:,:) = 0

!-----------save new element indices of every old elemenet after refinement----------
ALLOCATE(uniform_to_adaptive_pointer(4, number_of_element_old))
uniform_to_adaptive_pointer(:,:) = 0

!----------generate information matrix AHT after mesh refinement----------
nn = 0
DO k = 1, number_of_element_old
  IF (HT_IFE(row_number_of_HT_IFE, k)<=0) THEN
  
    nn = nn + 1
    AHT(1:4, nn) = HT(1:4, k)
    AHT(5, nn) = 0
    AHT(6, nn) = 0
    adaptive_to_uniform_pointer(1:2, nn) = [k, 0]
    uniform_to_adaptive_pointer(1:4, k) = [nn, 0, 0, 0]
    
    !LY add for particle location on PIC. 2021-12-7
    HT_flag(2,k) = 0
    HT_flag(3,k) = nn
    
    !==========LY Add for Multigrid Store, 2022-1-17==========
    !The old element information on new mesh after refinement.
    CellMesh_New(k)%Finalindex = nn
    CellMesh_New(k)%Nodeindex(1:4) = AHT(1:4, nn)
    !==========LY Add for Multigrid Store, 2022-1-17==========
    
  ELSE
  
    !left-bottom child element
    nn = nn + 1
    AHT(1:4,nn)=[DG_node,DG_node + 1,DG_node + 2,DG_node + 3]
    AHT(5, nn) = 1
    AHT(6, nn) = 1
    DG_node = DG_node + 4
    adaptive_to_uniform_pointer(1:2, nn) = [k, 1]
    !==========LY Add for Multigrid Store, 2022-1-17==========
    !We need set the Finalindex and Nodeindex value to 0 for element which isSplitted is 1,
    !because the Finalindex and Nodeindex value is meaningless after mesh refinement.
    CellMesh_New(k)%Finalindex = 0
    CellMesh_New(k)%Nodeindex(1:4) = 0
    
    !The new left-bottom child element
    count = count + 1
    CellMesh_New(k)%Child(1) = count    !The old element's left-bottom child element
    CellMesh_New(count)%Localindex = count
    CellMesh_New(count)%Tier = 2
    CellMesh_New(count)%Parent = k
    CellMesh_New(count)%Finalindex = nn
    CellMesh_New(count)%Nodeindex(1:4) = [DG_node-4, DG_node-3, DG_node-2, DG_node-1]
    !==========LY Add for Multigrid Store, 2022-1-17==========
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    CellMesh_New(count)%FinalParent = CellMesh_New(k)%FinalParent
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    
    
    
    !left-top child element
    nn = nn + 1
    AHT(1:4,nn)=[DG_node,DG_node + 1,DG_node + 2,DG_node + 3]
    AHT(5, nn) = 1
    AHT(6, nn) = 1
    DG_node = DG_node + 4
    adaptive_to_uniform_pointer(1:2, nn) = [k, 2]
    !==========LY Add for Multigrid Store, 2022-1-17==========
    !The new left-top child element
    count = count + 1
    CellMesh_New(k)%Child(2) = count    !The old element's left-top child element
    CellMesh_New(count)%Localindex = count
    CellMesh_New(count)%Tier = 2
    CellMesh_New(count)%Parent = k
    CellMesh_New(count)%Finalindex = nn
    CellMesh_New(count)%Nodeindex(1:4) = [DG_node-4, DG_node-3, DG_node-2, DG_node-1]
    !==========LY Add for Multigrid Store, 2022-1-17==========
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    CellMesh_New(count)%FinalParent = CellMesh_New(k)%FinalParent
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    
    
    
    !right-bottom child element
    nn = nn + 1
    AHT(1:4,nn)=[DG_node,DG_node + 1,DG_node + 2,DG_node + 3]
    AHT(5, nn) = 1
    AHT(6, nn) = 1
    DG_node = DG_node + 4
    adaptive_to_uniform_pointer(1:2, nn) = [k, 3]
    !==========LY Add for Multigrid Store, 2022-1-17==========
    !The new right-bottom child element
    count = count + 1
    CellMesh_New(k)%Child(3) = count    !The old element's right-bottom child element
    CellMesh_New(count)%Localindex = count
    CellMesh_New(count)%Tier = 2
    CellMesh_New(count)%Parent = k
    CellMesh_New(count)%Finalindex = nn
    CellMesh_New(count)%Nodeindex(1:4) = [DG_node-4, DG_node-3, DG_node-2, DG_node-1]
    !==========LY Add for Multigrid Store, 2022-1-17==========
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    CellMesh_New(count)%FinalParent = CellMesh_New(k)%FinalParent
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    
    
    
    !right-top child element
    nn = nn + 1
    AHT(1:4,nn)=[DG_node,DG_node + 1,DG_node + 2,DG_node + 3]
    AHT(5, nn) = 1
    AHT(6, nn) = 1
    DG_node = DG_node + 4
    adaptive_to_uniform_pointer(1:2, nn) = [k, 4]
    !==========LY Add for Multigrid Store, 2022-1-17==========
    !The new right-top child element
    count = count + 1
    CellMesh_New(k)%Child(4) = count    !The old element's right-top child element
    CellMesh_New(count)%Localindex = count
    CellMesh_New(count)%Tier = 2
    CellMesh_New(count)%Parent = k
    CellMesh_New(count)%Finalindex = nn
    CellMesh_New(count)%Nodeindex(1:4) = [DG_node-4, DG_node-3, DG_node-2, DG_node-1]
    !==========LY Add for Multigrid Store, 2022-1-17==========
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    CellMesh_New(count)%FinalParent = CellMesh_New(k)%FinalParent
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    
    uniform_to_adaptive_pointer(1:4, k) = [nn-3, nn-2, nn-1, nn]
    
    !LY add for particle location on PIC. 2021-12-3
    HT_flag(2, k) = 1
    HT_flag(3, k) = nn-3
    HT_flag(4, k) = nn-2
    HT_flag(5, k) = nn-1
    HT_flag(6, k) = nn
    
  END IF
END DO
!----------generate information matrix AHT after mesh refinement----------

!----------generate information matrix AHP after mesh refinement----------
number_of_element_adaptive = nn
ALLOCATE(AHP(2, MAXVAL(AHT)))
AHP(:,:) = 0

n = 0
i = 0
DO n = 1, number_of_element_adaptive
  IF (adaptive_to_uniform_pointer(2,n)==0) THEN
    DO i = 1, 4
      AHP(1:2, AHT(i,n)) = HP(1:2, HT(i,adaptive_to_uniform_pointer(1,n)))
    END DO
  ELSE IF (adaptive_to_uniform_pointer(2,n)==1) THEN
    AHP(1:2, AHT(1,n)) = HP(1:2, HT(1,adaptive_to_uniform_pointer(1,n)))
    AHP(1:2, AHT(2,n)) = (HP(1:2, HT(1,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(2,adaptive_to_uniform_pointer(1,n))))/2
    AHP(1:2, AHT(3,n)) = (HP(1:2, HT(1,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(2,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(3,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(4,adaptive_to_uniform_pointer(1,n))))/4
    AHP(1:2, AHT(4,n)) = (HP(1:2, HT(1,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(4,adaptive_to_uniform_pointer(1,n))))/2
  ELSE IF (adaptive_to_uniform_pointer(2,n)==2) THEN
    AHP(1:2, AHT(1,n)) = (HP(1:2, HT(1,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(4,adaptive_to_uniform_pointer(1,n))))/2
    AHP(1:2, AHT(2,n)) = (HP(1:2, HT(1,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(2,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(3,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(4,adaptive_to_uniform_pointer(1,n))))/4
    AHP(1:2, AHT(3,n)) = (HP(1:2, HT(3,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(4,adaptive_to_uniform_pointer(1,n))))/2
    AHP(1:2, AHT(4,n)) = HP(1:2, HT(4,adaptive_to_uniform_pointer(1,n)))
  ELSE IF (adaptive_to_uniform_pointer(2,n)==3) THEN
    AHP(1:2, AHT(1,n)) = (HP(1:2, HT(1,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(2,adaptive_to_uniform_pointer(1,n))))/2
    AHP(1:2, AHT(2,n)) = HP(1:2, HT(2,adaptive_to_uniform_pointer(1,n)))
    AHP(1:2, AHT(3,n)) = (HP(1:2, HT(2,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(3,adaptive_to_uniform_pointer(1,n))))/2
    AHP(1:2, AHT(4,n)) = (HP(1:2, HT(1,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(2,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(3,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(4,adaptive_to_uniform_pointer(1,n))))/4
  ELSE IF (adaptive_to_uniform_pointer(2,n)==4) THEN
    AHP(1:2, AHT(1,n)) = (HP(1:2, HT(1,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(2,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(3,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(4,adaptive_to_uniform_pointer(1,n))))/4
    AHP(1:2, AHT(2,n)) = (HP(1:2, HT(2,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(3,adaptive_to_uniform_pointer(1,n))))/2
    AHP(1:2, AHT(3,n)) = HP(1:2, HT(3,adaptive_to_uniform_pointer(1,n)))
    AHP(1:2, AHT(4,n)) = (HP(1:2, HT(3,adaptive_to_uniform_pointer(1,n)))+ &
                            HP(1:2, HT(4,adaptive_to_uniform_pointer(1,n))))/2
  END IF
END DO
!----------generate information matrix AHP after mesh refinement----------

!================LY Add for Multigrid Store, 2022-1-23===============
DO i = 1, SIZE(CellMesh_New, 1)
  IF (CellMesh_New(i)%Finalindex /= 0) THEN
    CellMesh_New(i)%Boundary(1) = AHP(1, AHT(1,CellMesh_New(i)%Finalindex))
    CellMesh_New(i)%Boundary(2) = AHP(1, AHT(2,CellMesh_New(i)%Finalindex))
    CellMesh_New(i)%Boundary(3) = AHP(2, AHT(1,CellMesh_New(i)%Finalindex))
    CellMesh_New(i)%Boundary(4) = AHP(2, AHT(4,CellMesh_New(i)%Finalindex))
  END IF
END DO
!================LY Add for Multigrid Store, 2022-1-23===============

!----------generate information matrix AHE after mesh refinement----------
ALLOCATE(temp_AHE(8, 8*number_of_element_adaptive))
temp_AHE(:,:) = 0
!Note: temp_AHE(1:2,:) denote global indices of node1 and node2 of this edge.
!      temp_AHE(3:4,:) denote local indices of node1 and node2 of this edge.
!      temp_AHE(5,:) denote global indices of this edge locate element.
!      temp_AHE(6,:) denote brother edge of this edge.(0: exterior edge, +: interior edge)
!      temp_AHE(7,:) denote type of this edge.

ALLOCATE(element_edge_pointer(8, number_of_element_adaptive))
element_edge_pointer(:,:) = 0
!Note: In first mesh refinement, one element almost have eight edge.(4*2=8)

nn = 0
!Note: nn have been denote number of edge after mesh refinement.

DO n = 1, number_of_element_adaptive

  k = adaptive_to_uniform_pointer(1, n)
  
  IF (HT_IFE(row_number_of_HT_IFE,k)>0) THEN
  
    nn = nn + 1
    temp_AHE(1:5, nn) = [AHT(1,n), AHT(2,n), 0, -1, n]
    temp_AHE(7, nn) = 1
    nn = nn + 1
    temp_AHE(1:5, nn) = [AHT(2,n), AHT(3,n), -1, 0, n]
    temp_AHE(7, nn) = 1
    nn = nn + 1
    temp_AHE(1:5, nn) = [AHT(3,n), AHT(4,n), 0, -1, n]
    temp_AHE(7, nn) = 1
    nn = nn + 1
    temp_AHE(1:5, nn) = [AHT(4,n), AHT(1,n), -1, 0, n]
    temp_AHE(7, nn) = 1
    element_edge_pointer(1:4, n) = [nn-3, nn-2, nn-1, nn]
    
  !ELSE
  !  !bottom edge
  !  IF (HE(6,4*k-3)==0)THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [AHT(1,n), AHT(2,n), 1, 2, n]
  !    temp_AHE(7, nn) = 0
  !    element_edge_pointer(1, n) = nn
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5,HE(6,4*k-3)))<=0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [AHT(1,n), AHT(2,n), 1, 2, n]
  !    temp_AHE(7, nn) = 0
  !    element_edge_pointer(1, n) = nn
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5,HE(6,4*k-3)))>0) THEN ! The element below the k_th element needs to be refinement
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [-100, -100, -100, -100, n]
  !    temp_AHE(7, nn) = 2
  !    element_edge_pointer(1, n) = nn
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [-100, -100, -100, -100, n]
  !    temp_AHE(7, nn) = 2
  !    element_edge_pointer(2, n) = nn
  !  END IF
  !  
  !  !right edge
  !  IF (HE(6,4*k-2)==0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [AHT(2,n), AHT(3,n), 2, 3, n]
  !    temp_AHE(7, nn) = 0
  !    element_edge_pointer(3, n) = nn
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5,HE(6,4*k-2)))<=0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [AHT(2,n), AHT(3,n), 2, 3, n]
  !    temp_AHE(7, nn) = 0
  !    element_edge_pointer(3, n) = nn
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5,HE(6,4*k-2)))>0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [-100, -100, -100, -100, n]
  !    temp_AHE(7, nn) = 2
  !    element_edge_pointer(3, n) = nn
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [-100, -100, -100, -100, n]
  !    temp_AHE(7, nn) = 2
  !    element_edge_pointer(4, n) = nn
  !  END IF
  !  
  !  !top edge
  !  IF (HE(6,4*k-1)==0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [AHT(3,n), AHT(4,n), 3, 4, n]
  !    temp_AHE(7, nn) = 0
  !    element_edge_pointer(5, n) = nn
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5,HE(6,4*k-1)))<=0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [AHT(3,n), AHT(4,n), 3, 4, n]
  !    temp_AHE(7, nn) = 0
  !    element_edge_pointer(5, n) = nn
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5,HE(6,4*k-1)))>0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [-100, -100, -100, -100, n]
  !    temp_AHE(7, nn) = 2
  !    element_edge_pointer(5, n) = nn
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [-100, -100, -100, -100, n]
  !    temp_AHE(7, nn) = 2
  !    element_edge_pointer(6, n) = nn
  !  END IF
  !  
  !  !left edge
  !  IF (HE(6,4*k)==0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [AHT(4,n), AHT(1,n), 4, 1, n]
  !    temp_AHE(7, nn) = 0
  !    element_edge_pointer(7, n) = nn
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5,HE(6,4*k)))<=0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [AHT(4,n), AHT(1,n), 4, 1, n]
  !    temp_AHE(7, nn) = 0
  !    element_edge_pointer(7, n) = nn
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5,HE(6,4*k)))>0) THEN
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [-100, -100, -100, -100, n]
  !    temp_AHE(7, nn) = 2
  !    element_edge_pointer(7, n) = nn
  !    nn = nn + 1
  !    temp_AHE(1:5, nn) = [-100, -100, -100, -100, n]
  !    temp_AHE(7, nn) = 2
  !    element_edge_pointer(8, n) = nn
  !  END IF
    
  END IF
  
END DO

nn = 0
DO n = 1, number_of_element_adaptive
  
  k = adaptive_to_uniform_pointer(1, n)
  
  IF (HT_IFE(row_number_of_HT_IFE, k)>0) THEN
    
    IF (adaptive_to_uniform_pointer(2, n)==1) THEN
      !bottom edge in left-lower element: generate temp_AHE(6,:)
      IF (element_edge_pointer_old(1, k) == 0) THEN
        nn = nn + 1
        temp_AHE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6 ,element_edge_pointer_old(1, k)))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(2, HE(6 ,element_edge_pointer_old(1, k)))
        temp_AHE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6 ,element_edge_pointer_old(1, k)))
        temp_AHE(6, nn) = neighbor_element_index
      END IF
      
      !right edge in left-lower element: generate temp_AHE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(3, k)
      temp_AHE(6, nn) = neighbor_element_index
      
      !top edge in left-lower element: generate temp_AHE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(2, k)
      temp_AHE(6, nn) = neighbor_element_index
      
      !left edge in left-lower element: generate temp_AHE(6,:)
      IF (element_edge_pointer_old(4, k) == 0) THEN
        nn = nn + 1
        temp_AHE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, element_edge_pointer_old(4, k)))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(3, HE(6, element_edge_pointer_old(4, k)))
        temp_AHE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, element_edge_pointer_old(4, k)))
        temp_AHE(6, nn) = neighbor_element_index
      END IF
      
    ELSEIF (adaptive_to_uniform_pointer(2, n)==2) THEN
      !bottom edge in left-upper element: generate temp_AHE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(1, k)
      temp_AHE(6, nn) = neighbor_element_index
      
      !right edge in left-upper element: generate temp_AHE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(4, k)
      temp_AHE(6, nn) = neighbor_element_index
      
      !top edge in left-upper element: generate temp_AHE(6,:)
      IF (element_edge_pointer_old(3, k) == 0) THEN
        nn = nn + 1
        temp_AHE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6 ,element_edge_pointer_old(3, k)))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6 ,element_edge_pointer_old(3, k)))
        temp_AHE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6 ,element_edge_pointer_old(3, k)))
        temp_AHE(6, nn) = neighbor_element_index
      END IF
      
      !left edge in left-upper element: generate temp_AHE(6,:)
      IF (element_edge_pointer_old(4, k) == 0) THEN
        nn = nn + 1
        temp_AHE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6 ,element_edge_pointer_old(4, k)))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(4, HE(6 ,element_edge_pointer_old(4, k)))
        temp_AHE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6 ,element_edge_pointer_old(4, k)))
        temp_AHE(6, nn) = neighbor_element_index
      END IF
      
    ELSE IF (adaptive_to_uniform_pointer(2, n)==3) THEN
      !bottom edge in right-lower element: generate temp_AHE(6,:)
      IF (element_edge_pointer_old(1, k) == 0) THEN
        nn = nn + 1
        temp_AHE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6 ,element_edge_pointer_old(1, k)))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(4, HE(6 ,element_edge_pointer_old(1, k)))
        temp_AHE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6 ,element_edge_pointer_old(1, k)))
        temp_AHE(6, nn) = neighbor_element_index
      END IF
      
      !right edge in right-lower element: generate temp_AHE(6,:)
      IF (element_edge_pointer_old(2, k)==0) THEN
        nn = nn + 1
        temp_AHE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6 ,element_edge_pointer_old(2, k)))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6 ,element_edge_pointer_old(2, k)))
        temp_AHE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6 ,element_edge_pointer_old(2, k)))
        temp_AHE(6, nn) = neighbor_element_index
      END IF
      
      !top edge in right-lower element: generate temp_AHE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(4, k)
      temp_AHE(6, nn) = neighbor_element_index
      
      !left edge in right-lower element: generate temp_AHE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(1, k)
      temp_AHE(6, nn) = neighbor_element_index
    
    ELSEIF (adaptive_to_uniform_pointer(2, n)==4) THEN
      !bottom edge in right-upper element: generate temp_AHE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(3, k)
      temp_AHE(6, nn) = neighbor_element_index
      
      !right edge in right-upper element: generate temp_AHE(6,:)
      IF (element_edge_pointer_old(2, k) == 0) THEN
        nn = nn + 1
        temp_AHE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6 ,element_edge_pointer_old(2, k)))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(2, HE(6 ,element_edge_pointer_old(2, k)))
        temp_AHE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6 ,element_edge_pointer_old(2, k)))
        temp_AHE(6, nn) = neighbor_element_index
      END IF
      
      !top edge in right-upper element: generate temp_AHE(6,:)
      IF (element_edge_pointer_old(3, k)==0) THEN
        nn = nn + 1
        temp_AHE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6 ,element_edge_pointer_old(3, k)))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(3, HE(6 ,element_edge_pointer_old(3, k)))
        temp_AHE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6 ,element_edge_pointer_old(3, k)))
        temp_AHE(6, nn) = neighbor_element_index
      END IF
      
      !left edge in right-upper element: generate temp_AHE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(2, k)
      temp_AHE(6, nn) = neighbor_element_index
      
    END IF
    
  !ELSE
  !  !bottom edge: renew temp_AHE(1:2,:) and generate temp_AHE(6,:)
  !  IF (HE(6, 4*k-3)==0) THEN
  !    nn = nn + 1
  !    temp_AHE(6, nn) = 0
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5, HE(6,4*k-3)))<=0) THEN
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(1, HE(5, HE(6, 4*k-3)))
  !    temp_AHE(6, nn) = element_edge_pointer(5, neighbor_element_index)
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5, HE(6,4*k-3)))>0) THEN
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(2, HE(5, HE(6, 4*k-3)))
  !    temp_AHE(6, nn) = element_edge_pointer(3, neighbor_element_index)
  !    temp_AHE(1:2, nn) = [AHT(4, neighbor_element_index), AHT(3, neighbor_element_index)]
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(4, HE(5, HE(6, 4*k-3)))
  !    temp_AHE(6, nn) = element_edge_pointer(3, neighbor_element_index)
  !    temp_AHE(1:2, nn) = [AHT(4, neighbor_element_index), AHT(3, neighbor_element_index)]
  !  END IF
  !  
  !  !right edge: renew temp_AHE(1:2,:) and generate temp_AHE(6,:)
  !  IF (HE(6, 4*k-2)==0) THEN
  !    nn = nn + 1
  !    temp_AHE(6, nn) = 0
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5, HE(6,4*k-2)))<=0) THEN
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(1, HE(5, HE(6, 4*k-2)))
  !    temp_AHE(6, nn) = element_edge_pointer(7, neighbor_element_index)
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5, HE(6,4*k-2)))>0) THEN
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(1, HE(5, HE(6, 4*k-2)))
  !    temp_AHE(6, nn) = element_edge_pointer(4, neighbor_element_index)
  !    temp_AHE(1:2, nn) = [AHT(1, neighbor_element_index), AHT(4, neighbor_element_index)]
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(2, HE(5, HE(6, 4*k-2)))
  !    temp_AHE(6, nn) = element_edge_pointer(4, neighbor_element_index)
  !    temp_AHE(1:2, nn) = [AHT(1, neighbor_element_index), AHT(4, neighbor_element_index)]
  !  END IF
  !  
  !  !top edge: renew temp_AHE(1:2,:) and generate temp_AHE(6,:)
  !  IF (HE(6, 4*k-1)==0) THEN
  !    nn = nn + 1
  !    temp_AHE(6, nn) = 0
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5, HE(6,4*k-1)))<=0) THEN
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(1, HE(5, HE(6, 4*k-1)))
  !    temp_AHE(6, nn) = element_edge_pointer(1, neighbor_element_index)
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5, HE(6,4*k-1)))>0) THEN
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(3, HE(5, HE(6, 4*k-1)))
  !    temp_AHE(6, nn) = element_edge_pointer(1, neighbor_element_index)
  !    temp_AHE(1:2, nn) = [AHT(2, neighbor_element_index), AHT(1, neighbor_element_index)]
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(1, HE(5, HE(6, 4*k-1)))
  !    temp_AHE(6, nn) = element_edge_pointer(1, neighbor_element_index)
  !    temp_AHE(1:2, nn) = [AHT(2, neighbor_element_index), AHT(1, neighbor_element_index)]
  !  END IF
  !  
  !  !left edge: renew temp_AHE(1:2,:) and generate temp_AHE(6,:)
  !  IF (HE(6, 4*k)==0) THEN
  !    nn = nn + 1
  !    temp_AHE(6, nn) = 0
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5, HE(6,4*k)))<=0) THEN
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(1, HE(5, HE(6, 4*k)))
  !    temp_AHE(6, nn) = element_edge_pointer(3, neighbor_element_index)
  !  ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(5, HE(6,4*k)))>0) THEN
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(4, HE(5, HE(6, 4*k)))
  !    temp_AHE(6, nn) = element_edge_pointer(2, neighbor_element_index)
  !    temp_AHE(1:2, nn) = [AHT(3, neighbor_element_index), AHT(2, neighbor_element_index)]
  !    nn = nn + 1
  !    neighbor_element_index = uniform_to_adaptive_pointer(3, HE(5, HE(6, 4*k)))
  !    temp_AHE(6, nn) = element_edge_pointer(2, neighbor_element_index)
  !    temp_AHE(1:2, nn) = [AHT(3, neighbor_element_index), AHT(2, neighbor_element_index)]
  !  END IF
    
  END IF
  
END DO

ALLOCATE(AHE(8, nn))
AHE(:,:) = 0

DO n = 1, nn
  !DO i = 1, 8
    AHE(:, n) = temp_AHE(:, n)
  !END DO
END DO

DEALLOCATE(temp_AHE, element_edge_pointer, adaptive_to_uniform_pointer, uniform_to_adaptive_pointer)

END SUBROUTINE