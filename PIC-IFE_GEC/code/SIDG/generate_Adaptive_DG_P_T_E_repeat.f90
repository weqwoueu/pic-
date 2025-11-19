SUBROUTINE generate_Adaptive_DG_P_T_E_repeat(HP, HT, HE, HT_IFE, repeat_index, AHP, AHT, ADGE, HT_flag, &
                                              CellMesh_Old, CellMesh_New, N_Tier_Old)
  
!------------------LY add AHP, AHT, ADGE for Adaptive Mesh 2021-9-25------------------
!Note: HP, HT, HE, HT_IFE--old, first refinement
!      AHP, AHT, ADGE--new, repeat refinement
!      repeat_index: refinement time

!==========LY Add for Multigrid Store, 2022-1-18==========
USE Cell_Data_2D
!==========LY Add for Multigrid Store, 2022-1-18==========
  
IMPLICIT NONE

!==========LY Add for Multigrid Store, 2022-1-17==========
INTEGER :: count, local_N, count_Tier
INTEGER :: N_Tier_Old     !The (N-2)st refinement generate number of CellMesh
INTEGER :: N_Tier_Pres    !The (N-1)st refinement generate number of CellMesh
TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh_Old   !After finish (N-1)st refinement, we will have CellMesh_Old
TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh_New   !In this, when finish Nst refinement, we could get CellMesh_New
TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh_temp
!==========LY Add for Multigrid Store, 2022-1-17==========

REAL(8), DIMENSION(:,:), POINTER :: HP, AHP
INTEGER, DIMENSION(:,:), POINTER :: HT, HE, HT_IFE, AHT, ADGE, HT_flag
INTEGER :: repeat_index

INTEGER, DIMENSION(:,:), POINTER :: temp_ADGE, element_edge_pointer, adaptive_to_uniform_pointer, &
                                    uniform_to_adaptive_pointer
REAL(8), DIMENSION(:,:), POINTER :: end_coordinates

INTEGER :: row_number_of_HT_IFE, number_of_element_old, nn, k, number_of_element_adaptive, n, i, &
            old_edge, max_old_edge, counter, wrong_information, wrong, &
            neighbor_element_index, neighbor_edge, DG_node, true_bro_edge, temp

REAL(8) :: standard_edge_index, print_standard_edge_index, delta

row_number_of_HT_IFE = SIZE(HT_IFE, 1)
number_of_element_old = SIZE(HT_IFE, 2)

nn = 0
k = 0
n = 0
i = 0
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


DG_node = nn * 4
DO k=1, SIZE(HT, 2)
    IF (HT(row_number_of_HT_IFE,k) > 0) THEN
         temp = MINVAL(HT(1:4,k))
         DG_node = MIN(DG_node, temp)
    END IF
END DO


!----------judge number of element in mesh refinement----------

ALLOCATE(AHT(6, nn))
AHT(:,:) = 0

ALLOCATE(HT_flag(6, SIZE(HT,2)))
HT_flag(:,:) = 0
DO k = 1, SIZE(HT, 2)
  HT_flag(1,k) = k
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
        IF (HT(5, k) == 0) THEN
            nn = nn + 1
            AHT(1:4,nn)= HT(1 : 4, k)
            AHT(5,nn) = 0
            AHT(6,nn) = 0
            adaptive_to_uniform_pointer(1:2, nn) = [k, 0]
            uniform_to_adaptive_pointer(1:4, k) = [nn, 0, 0, 0]
            
            HT_flag(2, k) = 0
            HT_flag(3, k) = nn
            
        ELSE
            nn = nn + 1
            AHT(1:4,nn)=[DG_node,DG_node + 1,DG_node + 2,DG_node + 3]
            AHT(5,nn) = HT(5, k)
            AHT(6,nn) = HT(6, k)
            DG_node = DG_node + 4
            adaptive_to_uniform_pointer(1:2, nn) = [k, 0]
            uniform_to_adaptive_pointer(1:4, k) = [nn, 0, 0, 0]
            
            HT_flag(2, k) = 0
            HT_flag(3, k) = nn
        ENDIF
    ELSE
        nn = nn + 1
        AHT(1:4,nn)=[DG_node,DG_node + 1,DG_node + 2,DG_node + 3]
        AHT(5,nn) = HT(5, k) + 1
        AHT(6,nn) = HT(6, k) + 1
        DG_node = DG_node + 4
        adaptive_to_uniform_pointer(1:2, nn) = [k, 1]
        
        nn = nn + 1
        AHT(1:4,nn)=[DG_node,DG_node + 1,DG_node + 2,DG_node + 3]
        AHT(5,nn) = HT(5, k) + 1
        AHT(6,nn) = HT(6, k) + 1
        DG_node = DG_node + 4
        adaptive_to_uniform_pointer(1:2, nn) = [k, 2]
        
        nn = nn + 1
        AHT(1:4,nn)=[DG_node,DG_node + 1,DG_node + 2,DG_node + 3]
        AHT(5,nn) = HT(5, k) + 1
        AHT(6,nn) = HT(6, k) + 1
        DG_node = DG_node + 4
        adaptive_to_uniform_pointer(1:2, nn) = [k, 3]
        
        nn = nn + 1
        AHT(1:4,nn)=[DG_node,DG_node + 1,DG_node + 2,DG_node + 3]
        AHT(5,nn) = HT(5, k) + 1
        AHT(6,nn) = HT(6, k) + 1
        DG_node = DG_node + 4
        adaptive_to_uniform_pointer(1:2, nn) = [k, 4]
        
        uniform_to_adaptive_pointer(1:4, k) = [nn-3, nn-2, nn-1, nn]
        
        HT_flag(2, k) = 1
        HT_flag(3, k) = nn-3
        HT_flag(4, k) = nn-2
        HT_flag(5, k) = nn-1
        HT_flag(6, k) = nn
    END IF
END DO
!----------generate information matrix AHT after mesh refinement----------

!===============LY Add for Multigrid Store, 2022-1-18===============
!We need to protect the 1st~(N-1)st CellMesh information.
local_N = SIZE(CellMesh_Old,1)

!IF (local_N /= (N_Tier_Old+N_Tier_Pres)) THEN
!  WRITE(6,*) 'The local_N /= (N_Tier_Old+N_Tier_Pres), check Program£¡'
!  STOP
!END IF

ALLOCATE(CellMesh_temp(local_N))
DO i = 1, local_N
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

DO i = 1, local_N
  CellMesh_temp(i) = CellMesh_Old(i)
END DO

N_Tier_Pres = 4*count
DEALLOCATE(CellMesh_New)
ALLOCATE(CellMesh_New(local_N+N_Tier_Pres))
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
  CellMesh_New(i)%FinalParent = 0  !LY modification for Speeding Particle Positioning, 2022-7-25
END DO

DO i = 1, local_N
  CellMesh_New(i) = CellMesh_temp(i)
END DO
!===============LY Add for Multigrid Store, 2022-1-18===============

!===============LY Add for Multigrid Store, 2022-1-18===============
!Befor generate the next Tier information, we need set the present Tier Splitted flag.
!DO i = N_Tier_Old+1, N_Tier_Old+N_Tier_Pres
DO i = Local_N-N_Tier_Old+1, Local_N
  IF (HT_IFE(row_number_of_HT_IFE, CellMesh_New(i)%Finalindex) > 0) THEN
    CellMesh_New(i)%isSplitted = 1
  END IF
END DO

!When we finish Splitted flag setting, we can generate refinement new element through refinement.
count = local_N
count_Tier = MAXVAL(CellMesh_New%Tier)
!DO i = N_Tier_Old+1, N_Tier_Old+N_Tier_Pres
DO i = Local_N-N_Tier_Old+1, Local_N
  
  IF (CellMesh_New(i)%isSplitted == 1) THEN
    
    !The new left-bottom child element.
    count = count + 1
    CellMesh_New(i)%Child(1) = count
    CellMesh_New(count)%Localindex = count
    CellMesh_New(count)%Tier = count_Tier+1
    CellMesh_New(count)%Parent = i
    CellMesh_New(count)%Finalindex = uniform_to_adaptive_pointer(1, CellMesh_New(i)%Finalindex)
    CellMesh_New(count)%Nodeindex(1:4) = AHT(1:4, CellMesh_New(count)%Finalindex)
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    CellMesh_New(count)%FinalParent = CellMesh_New(i)%FinalParent
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    
    !The new left-top child element.
    count = count + 1
    CellMesh_New(i)%Child(2) = count
    CellMesh_New(count)%Localindex = count
    CellMesh_New(count)%Tier = count_Tier+1
    CellMesh_New(count)%Parent = i
    CellMesh_New(count)%Finalindex = uniform_to_adaptive_pointer(2, CellMesh_New(i)%Finalindex)
    CellMesh_New(count)%Nodeindex(1:4) = AHT(1:4, CellMesh_New(count)%Finalindex)
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    CellMesh_New(count)%FinalParent = CellMesh_New(i)%FinalParent
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    
    !The new right-bottom child element.
    count = count + 1
    CellMesh_New(i)%Child(3) = count
    CellMesh_New(count)%Localindex = count
    CellMesh_New(count)%Tier = count_Tier+1
    CellMesh_New(count)%Parent = i
    CellMesh_New(count)%Finalindex = uniform_to_adaptive_pointer(3, CellMesh_New(i)%Finalindex)
    CellMesh_New(count)%Nodeindex(1:4) = AHT(1:4, CellMesh_New(count)%Finalindex)
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    CellMesh_New(count)%FinalParent = CellMesh_New(i)%FinalParent
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    
    !The new right-top child element.
    count = count + 1
    CellMesh_New(i)%Child(4) = count
    CellMesh_New(count)%Localindex = count
    CellMesh_New(count)%Tier = count_Tier+1
    CellMesh_New(count)%Parent = i
    CellMesh_New(count)%Finalindex = uniform_to_adaptive_pointer(4, CellMesh_New(i)%Finalindex)
    CellMesh_New(count)%Nodeindex(1:4) = AHT(1:4, CellMesh_New(count)%Finalindex)
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    CellMesh_New(count)%FinalParent = CellMesh_New(i)%FinalParent
    !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
    
    !We need set the Finalindex and Nodeindex value to 0 for element which isSplitted is 1,
    !because the Finalindex and Nodeindex value is meaningless after mesh refinement.
    CellMesh_New(i)%Finalindex = 0
    CellMesh_New(i)%Nodeindex(1:4) = 0
    
  END IF
  
END DO

DO i = 1, Local_N
  IF (CellMesh_New(i)%isSplitted == 0 .AND. CellMesh_New(i)%Finalindex /= 0) THEN
    !There is this element which isSplitted is 0 on old mesh, after refinement, the Finalindex is need update.
    CellMesh_New(i)%Finalindex = uniform_to_adaptive_pointer(1, CellMesh_New(i)%Finalindex)
    CellMesh_New(i)%Nodeindex(1:4) = AHT(1:4, CellMesh_New(i)%Finalindex)
  END IF
END DO

N_Tier_Old = N_Tier_Pres
!===============LY Add for Multigrid Store, 2022-1-18===============

!----------generate information matrix AHP after mesh refinement----------
number_of_element_adaptive = nn
ALLOCATE(AHP(2, MAXVAL(AHT)))
AHP(:,:) = 0.0

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

!----------generate information matrix ADGE after mesh refinement----------
ALLOCATE(temp_ADGE(8, 8*number_of_element_adaptive))
temp_ADGE(:,:) = 0
!Note: temp_ADGE(1:2,:) denote global indices of node1 and node2 of this edge.
!      temp_ADGE(3:4,:) denote local indices of node1 and node2 of this edge.
!      temp_ADGE(5,:) denote global indices of this edge locate element.
!      temp_ADGE(6,:) denote brother edge of this edge.(0: exterior edge, +: interior edge)
!      temp_ADGE(7,:) denote type of this edge.

ALLOCATE(element_edge_pointer((2**repeat_index)*4, number_of_element_adaptive))
element_edge_pointer(:,:) = 0
!Note: In repeating mesh refinement, the number that one element almost have edge will gradually increase.

ALLOCATE(end_coordinates(4, 8*number_of_element_adaptive))
end_coordinates(:,:) = 0.0
!Note: end_coordinates will be store two endpointer coordinates of every edge.

nn = 0
old_edge = 1
max_old_edge = SIZE(HE, 2)

DO n = 1, number_of_element_adaptive

  k = adaptive_to_uniform_pointer(1, n)
  
  IF (HT_IFE(row_number_of_HT_IFE, k)>0) THEN      
      
    nn = nn + 1
    temp_ADGE(1:5, nn) = [AHT(1,n), AHT(2,n), 0, -1, n]
    temp_ADGE(7, nn) = repeat_index
    end_coordinates(1:4, nn) = [AHP(1:2,AHT(1,n)), AHP(1:2,AHT(2,n))]
    nn = nn + 1
    temp_ADGE(1:5, nn) = [AHT(2,n), AHT(3,n), -1, 0, n]
    temp_ADGE(7, nn) = repeat_index
    end_coordinates(1:4, nn) = [AHP(1:2,AHT(2,n)), AHP(1:2,AHT(3,n))]
    nn = nn + 1
    temp_ADGE(1:5, nn) = [AHT(3,n), AHT(4,n), 0, -1, n]
    temp_ADGE(7, nn) = repeat_index
    end_coordinates(1:4, nn) = [AHP(1:2,AHT(3,n)), AHP(1:2,AHT(4,n))]
    nn = nn + 1
    temp_ADGE(1:5, nn) = [AHT(4,n), AHT(1,n), -1, 0, n]
    temp_ADGE(7, nn) = repeat_index
    end_coordinates(1:4, nn) = [AHP(1:2,AHT(4,n)), AHP(1:2,AHT(1,n))]
    element_edge_pointer(1:4, n) = [nn-3, nn-2, nn-1, nn]
    
    old_edge = old_edge + 1
    
  ELSE IF (HT_IFE(row_number_of_HT_IFE, k)<=0 .AND. HT(5, k) > 0) THEN ! non-interface DG element
  
      !IF (nn > 80) then
      !  write(*, *) 123    
      !END IF
      
    counter = 0
    standard_edge_index = 1
    DO WHILE ((old_edge<=max_old_edge).AND.(HE(5,old_edge)==k)) 
      
      !IF (HE(6, old_edge)==0 .OR. HT_IFE(row_number_of_HT_IFE, HE(6,old_edge))<=0) THEN
      IF (HE(6, old_edge)==0) THEN
        counter = counter + 1
        nn = nn + 1
        IF (HE(8, old_edge)==-100) THEN
          temp_ADGE(1:5, nn) = [-100, -100, -100, -100, n]
          temp_ADGE(3:4,nn) = HE(3:4, old_edge)
          temp_ADGE(7,nn) = HE(7, old_edge)
          temp_ADGE(8,nn) = -100 ! importance
        ELSE
          IF (ABS(standard_edge_index-1)<1.0D-12) THEN
            temp_ADGE(1:5, nn) = [AHT(1,n), AHT(2,n), 0, -1, n]
            temp_ADGE(7, nn) = 1
          ELSE IF (ABS(standard_edge_index-2)<1.0D-12) THEN
            temp_ADGE(1:5, nn) = [AHT(2,n), AHT(3,n), -1, 0, n]
            temp_ADGE(7, nn) = 1
          ELSE IF (ABS(standard_edge_index-3)<1.0D-12) THEN
            temp_ADGE(1:5, nn) = [AHT(3,n), AHT(4,n), 0, -1, n]
            temp_ADGE(7, nn) = 1
          ELSE IF (ABS(standard_edge_index-4)<1.0D-12) THEN
            temp_ADGE(1:5, nn) = [AHT(4,n), AHT(1,n), -1, 0, n]
            temp_ADGE(7, nn) = 1
          ELSE
            print_standard_edge_index = standard_edge_index
            wrong_information = 11111
            WRITE(*,*) 'print_standard_edge_index = ',standard_edge_index
            STOP
          END IF
        END IF
        end_coordinates(1:4, nn) = [HP(1:2,HE(1,old_edge)), HP(1:2,HE(2,old_edge))]
        element_edge_pointer(counter, n) = nn
        
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6,old_edge))<=0) THEN
        counter = counter + 1
        nn = nn + 1
        IF (HE(8, old_edge)==-100) THEN
          temp_ADGE(1:5, nn) = [-100, -100, -100, -100, n]
          temp_ADGE(3:4,nn) = HE(3:4, old_edge)
          temp_ADGE(7,nn) = HE(7, old_edge)
          temp_ADGE(8,nn) = -100 ! importance
        ELSE
          IF (ABS(standard_edge_index-1)<1.0D-12) THEN
            temp_ADGE(1:5, nn) = [AHT(1,n), AHT(2,n), 0, -1, n]
            temp_ADGE(7, nn) = 1
          ELSE IF (ABS(standard_edge_index-2)<1.0D-12) THEN
            temp_ADGE(1:5, nn) = [AHT(2,n), AHT(3,n), -1, 0, n]
            temp_ADGE(7, nn) = 1
          ELSE IF (ABS(standard_edge_index-3)<1.0D-12) THEN
            temp_ADGE(1:5, nn) = [AHT(3,n), AHT(4,n), 0, -1, n]
            temp_ADGE(7, nn) = 1
          ELSE IF (ABS(standard_edge_index-4)<1.0D-12) THEN
            temp_ADGE(1:5, nn) = [AHT(4,n), AHT(1,n), -1, 0, n]
            temp_ADGE(7, nn) = 1
          ELSE
            print_standard_edge_index = standard_edge_index
            wrong_information = 11111
            WRITE(*,*) 'print_standard_edge_index = ',standard_edge_index
            STOP
          END IF
        END IF
        end_coordinates(1:4, nn) = [HP(1:2,HE(1,old_edge)), HP(1:2,HE(2,old_edge))]
        element_edge_pointer(counter, n) = nn
        
      !ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge)) == 0 ) THEN
      !
      !  counter = counter + 1
      !  nn = nn + 1
      !  IF (HE(3, old_edge)==-100) THEN
      !    temp_ADGE(1:5, nn) = [-100, -100, -100, -100, n]
      !    temp_ADGE(7, nn) = 0
      !  ELSE
      !    IF (standard_edge_index==1) THEN
      !      temp_ADGE(1:5, nn) = [AHT(1,n), AHT(2,n), 1, 2, n]
      !      temp_ADGE(7, nn) = 0
      !    ELSE IF (standard_edge_index==2) THEN
      !      temp_ADGE(1:5, nn) = [AHT(2,n), AHT(3,n), 2, 3, n]
      !      temp_ADGE(7, nn) = 0
      !    ELSE IF (standard_edge_index==3) THEN
      !      temp_ADGE(1:5, nn) = [AHT(3,n), AHT(4,n), 3, 4, n]
      !      temp_ADGE(7, nn) = 0
      !    ELSE IF (standard_edge_index==4) THEN
      !      temp_ADGE(1:5, nn) = [AHT(4,n), AHT(1,n), 4, 1, n]
      !      temp_ADGE(7, nn) = 0
      !    ELSE
      !      print_standard_edge_index = standard_edge_index
      !      wrong_information = 11111
      !    END IF
      !  END IF
      !  end_coordinates(1:4, nn) = [HP(1:2,HE(1,old_edge)), HP(1:2,HE(2,old_edge))]
      !  element_edge_pointer(counter, n) = nn
        
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge))>0) THEN
        !one edge have been divide two samller edge
        counter = counter + 1
        nn = nn + 1
        temp_ADGE(1:5, nn) = [-100, -100, -100, -100, n]
        temp_ADGE(3:4,nn) = HE(3:4, old_edge)
        temp_ADGE(7, nn) = repeat_index
        temp_ADGE(8, nn) = -100
        end_coordinates(1:4, nn) = [HP(1:2,HE(1,old_edge)), (HP(1:2,HE(1,old_edge))+HP(1:2,HE(2,old_edge)))/2.0]
        element_edge_pointer(counter, n) = nn
        counter = counter + 1
        nn = nn + 1
        temp_ADGE(1:5, nn) = [-100, -100, -100, -100, n]
        temp_ADGE(3:4,nn) = HE(3:4, old_edge)
        temp_ADGE(7, nn) = repeat_index
        temp_ADGE(8, nn) = -100
        end_coordinates(1:4, nn) = [(HP(1:2,HE(1,old_edge))+HP(1:2,HE(2,old_edge)))/2.0, HP(1:2,HE(2,old_edge))]
        element_edge_pointer(counter, n) = nn
      
      END IF
      
      IF (ABS(HP(1,HE(1,old_edge))-HP(1,HE(2,old_edge)))<1.0D-12) THEN
        delta = ABS(HP(2,HE(1,old_edge))-HP(2,HE(2,old_edge)))/ABS(HP(2,HT(3,k))-HP(2,HT(2,k)))
      ELSE IF (ABS(HP(2,HE(1,old_edge))-HP(2,HE(2,old_edge)))<1.0D-12) THEN
        delta = ABS(HP(1,HE(1,old_edge))-HP(1,HE(2,old_edge)))/ABS(HP(1,HT(2,k))-HP(1,HT(1,k)))
      ELSE
        wrong = 1111111
        Write(*,*) 'Wrong =1111111--delta'
        Stop
      END IF
      !delta = 0.5 or 1.0
      
      standard_edge_index = standard_edge_index + delta
      old_edge = old_edge + 1
      
      IF (old_edge>max_old_edge) THEN
        EXIT
        !Prevents arrays from crossing bounds
      END IF
      
    END DO
    
  END IF
  
END DO

nn = 0
old_edge = 1
neighbor_edge = 0
DO n = 1, number_of_element_adaptive

  k = adaptive_to_uniform_pointer(1, n)
  
  IF (HT_IFE(row_number_of_HT_IFE, k)>0) THEN
  
    IF (adaptive_to_uniform_pointer(2, n)==1) THEN
      !bottom edge in left-lower element: generate temp_ADGE(6,:)
      IF (HE(6,old_edge)==0) THEN
        nn = nn + 1
        temp_ADGE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(2, HE(6, old_edge))
        temp_ADGE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge))
        temp_ADGE(6, nn) = neighbor_element_index
        !i = 1
        !neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !DO WHILE (neighbor_edge/=0)
        !  IF (end_coordinates(1,nn)==end_coordinates(3,neighbor_edge).AND.end_coordinates(2,nn)== &
        !        end_coordinates(4,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(1,neighbor_edge).AND. &
        !        end_coordinates(4,nn)==end_coordinates(2,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !        ELSEIF (end_coordinates(1,nn)==end_coordinates(1,neighbor_edge).AND.end_coordinates(2,nn)== &
        !                end_coordinates(2,neighbor_edge).AND.end_coordinates(3,nn)== &
        !                end_coordinates(3,neighbor_edge).AND.end_coordinates(4,nn)== &
        !                end_coordinates(4,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  END IF
        !  i = i + 1
        !  neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !END DO
      END IF
      
      !right edge in left-lower element: generate temp_ADGE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(3, k)
      temp_ADGE(6, nn) = neighbor_element_index
      
      !top edge in left-lower element: generate temp_ADGE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(2, k)
      temp_ADGE(6, nn) = neighbor_element_index
      
      !left edge in left-lower element: generate temp_ADGE(6,:)
      IF (HE(6,old_edge+3)==0) THEN
        nn = nn + 1
        temp_ADGE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge+3))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(3, HE(6, old_edge+3))
        temp_ADGE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge+3))
        temp_ADGE(6, nn) = neighbor_element_index
        !i = 1
        !neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !DO WHILE (neighbor_edge/=0)
        !  IF (end_coordinates(1,nn)==end_coordinates(3,neighbor_edge).AND.end_coordinates(2,nn)== &
        !        end_coordinates(4,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(1,neighbor_edge).AND. &
        !        end_coordinates(4,nn)==end_coordinates(2,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  ELSEIF (end_coordinates(1,nn)==end_coordinates(1,neighbor_edge).AND.end_coordinates(2,nn)== &
        !            end_coordinates(2,neighbor_edge).AND.end_coordinates(3,nn)== &
        !            end_coordinates(3,neighbor_edge).AND.end_coordinates(4,nn)== &
        !            end_coordinates(4,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  END IF
        !  i = i + 1
        !  neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !END DO
      END IF
      
    ELSE IF (adaptive_to_uniform_pointer(2,n)==2) THEN
      !bottom edge in left-upper element: generate temp_ADGE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(1, k)
      temp_ADGE(6, nn) = neighbor_element_index
      
      !right edge in left-upper element: generate temp_ADGE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(4, k)
      temp_ADGE(6, nn) = neighbor_element_index
      
      !top edge in left-upper element: generate temp_ADGE(6,:)
      IF (HE(6,old_edge+1)==0) THEN
        nn = nn + 1
        temp_ADGE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge+1))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge+1))
        temp_ADGE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge+1))
        temp_ADGE(6, nn) = neighbor_element_index
        !i = 1
        !neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !DO WHILE (neighbor_edge/=0)
        !  IF (end_coordinates(1,nn)==end_coordinates(3,neighbor_edge).AND.end_coordinates(2,nn)== &
        !        end_coordinates(4,neighbor_edge).AND.end_coordinates(3,nn)== &
        !        end_coordinates(1,neighbor_edge).AND.end_coordinates(4,nn)==end_coordinates(2,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  ELSEIF (end_coordinates(1,nn)==end_coordinates(1,neighbor_edge).AND.end_coordinates(2,nn)== &
        !            end_coordinates(2,neighbor_edge).AND.end_coordinates(3,nn)== &
        !            end_coordinates(3,neighbor_edge).AND.end_coordinates(4,nn)==end_coordinates(4,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  END IF
        !  i = i + 1
        !  neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !END DO
      END IF
      
      !left edge in left-upper element: generate temp_ADGE(6,:)
      IF (HE(6,old_edge+2)==0) THEN
        nn = nn + 1
        temp_ADGE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge+2))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(4, HE(6, old_edge+2))
        temp_ADGE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge+2))
        temp_ADGE(6, nn) = neighbor_element_index
        !i = 1
        !neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !DO WHILE (neighbor_edge/=0)
        !  IF (end_coordinates(1,nn)==end_coordinates(3,neighbor_edge).AND.end_coordinates(2,nn)== &
        !        end_coordinates(4,neighbor_edge).AND.end_coordinates(3,nn)== &
        !        end_coordinates(1,neighbor_edge).AND.end_coordinates(4,nn)==end_coordinates(2,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  ELSEIF (end_coordinates(1,nn)==end_coordinates(1,neighbor_edge).AND.end_coordinates(2,nn)== &
        !            end_coordinates(2,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(3,neighbor_edge).AND. &
        !            end_coordinates(4,nn)==end_coordinates(4,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  END IF
        !  i = i + 1
        !  neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !END DO
      END IF
      
    ELSE IF (adaptive_to_uniform_pointer(2,n)==3) THEN
      !bottom edge in right-lower element: generate temp_ADGE(6,:)
      IF (HE(6,old_edge-2)==0) THEN
        nn = nn + 1
        temp_ADGE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge-2))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(4, HE(6, old_edge-2))
        temp_ADGE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge-2))
        temp_ADGE(6, nn) = neighbor_element_index
        !i = 1
        !neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !DO WHILE (neighbor_edge/=0)
        !  IF (end_coordinates(1,nn)==end_coordinates(3,neighbor_edge).AND.end_coordinates(2,nn)== &
        !        end_coordinates(4,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(1,neighbor_edge).AND. &
        !        end_coordinates(4,nn)==end_coordinates(2,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  ELSEIF (end_coordinates(1,nn)==end_coordinates(1,neighbor_edge).AND.end_coordinates(2,nn)== &
        !            end_coordinates(2,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(3,neighbor_edge) &
        !            .AND.end_coordinates(4,nn)==end_coordinates(4,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  END IF
        !  i = i + 1
        !  neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !END DO
      END IF
      
      !right edge in right-lower element: generate temp_ADGE(6,:)
      IF (HE(6,old_edge-1)==0) THEN
        nn = nn + 1
        temp_ADGE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge-1))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge-1))
        temp_ADGE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge-1))
        temp_ADGE(6, nn) = neighbor_element_index
        !i = 1
        !neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !DO WHILE (neighbor_edge/=0)
        !  IF (end_coordinates(1,nn)==end_coordinates(3,neighbor_edge).AND.end_coordinates(2,nn)== &
        !        end_coordinates(4,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(1,neighbor_edge).AND. &
        !        end_coordinates(4,nn)==end_coordinates(2,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  ELSEIF (end_coordinates(1,nn)==end_coordinates(1,neighbor_edge).AND.end_coordinates(2,nn)==&
        !                end_coordinates(2,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(3,neighbor_edge) &
        !                .AND.end_coordinates(4,nn)==end_coordinates(4,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  END IF
        !  i = i + 1
        !  neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !END DO
      END IF
      
      !top edge in right-lower element: generate temp_ADGE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(4, k)
      temp_ADGE(6, nn) = neighbor_element_index
      
      !left edge in right-lower element: generate temp_ADGE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(1, k)
      temp_ADGE(6, nn) = neighbor_element_index
      
    ELSE IF (adaptive_to_uniform_pointer(2,n)==4) THEN
      !bottom edge in right-upper element: generate temp_ADGE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(3, k)
      temp_ADGE(6, nn) = neighbor_element_index
      
      !right edge in right-upper element: generate temp_ADGE(6,:)
      IF (HE(6,old_edge-2)==0) THEN
        nn = nn + 1
        temp_ADGE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge-2))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(2, HE(6, old_edge-2))
        temp_ADGE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge-2))
        temp_ADGE(6, nn) = neighbor_element_index
        !i = 1
        !neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !DO WHILE (neighbor_edge/=0)
        !  IF (end_coordinates(1,nn)==end_coordinates(3,neighbor_edge).AND.end_coordinates(2,nn)== &
        !        end_coordinates(4,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(1,neighbor_edge) &
        !        .AND.end_coordinates(4,nn)==end_coordinates(2,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  ELSEIF (end_coordinates(1,nn)==end_coordinates(1,neighbor_edge).AND.end_coordinates(2,nn)== &
        !            end_coordinates(2,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(3,neighbor_edge) &
        !            .AND.end_coordinates(4,nn)==end_coordinates(4,neighbor_edge)) THEN
        !    temp_ADGE(6, nn) = neighbor_edge
        !  END IF
        !  i = i + 1
        !  neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !END DO
      END IF
      
      !top edge in right-upper element: generate temp_ADGE(6,:)
      IF (HE(6,old_edge-1)==0) THEN
        nn = nn + 1
        temp_ADGE(6, nn) = 0
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge-1))>0) THEN
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(3, HE(6, old_edge-1))
        temp_ADGE(6, nn) = neighbor_element_index
      ELSE
        nn = nn + 1
        neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge-1))
        temp_ADGE(6, nn) = neighbor_element_index
        !i = 1
        !neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !DO WHILE (neighbor_edge/=0)
        !  IF (end_coordinates(1,nn)==end_coordinates(3,neighbor_edge).AND.end_coordinates(2,nn)== &
        !  end_coordinates(4,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(1,neighbor_edge).AND. &
        !  end_coordinates(4,nn)==end_coordinates(2,neighbor_edge)) THEN
        !      
        !    temp_ADGE(6, nn) = neighbor_edge
        !  ELSEIF (end_coordinates(1,nn)==end_coordinates(1,neighbor_edge).AND.end_coordinates(2,nn)== &
        !      end_coordinates(2,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(3,neighbor_edge) &
        !      .AND.end_coordinates(4,nn)==end_coordinates(4,neighbor_edge)) THEN
        !      
        !    temp_ADGE(6, nn) = neighbor_edge
        !  END IF
        !  i = i + 1
        !  neighbor_edge = element_edge_pointer(i, neighbor_element_index)
        !END DO
      END IF
      
      !left edge in right-upper element: generate temp_ADGE(6,:)
      nn = nn + 1
      neighbor_element_index = uniform_to_adaptive_pointer(2, k)
      temp_ADGE(6, nn) = neighbor_element_index
      
    END IF
    
    old_edge = old_edge + 1
    
  ELSE IF (HT_IFE(row_number_of_HT_IFE,k) <= 0 .AND. HT(5 ,k) > 0) THEN
    
    standard_edge_index = 1
    
    DO WHILE ((old_edge<=max_old_edge).AND.(HE(5, old_edge)==k))
    
      IF (HE(6, old_edge)==0) THEN
      
        nn = nn + 1
        temp_ADGE(6, nn) = 0
        
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge))<=0) THEN ! bro_ele is non-interface DG element
      
        counter = counter + 1
        nn = nn + 1
        IF (HE(8, old_edge)==-100) THEN
          neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge))
          i = 1
          neighbor_edge = element_edge_pointer(i, neighbor_element_index)
          true_bro_edge = 0
          DO WHILE (neighbor_edge/=0)
            IF (ABS(end_coordinates(1,nn)-end_coordinates(3,neighbor_edge))<1.0D-12.AND.ABS(end_coordinates(2,nn)- &
                end_coordinates(4,neighbor_edge))<1.0D-12.AND.&
                ABS(end_coordinates(3,nn)-end_coordinates(1,neighbor_edge))<1.0D-12.AND. &
                ABS(end_coordinates(4,nn)-end_coordinates(2,neighbor_edge))<1.0D-12) THEN
                
            true_bro_edge = neighbor_edge
            ELSEIF (ABS(end_coordinates(1,nn)-end_coordinates(1,neighbor_edge))<1.0D-12.AND.ABS(end_coordinates(2,nn)- &
                end_coordinates(2,neighbor_edge))<1.0D-12.AND.&
                ABS(end_coordinates(3,nn)-end_coordinates(3,neighbor_edge))<1.0D-12 &
                .AND.ABS(end_coordinates(4,nn)-end_coordinates(4,neighbor_edge))<1.0D-12) THEN
                    
            true_bro_edge = neighbor_edge
          END IF
            i = i + 1
            neighbor_edge = element_edge_pointer(i, neighbor_element_index)
            
          END DO
          
          temp_ADGE(1:2, nn) = [temp_ADGE(2, true_bro_edge), temp_ADGE(1, true_bro_edge)]
          temp_ADGE(6, nn) = neighbor_element_index
        ELSE
          neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge))
          temp_ADGE(6, nn) = neighbor_element_index
          !i = 1
          !neighbor_edge = element_edge_pointer(i, neighbor_element_index)
          !DO WHILE (neighbor_edge/=0)
          !  IF (end_coordinates(1,nn)==end_coordinates(3,neighbor_edge).AND.end_coordinates(2,nn)== &
          !      end_coordinates(4,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(1,neighbor_edge) &
          !      .AND.end_coordinates(4,nn)==end_coordinates(2,neighbor_edge)) THEN
          !  temp_ADGE(6, nn) = neighbor_edge
          !  ELSEIF (end_coordinates(1,nn)==end_coordinates(1,neighbor_edge).AND.end_coordinates(2,nn)== &
          !          end_coordinates(2,neighbor_edge).AND.end_coordinates(3,nn)==end_coordinates(3,neighbor_edge) &
          !          .AND.end_coordinates(4,nn)==end_coordinates(4,neighbor_edge)) THEN
          !  temp_ADGE(6, nn) = neighbor_edge
          !END IF
          !  i = i + 1
          !  neighbor_edge = element_edge_pointer(i, neighbor_element_index)
          !END DO
        END IF
        
      ELSE IF (HT_IFE(row_number_of_HT_IFE, HE(6, old_edge))>0) THEN
      
        nn = nn + 1
        IF ((standard_edge_index>=0.9).AND.(standard_edge_index<1.9)) THEN
            
          neighbor_element_index = uniform_to_adaptive_pointer(2, HE(6, old_edge))
          neighbor_edge = element_edge_pointer(3,neighbor_element_index)
          temp_ADGE(6, nn) = neighbor_element_index
          
        ELSE IF ((standard_edge_index>=1.9).AND.(standard_edge_index<2.9)) THEN
            
          neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge))
          neighbor_edge= element_edge_pointer(4,neighbor_element_index)
          temp_ADGE(6, nn) = neighbor_element_index
          
        ELSE IF ((standard_edge_index>=2.9).AND.(standard_edge_index<3.9)) THEN
            
          neighbor_element_index = uniform_to_adaptive_pointer(3, HE(6, old_edge))
          neighbor_edge= element_edge_pointer(1,neighbor_element_index)
          temp_ADGE(6, nn) = neighbor_element_index
          
        ELSE IF ((standard_edge_index>=3.9).AND.(standard_edge_index<4.9)) THEN
            
          neighbor_element_index = uniform_to_adaptive_pointer(4, HE(6, old_edge))
          neighbor_edge = element_edge_pointer(2,neighbor_element_index)
          temp_ADGE(6, nn) = neighbor_element_index
          
        ELSE
          wrong_information = 22222
          WRITE(*,*) 'wrong_information = 22222'
          STOP
        END IF
        temp_ADGE(1:2, nn) = [temp_ADGE(2, neighbor_edge), temp_ADGE(1, neighbor_edge)]
        
        nn = nn + 1
        IF ((standard_edge_index>=0.9).AND.(standard_edge_index<1.9)) THEN
          neighbor_element_index = uniform_to_adaptive_pointer(4, HE(6, old_edge))
          neighbor_edge = element_edge_pointer(3,neighbor_element_index)
          temp_ADGE(6, nn) = neighbor_element_index
        ELSE IF ((standard_edge_index>=1.9).AND.(standard_edge_index<2.9)) THEN
          neighbor_element_index = uniform_to_adaptive_pointer(2, HE(6, old_edge))
          neighbor_edge= element_edge_pointer(4,neighbor_element_index)
          temp_ADGE(6, nn) = neighbor_element_index
        ELSE IF ((standard_edge_index>=2.9).AND.(standard_edge_index<3.9)) THEN
          neighbor_element_index = uniform_to_adaptive_pointer(1, HE(6, old_edge))
          neighbor_edge = element_edge_pointer(1,neighbor_element_index)
          temp_ADGE(6, nn) = neighbor_element_index
        ELSE IF ((standard_edge_index>=3.9).AND.(standard_edge_index<4.9)) THEN
          neighbor_element_index = uniform_to_adaptive_pointer(3, HE(6, old_edge))
          neighbor_edge = element_edge_pointer(2,neighbor_element_index)
          temp_ADGE(6, nn) = neighbor_element_index
        ELSE
          wrong_information = 33333
          WRITE(*,*) 'wrong_information = 33333'
          STOP
        END IF
        temp_ADGE(1:2, nn) = [temp_ADGE(2, neighbor_edge), temp_ADGE(1, neighbor_edge)]
        
      END IF
      
      IF (ABS(HP(1,HE(1,old_edge))-HP(1,HE(2,old_edge)))<1.0D-12) THEN
        delta = ABS(HP(2,HE(1,old_edge))-HP(2,HE(2,old_edge)))/ABS(HP(2,HT(3,k))-HP(2,HT(2,k)))
      ELSE IF (ABS(HP(2,HE(1,old_edge))-HP(2,HE(2,old_edge)))<1.0D-12) THEN
        delta = ABS(HP(1,HE(1,old_edge))-HP(1,HE(2,old_edge)))/ABS(HP(1,HT(2,k))-HP(1,HT(1,k)))
      ELSE
        wrong = 1111111
        WRITE(*,*) 'wrong = 1111111---delta'
        STOP
      END IF
      !delta = 0.5 or 1.0
      standard_edge_index = standard_edge_index + delta
      old_edge = old_edge + 1
      
      IF (old_edge>max_old_edge) THEN
        EXIT
        !Prevents arrays from crossing bounds
      END IF
      
    END DO
    
  END IF
  
    IF (old_edge > max_old_edge) THEN  ! importance revise
        EXIT
    END IF
  
  
END DO

ALLOCATE(ADGE(8, nn))
ADGE(:,:) = 0

DO n = 1, nn
  DO i = 1, 8
    ADGE(i, n) = temp_ADGE(i, n)
  END DO
END DO

DEALLOCATE(temp_ADGE,element_edge_pointer, end_coordinates, adaptive_to_uniform_pointer, uniform_to_adaptive_pointer)

END SUBROUTINE