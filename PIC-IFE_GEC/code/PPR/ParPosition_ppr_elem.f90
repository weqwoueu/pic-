SUBROUTINE ParPosition_ppr_elem(node_corx, node_cory, elem_index_1, elem_index_2, elem_index_3, elem_index_4)
 
  
USE IFE_Data
USE Particle_2D
USE Field_2D
USE IFE_MAIN_PARAM
USE Domain_2D
USE Constant_Variable_2D
IMPLICIT NONE

INTEGER :: i_part
INTEGER :: i_1, j_1, i, j, num, k
INTEGER :: n_element_old
INTEGER :: nnx, nny, nnz
INTEGER :: N_part_tot_ppr
INTEGER :: elem_index_1, elem_index_2, elem_index_3, elem_index_4
INTEGER, DIMENSION(:), ALLOCATABLE    :: part_ele_location_ppr
REAL(8) :: node_corx, node_cory
REAL(8) :: bias
REAL(8) :: zmin, zmax
REAL(8) :: rxp_1, ryp_1, part_ele_left, part_ele_right, part_ele_bottom, part_ele_top
REAL(8) :: bottom_top_centre, left_right_centre
REAL(8) :: delta_1, delta_2
REAL(8) :: x, y, xmin, xmax, ymin, ymax
REAL(8), DIMENSION(:,:), ALLOCATABLE	:: part_ppr


OPEN(1, ACTION = 'READ', FILE = './INPUT/mesh.inp')
READ(1,*) (Vert_o(i), i = 1,3)
READ(1,*)  xmax, ymax, zmax
! nx, ny, nz							: Number of mesh points in x, y, z
READ(1,*) nnx, nny, nnz
! hx(1:3)								: Grid resolution in each direction
READ(1,*) (hx(i), i = 1, 2)
READ(1,*) hz
CLOSE(1)
!==============================================================zyz============
xmin = Vert_o(1)
ymin = Vert_o(2)
!==============================================================zyz============
nx = nnx
ny = nny
nz = nnz

hxi = One/hx
hzi = One/hz
! Origin for the PIC domain
Vert_o(1:2) = Vert_o(1:2) - hx
Vert_o(3) = Vert_o(3) - hz
!===============zyz vert_o===============================================
N_part_tot_ppr = 4
ALLOCATE(part_ele_location_ppr(N_part_tot_ppr))
ALLOCATE(part_ppr(4,2))
part_ppr= 0.D0
part_ele_location_ppr=0.D0

OPEN(1,ACTION='READ',FILE='./INPUT/IDG_inf.inp')
    READ(1, *)
    READ(1, *) repeat_refinement
    READ(1, *) 
CLOSE(1)
bias=hx(1)/(2**(repeat_refinement+2))
!bias changes when the repeat refinment times change,so the particles spread to four directions do not go to the lines,they only go inside meshes  

IF (ABS(node_corx-xmin)<SmallValue .AND. ABS(node_cory-ymin)<SmallValue) THEN!bottom left node
  
  part_ppr(1,1) = xmin + bias
  part_ppr(1,2) = ymin + bias
  
ELSEIF (ABS(node_corx-xmin)<SmallValue .AND. ABS(node_cory-ymax)<SmallValue) THEN!top left node
  
  part_ppr(1,1) = xmin + bias
  part_ppr(1,2) = ymax - bias
  
ELSEIF (ABS(node_corx-xmax)<SmallValue .AND. ABS(node_cory-ymax)<SmallValue) THEN!top right node
  
  part_ppr(1,1) = xmax - bias
  part_ppr(1,2) = ymax - bias
  
ELSEIF (ABS(node_corx-xmax)<SmallValue .AND. ABS(node_cory-ymin)<SmallValue) THEN!bottom right node
 
  part_ppr(1,1) = xmax - bias
  part_ppr(1,2) = ymin + bias
  
ELSEIF (ABS(node_corx-xmin)<SmallValue) THEN!left 
  
  part_ppr(2,1) = node_corx + bias
  part_ppr(2,2) = node_cory - bias
  part_ppr(3,1) = node_corx + bias
  part_ppr(3,2) = node_cory + bias
  
ELSEIF (ABS(node_cory-ymax)<SmallValue)THEN!top
  
  part_ppr(1,1) = node_corx - bias
  part_ppr(1,2) = node_cory - bias
  part_ppr(2,1) = node_corx + bias
  part_ppr(2,2) = node_cory - bias

ELSEIF (ABS(node_corx-xmax)<SmallValue) THEN!right
  
  part_ppr(1,1) = node_corx - bias
  part_ppr(1,2) = node_cory - bias
  part_ppr(4,1) = node_corx - bias
  part_ppr(4,2) = node_cory + bias
  

ELSEIF (ABS(node_cory-ymin)<SmallValue)THEN!bottom
  
  part_ppr(3,1) = node_corx + bias
  part_ppr(3,2) = node_cory + bias
  part_ppr(4,1) = node_corx - bias
  part_ppr(4,2) = node_cory + bias
  
ELSE!inside
  
  part_ppr(1,1) = node_corx - bias
  part_ppr(1,2) = node_cory - bias

  part_ppr(2,1) = node_corx + bias
  part_ppr(2,2) = node_cory - bias

  part_ppr(3,1) = node_corx + bias
  part_ppr(3,2) = node_cory + bias

  part_ppr(4,1) = node_corx - bias
  part_ppr(4,2) = node_cory + bias
ENDIF

!========== Particle location==========
!zyz ========================================================================
  
DO i_part = 1, N_part_tot_ppr

  !Firstly, we need to locate particle on the initial mesh: 1st mesh(no refinement).
  rxp_1 = (part_ppr(i_part,1) - Vert_o(1)) * hxi(1)
  i_1 = INT(rxp_1)
  
  ryp_1 = (part_ppr(i_part,2) - Vert_o(2)) * hxi(2)
  j_1 = INT(ryp_1)
  
  IF (j_1 == ny) THEN
    j_1 = j_1 - 1
  END IF
  
  IF (i_1 == nx) THEN
    i_1 = i_1 - 1
  END IF
  
  n_element_old = (j_1 + (i_1 - 1)*(ny-1))    !Initial element on Tier1: no refinement
  
  If (i_1==0 .AND. j_1==0) Then
    n_element_old = 1
  End If
  
  DO WHILE (CellMesh(n_element_old)%isSplitted == 1)
  !There denote the initial element is DG element, it will be refine.
    
    !The following four variables denote the distance between particle and DG element boundary. 
    part_ele_left = part_ppr(i_part,1) - CellMesh(n_element_old)%Boundary(1)
    part_ele_right = part_ppr(i_part,1) - CellMesh(n_element_old)%Boundary(2)
    part_ele_bottom = part_ppr(i_part,2) - CellMesh(n_element_old)%Boundary(3)
    part_ele_top = part_ppr(i_part,2) - CellMesh(n_element_old)%Boundary(4)
    
    !The following two variables 'delta_1' and 'delta_2' denote the distance between particle and DG element centre. 
    left_right_centre = (CellMesh(n_element_old)%Boundary(1)+CellMesh(n_element_old)%Boundary(2)) / 2.0
    bottom_top_centre = (CellMesh(n_element_old)%Boundary(3)+CellMesh(n_element_old)%Boundary(4)) / 2.0
    delta_1 = part_ppr(i_part,2) - bottom_top_centre
    delta_2 = part_ppr(i_part,1) - left_right_centre

      !There denote the particle locate inside the DG element.
      !We need to discuss the particle locate inside which one child element.     
      IF (delta_1 > SmallValue) THEN   !Particle locate inside left-top child element or right-top child element.
        IF (delta_2 > SmallValue) THEN   !Particle locate inside right-top child element.
          n_element_old = CellMesh(n_element_old)%Child(4)
        ELSEIF (delta_2 < -SmallValue) THEN   !Particle locate inside left-top child element.
          n_element_old = CellMesh(n_element_old)%Child(2)
        END IF
      ELSEIF (delta_1 < -SmallValue) THEN   !Particle locate inside left-bottom child element or right-bottom child element.
        IF (delta_2 > SmallValue) THEN   !Particle locate inside right-bottom child element.
          n_element_old = CellMesh(n_element_old)%Child(3)
        ELSEIF (delta_2 < -SmallValue) THEN   !Particle locate inside left-bottom child element.
          n_element_old = CellMesh(n_element_old)%Child(1)
        END IF
      END IF 
    !END IF
  END DO
  
  IF (CellMesh(n_element_old)%isSplitted==0 .AND. CellMesh(n_element_old)%Finalindex/=0) THEN
  !There denote the initial element is CG element or the updated element is final child DG element.
    n_element_old = CellMesh(n_element_old)%Finalindex
    part_ele_location_ppr(i_part) = n_element_old
  END IF
END DO
!==========Part1 : Particle location==========

IF (ABS(node_corx-xmin)<SmallValue .AND. ABS(node_cory-ymin)<SmallValue) THEN!bottom left
    
  elem_index_1 = 1
  elem_index_2 = 0
  elem_index_3 = 0 
  elem_index_4 = 0
  
ELSEIF (ABS(node_corx-xmin)<SmallValue .AND. ABS(node_cory-ymax)<SmallValue) THEN!top left
    
  elem_index_1 = part_ele_location_ppr(1)
  elem_index_2 = 0 
  elem_index_3 = 0 
  elem_index_4 = 0
 
ELSEIF (ABS(node_corx-xmax)<SmallValue .AND. ABS(node_cory-ymax)<SmallValue) THEN!top right
  
  elem_index_1 = part_ele_location_ppr(1)
  elem_index_2 = 0 
  elem_index_3 = 0 
  elem_index_4 = 0
 
ELSEIF (ABS(node_corx-xmax)<SmallValue .AND. ABS(node_cory-ymin)<SmallValue) THEN!bottom right
  
  elem_index_1 = part_ele_location_ppr(1)
  elem_index_2 = 0
  elem_index_3 = 0
  elem_index_4 = 0
  
ELSEIF (ABS(node_corx-xmin)<SmallValue) THEN!left 
  
  elem_index_1 = 0
  elem_index_4 = 0
  elem_index_2=part_ele_location_ppr(2)
  elem_index_3=part_ele_location_ppr(3)
  
ELSEIF (ABS(node_cory-ymax)<SmallValue)THEN!top
 
  elem_index_3 = 0
  elem_index_4 = 0
  elem_index_1=part_ele_location_ppr(1)
  elem_index_2=part_ele_location_ppr(2)
  
ELSEIF (ABS(node_corx-xmax)<SmallValue) THEN!right
 
  elem_index_2 = 0
  elem_index_3 = 0
  elem_index_1=part_ele_location_ppr(1)
  elem_index_4=part_ele_location_ppr(4)
  
ELSEIF (ABS(node_cory-ymin)<SmallValue)THEN!bottom
  
  elem_index_1 = 0
  elem_index_2 = 0
  elem_index_3=part_ele_location_ppr(3)
  elem_index_4=part_ele_location_ppr(4)
  
ELSE!inside
  
  elem_index_1=part_ele_location_ppr(1)
  elem_index_2=part_ele_location_ppr(2)
  elem_index_3=part_ele_location_ppr(3)
  elem_index_4=part_ele_location_ppr(4)
  
ENDIF

END SUBROUTINE 
    
    
!SUBROUTINE ParPosition_ppr_elem( node_corx, node_cory, elem_index_1, elem_index_2, elem_index_3, elem_index_4)
!  
!USE IFE_Data
!USE Particle_2D
!USE Field_2D
!USE IFE_MAIN_PARAM
!USE Domain_2D
!USE Constant_Variable_2D
!USE Object_Data_2D
!IMPLICIT NONE
!
!INTEGER :: i_part, repeat_refinement, count
!REAL(8) :: rxp_1, ryp_1, bottom_top_centre, left_right_centre, delta_1, delta_2
!REAL(8) :: rxp, ryp
!REAL(8) :: xmin, xmax, ymin, ymax
!INTEGER :: i_1, j_1, i, j, x
!INTEGER :: n_element_old, n_element_new, n_element_bro
!REAL(8) :: min_length, min_edge_length
!REAL(8), DIMENSION(:,:), ALLOCATABLE :: Sreal
!REAL(8), DIMENSION(:,:), ALLOCATABLE :: HE_particle
!CHARACTER*20	fname, sname, filename, num_n, num_e
!
!REAL(8), DIMENSION(:), ALLOCATABLE :: part_random_factor
!DOUBLE PRECISION :: ranum
!
!INTEGER, DIMENSION(:), ALLOCATABLE :: part_edge_count
!INTEGER :: edge_count
!INTEGER, DIMENSION(:), ALLOCATABLE :: edge_bro_element
!INTEGER :: edge_bro_count, edge_temp_count, edge_self_element, edge_self_element_min
!INTEGER :: edge_boundary
!
!!===============zyz===========================================
!REAL(8)                               :: node_corx, node_cory
!REAL(8)                               :: bias
!REAL(8)                               :: zmin, zmax
!INTEGER                               :: nnx, nny, nnz
!INTEGER                               :: N_part_tot_ppr
!INTEGER, DIMENSION(:), ALLOCATABLE    :: part_ele_location_ppr
!REAL(8), DIMENSION(:,:), ALLOCATABLE	:: part_ppr
!INTEGER, DIMENSION(:), ALLOCATABLE    :: elem_cg
!INTEGER, DIMENSION(:), ALLOCATABLE    :: error
!INTEGER                               :: elem_index_1, elem_index_2, elem_index_3, elem_index_4
!REAL(8)	                                    :: dimensions(2,2)
!!=============zyz=============================================
!
!!===============zyz vert_o赋值===============================================
!
!
!OPEN(1, ACTION = 'READ', FILE = 'mesh.inp')
!READ(1,*) (Vert_o(i), i = 1,3)
!READ(1,*)  xmax, ymax, zmax
!! nx, ny, nz							: Number of mesh points in x, y, z
!READ(1,*) nnx, nny, nnz
!
!! hx(1:3)								: Grid resolution in each direction
!READ(1,*) (hx(i), i = 1, 2)
!READ(1,*) hz
!CLOSE(1)
!! xmin, ymin, xmax, ymax				: Domain limits
!!Vert_o(1) = xmin
!!Vert_o(2) = ymin
!!Vert_o(3) = zmin!zyz 这里为什么要重新赋值不理解
!!==============================================================zyz============
!xmin = Vert_o(1)
!ymin = Vert_o(2)
!!==============================================================zyz============
!nx = nnx
!ny = nny
!nz = nnz
!
!hxi = One/hx
!hzi = One/hz
!! Origin for the PIC domain
!Vert_o(1:2) = Vert_o(1:2) - hx
!Vert_o(3) = Vert_o(3) - hz
!!===============zyz vert_o赋值===============================================
!ALLOCATE(Sreal(2, SIZE(P_average,2)))
!ALLOCATE(HE_particle(5, SIZE(HE,2)))
!ALLOCATE(part_edge_count(9))
!ALLOCATE(edge_bro_element(4))
!
!Sreal(1:2, 1:SIZE(P_average,2)) = Zero
!HE_particle(1:5, 1:SIZE(HE,2)) = Zero
!edge_bro_element(1:4) = 0
!edge_boundary = 0
!
!!repeat_refinement = SIZE(HT_flag_MESH,1) - 1
!!===============================LY ADD FOR DEBUG, 2021-12-23===================================
!!=================================找到四个点进行粒子定位 ZYZ
!!=================================!=================================
!N_part_tot_ppr = 4
!ALLOCATE(part_ele_location_ppr(N_part_tot_ppr))
!ALLOCATE(part_ppr(4,2))
!part_ppr= 0.D0
!part_ele_location_ppr=0.D0
!
!OPEN(1,ACTION='READ',FILE='IDG_inf.inp')
!    READ(1, *)
!    READ(1, *) repeat_refinement
!    READ(1, *) 
!CLOSE(1)
!bias=hx(1)/(2**(repeat_refinement+2))
!!IF (node_corx = xmin .AND. node_cory = ymin) THEN!左下角
!!  elem_index_1 = 1
!!  elem_index_2 = 1 
!!  elem_index_3 = 1 
!!  elem_index_4 = 1
!!ELSEIF (node_corx = xmin .AND. node_cory = ymax) THEN!左上角
!!  elem_index_1 = nny-1
!!  elem_index_2 = nny-1 
!!  elem_index_3 = nny-1 
!!  elem_index_4 = nny-1
!!ELSEIF (node_corx = xmin .AND. node_cory = ymax) THEN!右上角
!!  elem_index_1 = (nnx-1)*(nny-1)
!!  elem_index_2 = (nnx-1)*(nny-1) 
!!  elem_index_3 = (nnx-1)*(nny-1) 
!!  elem_index_4 = (nnx-1)*(nny-1)  
!!ELSEIF (node_corx = xmin .AND. node_cory = ymax) THEN!右下角
!!  elem_index_1 = (nnx-1)*(nny-2) + 1
!!  elem_index_2 = (nnx-1)*(nny-2) + 1 
!!  elem_index_3 = (nnx-1)*(nny-2) + 1 
!!  elem_index_4 = (nnx-1)*(nny-2) + 1  
!IF (node_corx == xmin .AND. node_cory == ymin) THEN!左下角
!  
!  part_ppr(1,1) = xmin + bias
!  part_ppr(1,2) = ymin + bias
!  part_ppr(2,1) = xmin
!  part_ppr(2,2) = ymin
!  part_ppr(3,1) = xmin
!  part_ppr(3,2) = ymin
!  part_ppr(4,1) = xmin
!  part_ppr(4,2) = ymin
!  
!ELSEIF (node_corx == xmin .AND. node_cory == ymax) THEN!左上角
!  
!  part_ppr(1,1) = xmin + bias
!  part_ppr(1,2) = ymax - bias
!  part_ppr(2,1) = xmin
!  part_ppr(2,2) = ymin
!  part_ppr(3,1) = xmin
!  part_ppr(3,2) = ymin
!  part_ppr(4,1) = xmin
!  part_ppr(4,2) = ymin
!  
!ELSEIF (node_corx == xmax .AND. node_cory == ymax) THEN!右上角
!  
!  part_ppr(1,1) = xmax - bias
!  part_ppr(1,2) = ymax - bias
!  part_ppr(2,1) = xmin
!  part_ppr(2,2) = ymin
!  part_ppr(3,1) = xmin
!  part_ppr(3,2) = ymin
!  part_ppr(4,1) = xmin
!  part_ppr(4,2) = ymin
!  
!ELSEIF (node_corx == xmax .AND. node_cory == ymin) THEN!右下角
! 
!  part_ppr(1,1) = xmax - bias
!  part_ppr(1,2) = ymin + bias
!  part_ppr(2,1) = xmax
!  part_ppr(2,2) = ymin
!  part_ppr(3,1) = xmax
!  part_ppr(3,2) = ymin
!  part_ppr(4,1) = xmax
!  part_ppr(4,2) = ymin
!  
!ELSEIF (node_corx == xmin) THEN!左边界
!  
!  part_ppr(2,1) = node_corx + bias
!  part_ppr(2,2) = node_cory - bias
!  part_ppr(3,1) = node_corx + bias
!  part_ppr(3,2) = node_cory + bias
!  
!  part_ppr(1,1) = xmin
!  part_ppr(1,2) = ymin
!  part_ppr(4,1) = xmin
!  part_ppr(4,2) = ymin
!  
!ELSEIF(node_cory == ymax)THEN!上边界
!  
!  part_ppr(1,1) = node_corx - bias
!  part_ppr(1,2) = node_cory - bias
!  part_ppr(2,1) = node_corx + bias
!  part_ppr(2,2) = node_cory - bias
!
!  part_ppr(3,1) = xmin
!  part_ppr(3,2) = ymin
!  part_ppr(4,1) = xmin
!  part_ppr(4,2) = ymin
!  
!ELSEIF(node_corx == xmax) THEN!右边界
!  
!  part_ppr(1,1) = node_corx - bias
!  part_ppr(1,2) = node_cory - bias
!  part_ppr(4,1) = node_corx - bias
!  part_ppr(4,2) = node_cory + bias
!  
!  part_ppr(2,1) = xmin
!  part_ppr(2,2) = ymin
!  part_ppr(3,1) = xmin
!  part_ppr(3,2) = ymin
!  
!ELSEIF(node_cory == ymin)THEN!下边界
!  
!  part_ppr(3,1) = node_corx + bias
!  part_ppr(3,2) = node_cory + bias
!  part_ppr(4,1) = node_corx - bias
!  part_ppr(4,2) = node_cory + bias
!  
!  part_ppr(1,1) = xmin
!  part_ppr(1,2) = ymin
!  part_ppr(2,1) = xmin
!  part_ppr(2,2) = ymin
!  
!ELSE!内部点
!  
!  part_ppr(1,1) = node_corx - bias
!  part_ppr(1,2) = node_cory - bias
!
!  part_ppr(2,1) = node_corx + bias
!  part_ppr(2,2) = node_cory - bias
!
!  part_ppr(3,1) = node_corx + bias
!  part_ppr(3,2) = node_cory + bias
!
!  part_ppr(4,1) = node_corx - bias
!  part_ppr(4,2) = node_cory + bias
!ENDIF
!!=================================!=================================zyz
!
!
!!==========Part1 : Particle location==========
!!Before we start distribute every particle into any mesh nodes, we need to locate particle.
!
!IF(repeat_refinement == 0)THEN
!  
!  DO i_part = 1, N_part_tot_ppr
!
!! 1st find out the species
!	
!
!	rxp = (part_ppr(i_part,1) - Vert_o(1))*hxi(1)
!	i = rxp
!	
!
!	ryp = (part_ppr(i_part,2) - Vert_o(2))*hxi(2)
!	j = ryp
!  
!  n_element_old = (i-1)*(nny-1) + j
!  
!  part_ele_location_ppr(i_part) = n_element_old
!  ENDDO
!  
!  
!  
!	
!  
!ELSE 
!  DO i_part = 1, N_part_tot_ppr
!  
!  count = 1   !The variable 'count' is record researching layer of mesh.
!  part_edge_count(1:9) = 0   !The variable 'part_edge_count' is record intersection number between particle and edge.
!  edge_count = 2    !The variable 'edge_count'
!  
!  !Firstly, we need to locate particle on the initial mesh: 1st mesh(no refinement).
!  rxp_1 = (part_ppr(i_part,1) - Vert_o(1)) * hxi(1)
!  i_1 = rxp_1
!  
!  ryp_1 = (part_ppr(i_part,2) - Vert_o(2)) * hxi(2)
!  j_1 = ryp_1
!  
!  IF (j_1 == ny) THEN
!    j_1 = j_1 - 1
!  END IF
!  
!  IF (i_1 == nx) THEN
!    i_1 = i_1 - 1
!  END IF
!  
!  n_element_old = (j_1 + (i_1 - 1)*(ny-1))
!  
!  DO j = 1, SIZE(HE,2)
!    HE_particle(1, j) = SQRT((part_ppr(i_part,1)-HP(1,HE(1,j)))**2 + (part_ppr(i_part,2)-HP(2,HE(1,j)))**2)
!    HE_particle(2, j) = SQRT((part_ppr(i_part,1)-HP(1,HE(2,j)))**2 + (part_ppr(i_part,2)-HP(2,HE(2,j)))**2)
!    HE_particle(3, j) = HE_particle(1, j) + HE_particle(2, j)
!    HE_particle(4, j) = HE(5, j)
!    HE_particle(5, j) = SQRT((HP(1,HE(2,j))-HP(1,HE(1,j)))**2 + (HP(2,HE(2,j))-HP(2,HE(1,j)))**2)
!    
!    IF (HE_particle(3,j) == HE_particle(5,j)) THEN
!      part_edge_count(1) = part_edge_count(1) + 1
!      part_edge_count(edge_count)  = j
!      edge_count = edge_count + 1
!    END IF
!  END DO
!  !HE_particle(3,:) store distance between the particle and node of every edge.
!  !HE_particle(5,:) store length of every edge.
!  
!  min_length = MINVAL(HE_particle(3, 1:SIZE(HE,2)))
!  
!  !========================particle locate on CG and DG intersect edge: CG, DG coupling===========================
!  IF (part_edge_count(1) == 1) THEN
!  !The particle locate the edge between CG(0.5) and DG(0.5).√
!  
!    IF (HE(6, part_edge_count(2)) == 0) THEN    !There denote the DG edge is exterior boundary edge. We need protect program from array.
!    
!      n_element_old = HE(5, part_edge_count(2))
!      count = repeat_refinement + 1
!      
!    ELSEIF (HE(6, part_edge_count(2)) /= 0) THEN  !There denote the DG edge is interior edge.
!    
!      IF (HT_MESH(repeat_refinement+1, 5, HE(5,part_edge_count(2)))/=0 .AND. &
!          HT_MESH(repeat_refinement+1, 5, HE(6,part_edge_count(2)))==0) THEN
!        count = repeat_refinement + 1
!        n_element_old = HE(5, part_edge_count(2))   !The 0.5 particle distribute DG element
!        n_element_bro = HE(6, part_edge_count(2))   !The 0.5 particle distribute CG element
!
!        !===============================LY ADD FOR DEBUG, 2021-12-23===================================
!
!        CYCLE
!
!      ELSE
!        WRITE(6,*) 'Program error, pelase check : part_edge_count(1) == 1!'
!        STOP
!      END IF
!    END IF
!  
!    
!  ELSEIF (part_edge_count(1) == 2) THEN
!  !The particle locate the edge between CG(0.0) and DG(1.0).√
!  !The particle locate the edge between DG(0.5) and DG(0.5).√
!    
!    IF (HE(5, part_edge_count(2)) == HE(5, part_edge_count(3))) THEN
!    !The particle locate the edge between CG(0.0) and DG(1.0).(3CG and 1DG, the particle locate the corner node)
!      n_element_old = HE(5,part_edge_count(2))
!      count = repeat_refinement + 1
!      
!    ELSEIF (HE(5, part_edge_count(2)) /= HE(5, part_edge_count(3))) THEN
!    !The particle locate the edge between DG(0.5) and DG(0.5).  
!      count = repeat_refinement + 1
!      n_element_old = HE(5, part_edge_count(2))   !The 0.5 particle distribute DG element
!      n_element_bro = HE(5, part_edge_count(3))   !The 0.5 particle distribute DG element
!
!      !===============================LY ADD FOR DEBUG, 2021-12-23===================================
!      
!      CYCLE
!    END IF
!  
!    
!  ELSEIF (part_edge_count(1) == 4) THEN
!  !The particle locate the edge between CG(0.0) and DG(1.0). (condition 1)√
!  !The particle locate the edge between CG(0.5) and DG(0.5). (condition 2)√
!  
!    edge_boundary = 0
!    DO j = 2, 5
!      IF (HE(6, part_edge_count(j)) == 0) THEN
!        edge_boundary = edge_boundary + 1
!      END IF
!    END DO
!    
!    IF (edge_boundary == 2) THEN
!      n_element_old = HE(5, part_edge_count(2))   !There is to pretect array from outing of bound.
!      count = repeat_refinement + 1
!    ELSE
!    
!      edge_bro_count = 0
!      DO j = 2, 5
!        IF (HT_MESH(repeat_refinement+1, 5, HE(5, part_edge_count(j)))>0 .AND. &
!            HT_MESH(repeat_refinement+1, 5, HE(6, part_edge_count(j)))==0) THEN
!          edge_bro_element(1+edge_bro_count) = HE(6, part_edge_count(j))
!          edge_bro_count = edge_bro_count + 1
!        END IF
!      END DO
!
!      IF (edge_bro_count==2 .AND. edge_bro_element(1)==edge_bro_element(2)) THEN
!        !The particle locate the edge between CG(0.5) and DG(0.5). (condition 2: 2CG and 2DG--opposite, CG element is same)
!        count = repeat_refinement + 1
!        n_element_old = HE(5, part_edge_count(2))   !The 0.5 particle distribute DG element
!        n_element_bro = edge_bro_element(1)   !The 0.5 particle distribute CG element
!        !===============================LY ADD FOR DEBUG, 2021-12-23===================================
!
!        CYCLE
!
!      ELSEIF (edge_bro_count==2 .AND. edge_bro_element(1)/=edge_bro_element(2)) THEN
!        !The particle locate the edge between CG(0.0) and DG(1.0). (condition 1: 2CG and 2DG--opposite, CG element is different)
!        n_element_old = HE(5, part_edge_count(2))
!        count = repeat_refinement + 1
!
!      ELSEIF (edge_bro_count==4) THEN
!        !The particle locate the edge between CG(0.0) and DG(1.0). (condition 1: 2CG and 2DG--diagonal, CG element is different)
!        n_element_old = HE(5, part_edge_count(2))
!        count = repeat_refinement + 1
!      END IF
!    END IF
!  
!    
!  ELSEIF (part_edge_count(1) == 6) THEN
!  !The particle locate the edge between CG(1.0) and DG(0.0).(1CG and 3DG, the particle locate the corner node)√
!  !The particle locate the edge between DG(0.5) and DG(0.5).(2DG and 1DG, the particle locate the hand node)√
!  
!    edge_bro_count = 0
!    edge_temp_count = 0
!    edge_self_element = 0    
!    edge_self_element_min = MINVAL(HT_MESH(repeat_refinement+1, 5, HE(5,part_edge_count(2:7))))
!    
!    DO j = 2, 7
!      IF (HT_MESH(repeat_refinement+1, 5, HE(5,part_edge_count(j)))>0 .AND. &
!          HT_MESH(repeat_refinement+1, 5, HE(6,part_edge_count(j)))==0) THEN
!        edge_bro_element(1+edge_bro_count) = HE(6, part_edge_count(j))
!        edge_bro_count = edge_bro_count + 1
!      ELSEIF (HT_MESH(repeat_refinement+1, 5, HE(5,part_edge_count(j)))>0 .AND. &
!              HT_MESH(repeat_refinement+1, 5, HE(6,part_edge_count(j)))>0) THEN
!        edge_temp_count = edge_temp_count + 1
!        IF (edge_self_element_min == HT_MESH(repeat_refinement+1, 5, HE(5,part_edge_count(j)))) THEN
!          edge_self_element = part_edge_count(j)
!        END IF
!      END IF
!    END DO
!    
!    IF (edge_bro_count==2 .AND. edge_bro_element(1)==edge_bro_element(2)) THEN
!    !The particle locate the edge between CG(1.0) and DG(0.0).(1CG and 3DG, the particle locate the corner node)
!      n_element_old = edge_bro_element(1)
!      count = repeat_refinement + 1
!      
!    ELSEIF (edge_temp_count == 6) THEN
!    !The particle locate the edge between DG(0.5) and DG(0.5).(2DG and 1DG, the particle locate the hand node)
!      count = repeat_refinement + 1
!      n_element_old = HE(5, edge_self_element)
!      n_element_bro = HE(6, edge_self_element)   !The 0.5 particle distribute DG element
!      !===============================LY ADD FOR DEBUG, 2021-12-23===================================
!      
!      CYCLE
!    END IF
!    
!    
!  ELSEIF (part_edge_count(1) == 8) THEN
!  !The particle locate the edge between DG(1.0) and DG(0.0).(4DG, the particle locate the corner node)√
!    n_element_old = HE(5, part_edge_count(2))
!    count = repeat_refinement + 1
!  END IF
! 
!  !========================particle locate on CG and DG intersect edge: Random===========================
!  
!  DO WHILE(HT_flag_MESH(count, 2, n_element_old) == 1)
!    !There denote we need continually research particle location.
!    bottom_top_centre = HP_MESH(count+1, 2, HT_MESH(count+1, 3, HT_flag_MESH(count, 3, n_element_old)))
!    left_right_centre = HP_MESH(count+1, 1, HT_MESH(count+1, 3, HT_flag_MESH(count, 3, n_element_old)))
!    
!    delta_1 = part_ppr(i_part,2) - bottom_top_centre
!    delta_2 = part_ppr(i_part,1) - left_right_centre
!    
!  
!    
!    IF (delta_1 > 0) THEN
!      IF (delta_2 > 0) THEN       !right-upper children  element
!        n_element_new = HT_flag_MESH(count, 6, n_element_old)
!      ELSEIF (delta_2 < 0) THEN   !left-upper children element
!        n_element_new = HT_flag_MESH(count, 4, n_element_old)
!      ELSEIF (delta_2 == 0) THEN  !left-right_centre line
!        !when particle locate element-element intersect line, we directly traverse particle and every edge.
!        DO j = 1, SIZE(HE,2)
!          IF (HE_particle(3,j) == min_length) THEN  
!            IF (HT_MESH(repeat_refinement+1, 5, HE(5,j)) > HT_MESH(repeat_refinement+1, 5, HE(6,j))) THEN
!              n_element_old = HE(5,j)
!            ELSEIF (HT_MESH(repeat_refinement+1, 5, HE(5,j)) < HT_MESH(repeat_refinement+1, 5, HE(6,j))) THEN
!              n_element_old = HE(6,j)
!            ELSEIF (HT_MESH(repeat_refinement+1, 5, HE(5,j)) == HT_MESH(repeat_refinement+1, 5, HE(6,j))) THEN
!              n_element_old = HE(5,j)
!            END IF
!          END IF
!        END DO
!        count = repeat_refinement + 1
!        EXIT
!        !Here, we directly determine the final element of particle location, so we don't need traverse other mesh information.
!      END IF
!      
!    ELSEIF (delta_1 < 0) THEN
!      IF (delta_2 > 0) THEN       !right-lower children  element
!        n_element_new = HT_flag_MESH(count, 5, n_element_old)
!      ELSEIF (delta_2 < 0) THEN   !left-lower children element
!        n_element_new = HT_flag_MESH(count, 3, n_element_old)
!      ELSEIF (delta_2 == 0) THEN  !left_right_centre line
!        !when particle locate element-element intersect line, we directly traverse particle and every edge.
!        DO j = 1, SIZE(HE,2)
!          IF (HE_particle(3,j) == min_length) THEN  
!            IF (HT_MESH(repeat_refinement+1, 5, HE(5,j)) > HT_MESH(repeat_refinement+1, 5, HE(6,j))) THEN
!              n_element_old = HE(5,j)
!            ELSEIF (HT_MESH(repeat_refinement+1, 5, HE(5,j)) < HT_MESH(repeat_refinement+1, 5, HE(6,j))) THEN
!              n_element_old = HE(6,j)
!            ELSEIF (HT_MESH(repeat_refinement+1, 5, HE(5,j)) == HT_MESH(repeat_refinement+1, 5, HE(6,j))) THEN
!              n_element_old = HE(5,j)
!            END IF
!          END IF
!        END DO
!        count = repeat_refinement + 1
!        EXIT
!        !Here, we directly determine the final element of particle location, so we don't need traverse other mesh information.
!      END IF
!      
!    ELSEIF (delta_1 == 0) THEN    !particle locate bottom_top_centre
!      IF (delta_2 == 0) THEN
!        n_element_new = HT_flag_MESH(count, 3, n_element_old)
!      ELSE
!        !when particle locate element-element intersect line, we directly traverse particle and every edge.
!        DO j = 1, SIZE(HE,2)
!          IF (HE_particle(3,j) == min_length) THEN  
!            IF (HT_MESH(repeat_refinement+1, 5, HE(5,j)) > HT_MESH(repeat_refinement+1, 5, HE(6,j))) THEN
!              n_element_old = HE(5,j)
!            ELSEIF (HT_MESH(repeat_refinement+1, 5, HE(5,j)) < HT_MESH(repeat_refinement+1, 5, HE(6,j))) THEN
!              n_element_old = HE(6,j)
!            ELSEIF (HT_MESH(repeat_refinement+1, 5, HE(5,j)) == HT_MESH(repeat_refinement+1, 5, HE(6,j))) THEN
!              n_element_old = HE(5,j)
!            END IF
!          END IF
!        END DO
!        count = repeat_refinement + 1
!        EXIT
!        !Here, we directly determine the final element of particle location, so we don't need traverse other mesh information.
!      END IF
!    END IF
!    
!    count = count + 1
!    n_element_old = n_element_new
!    
!  END DO
!  
!  DO WHILE (count /= (repeat_refinement+1))
!    n_element_new = HT_flag_MESH(count, 3, n_element_old)
!    count = count + 1
!    n_element_old = n_element_new
!  END DO
!  
!  n_element_old = HT_flag_MESH(count, 3, n_element_old)
!  
!  part_ele_location_ppr(i_part) = n_element_old
!  
!END DO
!ENDIF
!
!!==========Part1 : Particle location==========
!!=================================!=================================zyz
!!=================================!=================================各种类型的点汇总
!
!IF (node_corx == xmin .AND. node_cory == ymin) THEN!左下角
!    
!  elem_index_1 = 1
!  elem_index_2 = 0
!  elem_index_3 = 0 
!  elem_index_4 = 0
!  
!ELSEIF (node_corx == xmin .AND. node_cory == ymax) THEN!左上角
!    
!  elem_index_1 = part_ele_location_ppr(1)
!  elem_index_2 = 0 
!  elem_index_3 = 0 
!  elem_index_4 = 0
! 
!ELSEIF (node_corx == xmax .AND. node_cory == ymax) THEN!右上角
!  
!  elem_index_1 = part_ele_location_ppr(1)
!  elem_index_2 = 0 
!  elem_index_3 = 0 
!  elem_index_4 = 0
! 
!ELSEIF (node_corx == xmax .AND. node_cory == ymin) THEN!右下角
!  
!  elem_index_1 = part_ele_location_ppr(1)
!  elem_index_2 = 0
!  elem_index_3 = 0
!  elem_index_4 = 0
!  
!ELSEIF (node_corx == xmin) THEN!左边界
!  
!  elem_index_1 = 0
!  elem_index_4 = 0
!  elem_index_2=part_ele_location_ppr(2)
!  elem_index_3=part_ele_location_ppr(3)
!  
!ELSEIF(node_cory == ymax)THEN!上边界
! 
!  elem_index_3 = 0
!  elem_index_4 = 0
!  elem_index_1=part_ele_location_ppr(1)
!  elem_index_2=part_ele_location_ppr(2)
!  
!ELSEIF(node_corx == xmax) THEN!右边界
! 
!  elem_index_2 = 0
!  elem_index_3 = 0
!  elem_index_1=part_ele_location_ppr(1)
!  elem_index_4=part_ele_location_ppr(4)
!  
!ELSEIF(node_cory == ymin)THEN!下边界
!  
!  elem_index_1 = 0
!  elem_index_2 = 0
!  elem_index_3=part_ele_location_ppr(3)
!  elem_index_4=part_ele_location_ppr(4)
!  
!ELSE!内部点
!  
!  elem_index_1=part_ele_location_ppr(1)
!  elem_index_2=part_ele_location_ppr(2)
!  elem_index_3=part_ele_location_ppr(3)
!  elem_index_4=part_ele_location_ppr(4)
!  
!ENDIF
!
!END SUBROUTINE