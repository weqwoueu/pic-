SUBROUTINE GetParChg_2D(delta)

! Updated:	11/28/2004 06:00 PM
! Purpose:	Get the cahrge desnity for particles.
!			Changed affp defination. change GetParChg accordingly
!			get the charge density for each species and the total charge density
!			modified to account for hx,hy,hz .NE. 1

!USE PIC_Main_Param_2D
USE Domain_2D
USE Field_2D
USE Particle_2D
USE Cell_Volume_Data
USE IFE_Data
USE Constant_Variable_2D
Use ModuleMCCInterface,ONLY:ControlFlowGlobal, ParticleGlobal,JtoeV   !$ ab.ZWZ for using JW's particle data structure
Use IFE_INTERFACE, Only: ParticlePositioning
!Use Constants
IMPLICIT NONE

INTEGER		delta

INTEGER		i, j, nden_var, i_part, isp, ii, jj, l
REAL(8)		rxp, ryp, dx, dy, xcellmdx, ycellmdy, R1, R2, den
INTEGER		ntemp
INTEGER   n, n_element, n_node
INTEGER, DIMENSION(4)    ::  node_index_el
REAL(8)        rho_temp
REAL(8), DIMENSION(2,4)  ::  VertX_Cell
REAL(8)                    ::  dist1, dist2, dist3, dist4, dist_total
REAL(8)                    ::  rho_temp1, rho_temp2, rho_temp3, rho_temp4, rho_temp_total
REAL(8)                    ::  x_part, y_part, x_rec, y_rec, Cell_volume_bjw

!$ =============== mb.ZWZ 2021/7/11 ========================== \\
REAL(8) :: P1, P2, P3, P4 
REAL(8) :: R
REAL(8) :: YFACTOR
REAL(8) :: BETA  = 1.!$ 参靠操慧珺论文，添加Larson的修正因子
!$ =============== mb.ZWZ 2021/7/11 ========================== //
REAL(8) :: V2 !$ ab.ZWZ 2022/2/29 for output energy distribution
REAL(8) :: Energy !$ ab.ZWZ 2022/5/26 for output energy distribution

!=========LY modification for Multi-Layer-Grid, 2022-7-25=========
Real(8) :: part_x, part_y
Integer :: traverse_flag, n_element_old, n_element_bro
Integer :: n_element_bro_initial, number_ele_column
Real(8) :: part_ele_left, part_ele_right, part_ele_bottom, part_ele_top
Real(8) :: left_right_centre, bottom_top_centre, delta_1, delta_2
Integer :: part_edge_count(9)
Integer :: edge_count
Integer :: num, k
Real(8) :: hx_partition, hy_partition
Real(8) :: W1, W2, W3, W4
Integer :: edge_boundary
Integer :: edge_bro_count, edge_temp_count, edge_self_element, edge_self_element_min
Integer :: edge_bro_element(4)
Integer :: temp_count
Integer, Dimension(:), Allocatable :: EdgeArray
Real(8), Dimension(:,:), Allocatable :: HE_particle
Real(8), Dimension(:,:,:), Allocatable :: rho_s_sum, rho_sum, rho_p_sum
Real(8), Dimension(:,:,:), Allocatable :: Ek_tot_sum, Ek_tot_mesh_sum, rho_p_mesh_sum
Real(8), Dimension(:,:,:), Allocatable :: Ek_s_mesh_sum
Real(8), Dimension(:,:,:), Allocatable :: rho_n_mesh_sum
Real(8), Dimension(:,:), Allocatable :: rho_test
Real(8), Dimension(:,:,:), Allocatable :: rho_s_test, rho_p_test
Real(8), Dimension(:,:,:), Allocatable :: Ek_tot_test
!=========LY modification for Multi-Layer-Grid, 2022-7-25=========

!=========LY modification for Multi-Layer-Grid, 2022-7-25=========
!======Part 0: Initialization======
Write(6,*) 'GetParChg_2D'
Write(6,*) 'delta = ', delta
number_ele_column = ny - 1

Deallocate(rho, rho_s, rho_p, Ek_tot)
Allocate(rho(Size(P_average,2), 1))
Allocate(rho_s(4, Size(HT,2), ispe_tot))
Allocate(rho_p(4, Size(HT,2), ispe_tot))
Allocate(Ek_tot(4, Size(HT,2), ispe_tot))

rho = 0.0
rho_s = 0.0
rho_p = 0.0
Ek_tot = 0.0

If (.Not. Allocated(HE_particle))   Allocate(HE_particle(5, Size(HE,2)))
If (.Not. Allocated(rho_s_sum))     Allocate(rho_s_sum(1, Size(HP,2), ispe_tot))
If (.Not. Allocated(rho_p_sum))     Allocate(rho_p_sum(1, Size(HP,2), ispe_tot))
If (.Not. Allocated(Ek_tot_sum))    Allocate(Ek_tot_sum(1, Size(HP,2), ispe_tot))
HE_particle = 0.0
rho_s_sum = 0.0
rho_p_sum = 0.0
Ek_tot_sum = 0.0

If (.Not. Allocated(rho_sum))           Allocate(rho_sum(Size(P_average,2), 1, ispe_tot))
If (.Not. Allocated(rho_p_mesh_sum))    Allocate(rho_p_mesh_sum(Size(P_average,2), 1, ispe_tot))
If (.Not. Allocated(rho_n_mesh_sum))    Allocate(rho_n_mesh_sum(Size(P_average,2), 1, ispe_tot))
If (.Not. Allocated(Ek_tot_mesh_sum))   Allocate(Ek_tot_mesh_sum(Size(P_average,2), 1, ispe_tot))
If (.Not. Allocated(Ek_s_mesh_sum))     Allocate(Ek_s_mesh_sum(Size(P_average,2), 1, ispe_tot))
rho_sum = 0.0
rho_p_mesh_sum = 0.0
rho_n_mesh_sum = 0.0
Ek_tot_mesh_sum = 0.0
Ek_s_mesh_sum = 0.0

If (.Not. Allocated(rho_test))      Allocate(rho_test(Size(HP,2), 1))
If (.Not. Allocated(rho_s_test))    Allocate(rho_s_test(Size(HP,2), 1, ispe_tot))
If (.Not. Allocated(rho_p_test))    Allocate(rho_p_test(Size(HP,2), 1, ispe_tot))
If (.Not. Allocated(Ek_tot_test))   Allocate(Ek_tot_test(Size(HP,2), 1, ispe_tot))
rho_test = 0.0
rho_s_test = 0.0
rho_p_test = 0.0
rho_n = 0.0
Ek_tot_test = 0.0
Ek_s = 0.0


Mapxmin = hx(1) / (2**repeat_refinement)
Mapymin = hx(2) / (2**repeat_refinement)

!==========Part 1 : Particle location and charge distribution==========
Do isp=0, ControlFlowGlobal%Ns
  Do i_part = 1, ParticleGlobal(isp)%NPar
    
    part_x = ParticleGlobal(isp)%PO(i_part)%X
    part_y = ParticleGlobal(isp)%PO(i_part)%Y
    
    Energy = ParticleGlobal(isp)%PO(i_part)%Energy(ParticleGlobal(isp)%Mass,ParticleGlobal(isp)%VFactor)/JtoeV
    !Ek_one(isp,i_part)=Energy
    
    traverse_flag = 0
    
    If (Allocated(EdgeArray)) Then
      Deallocate(EdgeArray)
    End If
      
    rxp = (part_x - Vert_o(1)) * hxi(1)
    i = INT(rxp)
    
    ryp = (part_y - Vert_o(2)) * hxi(2)
    j = INT(ryp)
    
    If (j == ny) Then
      j = j - 1
    End If
    
    If (i == nx) Then
      i = i - 1
    End If
    
    n_element_old = (j + (i-1)*(ny-1))
    
    !Write(*,*) Mod(rxp, Mapxmin), Mod(ryp, Mapymin)
    
    If ( Mod(rxp, Mapxmin) > (1D-10) .And. Mod(ryp, Mapymin) > (1D-10) ) Then
        traverse_flag = 0
        Call InjectPositioning(part_x, part_y, n_element_old)
    Else
        Call ParticlePositioning(part_x, part_y, n_element_old, traverse_flag, EdgeArray)
    End If
    
    
    Call ParticleGlobal(isp)%PO(i_part)%ParLocate(n_element_old)
    
    IF (traverse_flag == 1) THEN
      !There denote the particle locate on edge. We need to traverse between particle and every edge.
      part_edge_count(1:9) = 0   !The variable 'part_edge_count' is record intersection number between particle and edge.
      edge_count = 2    !The variable 'edge_count' is temp counter.
      
      If (Allocated(HE_particle)) Then
        Deallocate(HE_particle)
      End If
      Allocate(HE_particle(5,Size(EdgeArray)))
      HE_particle = 0.0

      Do j = 1, Size(EdgeArray)
        HE_particle(1,j) = DSQRT((part_x-HP(1,HE(1,EdgeArray(j))))**2 + &
                                 (part_y-HP(2,HE(1,EdgeArray(j))))**2)
        HE_particle(2,j) = DSQRT((part_x-HP(1,HE(2,EdgeArray(j))))**2 + &
                                 (part_y-HP(2,HE(2,EdgeArray(j))))**2)
        HE_particle(5,j) = DSQRT((HP(1,HE(2,EdgeArray(j)))-HP(1,HE(1,EdgeArray(j))))**2 + &
                                 (HP(2,HE(2,EdgeArray(j)))-HP(2,HE(1,EdgeArray(j))))**2)

        HE_particle(3,j) = HE_particle(1,j) + HE_particle(2,j)
        HE_particle(4,j) = HE(5,EdgeArray(j))
        !HE_particle(3,:) store distance between the particle and node of every edge.
        !HE_particle(5,:) store length of every edge.

        If (ABS(HE_particle(3,j)-HE_particle(5,j))<SmallValue) Then
          part_edge_count(1) = part_edge_count(1) + 1
          part_edge_count(edge_count)  = EdgeArray(j)
          edge_count = edge_count + 1
        End If
      End Do
      !=========LY modification for Speeding Particle Positioning, 2022-6-22=========

      !========================particle locate on CG and DG intersect edge: CG, DG coupling===========================
      IF (part_edge_count(1) == 1) THEN
        !The particle locate the edge between CG(0.5) and DG(0.5).

        IF (HE(6, part_edge_count(2)) == 0) THEN    !There denote the DG edge is exterior boundary edge. We need protect program from array.

          n_element_old = HE(5, part_edge_count(2))

        ELSEIF (HE(6, part_edge_count(2)) /= 0) THEN  !There denote the DG edge is interior edge.

          IF (HT(5, HE(5,part_edge_count(2)))/=0 .AND. HT(5, HE(6,part_edge_count(2)))==0) THEN

            n_element_old = HE(5, part_edge_count(2))   !The 0.5 particle distribute DG element

            hx_partition = HP(1, HT(2, n_element_old)) - HP(1, HT(1, n_element_old))
            hy_partition = HP(2, HT(4, n_element_old)) - HP(2, HT(1, n_element_old))

            dx = (part_x - HP(1, HT(1, n_element_old))) / hx_partition
            dy = (part_y - HP(2, HT(1, n_element_old))) / hy_partition
            xcellmdx = 1.0 - dx
            ycellmdy = 1.0 - dy

            IF (delta == 0) THEN
              W1 = xcellmdx * ycellmdy  !local 1
              W2 = dx       * ycellmdy  !local 2
              W3 = dx       * dy        !local 3
              W4 = xcellmdx * dy        !local 4
            ELSEIF (delta == 1) THEN

              IF (ABS(HP(2,HT(1,n_element_old))-dymin) < SmallValue) THEN
                !There denote the element locate on the bottom edge, it equal to 'j==1' on Program 'hall_zwz'
                BETA = 0.75
              ELSE
                BETA = 1.0
              END IF

              R1 = dymin + HP(2, HT(1, n_element_old))
              R2 = dymin + HP(2, HT(4, n_element_old))
              R = part_y
              den = R2*R2 - R1*R1

              W1 = xcellmdx * (R2*R2-R*R) / den * BETA  !local 1
              W2 = dx       * (R2*R2-R*R) / den * BETA  !local 2
              W3 = dx       * (R*R-R1*R1) / den         !local 3
              W4 = xcellmdx * (R*R-R1*R1) / den         !local 4
            END IF
            
            rho_p(1, n_element_old, isp+1) = rho_p(1, n_element_old, isp+1) + W1 * 0.5
            rho_p(2, n_element_old, isp+1) = rho_p(2, n_element_old, isp+1) + W2 * 0.5
            rho_p(3, n_element_old, isp+1) = rho_p(3, n_element_old, isp+1) + W3 * 0.5
            rho_p(4, n_element_old, isp+1) = rho_p(4, n_element_old, isp+1) + W4 * 0.5

            If (ParticleGlobal(isp)%UnequalWeightFlag) Then
              W1 = W1 * ParticleGlobal(isp)%PO(i_part)%WQ
              W2 = W2 * ParticleGlobal(isp)%PO(i_part)%WQ
              W3 = W3 * ParticleGlobal(isp)%PO(i_part)%WQ
              W4 = W4 * ParticleGlobal(isp)%PO(i_part)%WQ
            End If
            
            rho_s(1, n_element_old, isp+1) = rho_s(1, n_element_old, isp+1) + W1 * 0.5
            rho_s(2, n_element_old, isp+1) = rho_s(2, n_element_old, isp+1) + W2 * 0.5
            rho_s(3, n_element_old, isp+1) = rho_s(3, n_element_old, isp+1) + W3 * 0.5
            rho_s(4, n_element_old, isp+1) = rho_s(4, n_element_old, isp+1) + W4 * 0.5
            
            Ek_tot(1, n_element_old, isp+1) = Ek_tot(1, n_element_old, isp+1) + W1 * Energy * 0.5
            Ek_tot(2, n_element_old, isp+1) = Ek_tot(2, n_element_old, isp+1) + W2 * Energy * 0.5
            Ek_tot(3, n_element_old, isp+1) = Ek_tot(3, n_element_old, isp+1) + W3 * Energy * 0.5
            Ek_tot(4, n_element_old, isp+1) = Ek_tot(4, n_element_old, isp+1) + W4 * Energy * 0.5

            n_element_bro = HE(6, part_edge_count(2))   !The 0.5 particle distribute CG element

            hx_partition = HP(1, HT(2, n_element_bro)) - HP(1, HT(1, n_element_bro))
            hy_partition = HP(2, HT(4, n_element_bro)) - HP(2, HT(1, n_element_bro))

            dx = (part_x - HP(1, HT(1, n_element_bro))) / hx_partition
            dy = (part_y - HP(2, HT(1, n_element_bro))) / hy_partition
            xcellmdx = 1.0 - dx
            ycellmdy = 1.0 - dy

            IF (delta == 0) THEN
              W1 = xcellmdx * ycellmdy  !local 1
              W2 = dx       * ycellmdy  !local 2
              W3 = dx       * dy        !local 3
              W4 = xcellmdx * dy        !local 4
            ELSEIF (delta == 1) THEN

              IF (ABS(HP(2,HT(1,n_element_bro))-dymin) < SmallValue) THEN
                !There denote the element locate on the bottom edge, it equal to 'j==1' on Program 'hall_zwz'
                BETA = 0.75
              ELSE
                BETA = 1.0
              END IF

              R1 = dymin + HP(2, HT(1, n_element_bro))
              R2 = dymin + HP(2, HT(4, n_element_bro))
              R = part_y
              den = R2*R2 - R1*R1

              W1 = xcellmdx * (R2*R2-R*R) / den * BETA
              W2 = dx       * (R2*R2-R*R) / den * BETA
              W3 = dx       * (R*R-R1*R1) / den
              W4 = xcellmdx * (R*R-R1*R1) / den
            END IF

            rho_p(1, n_element_bro, isp+1) = rho_p(1, n_element_bro, isp+1) + W1 * 0.5
            rho_p(2, n_element_bro, isp+1) = rho_p(2, n_element_bro, isp+1) + W2 * 0.5
            rho_p(3, n_element_bro, isp+1) = rho_p(3, n_element_bro, isp+1) + W3 * 0.5
            rho_p(4, n_element_bro, isp+1) = rho_p(4, n_element_bro, isp+1) + W4 * 0.5
            
            If (ParticleGlobal(isp)%UnequalWeightFlag) Then
              W1 = W1 * ParticleGlobal(isp)%PO(i_part)%WQ
              W2 = W2 * ParticleGlobal(isp)%PO(i_part)%WQ
              W3 = W3 * ParticleGlobal(isp)%PO(i_part)%WQ
              W4 = W4 * ParticleGlobal(isp)%PO(i_part)%WQ
            End If
            
            rho_s(1, n_element_bro, isp+1) = rho_s(1, n_element_bro, isp+1) + W1 * 0.5
            rho_s(2, n_element_bro, isp+1) = rho_s(2, n_element_bro, isp+1) + W2 * 0.5
            rho_s(3, n_element_bro, isp+1) = rho_s(3, n_element_bro, isp+1) + W3 * 0.5
            rho_s(4, n_element_bro, isp+1) = rho_s(4, n_element_bro, isp+1) + W4 * 0.5
            
            Ek_tot(1, n_element_bro, isp+1) = Ek_tot(1, n_element_bro, isp+1) + W1 * Energy * 0.5
            Ek_tot(2, n_element_bro, isp+1) = Ek_tot(2, n_element_bro, isp+1) + W2 * Energy * 0.5
            Ek_tot(3, n_element_bro, isp+1) = Ek_tot(3, n_element_bro, isp+1) + W3 * Energy * 0.5
            Ek_tot(4, n_element_bro, isp+1) = Ek_tot(4, n_element_bro, isp+1) + W4 * Energy * 0.5

            CYCLE   !Jump the big 'Do-End Do' circle, start interpolate the next particle's charge.

          ELSE
            WRITE(6,*) 'Program error, pelase check : part_edge_count(1) == 1!'
            STOP
          END IF
        END IF


      ELSEIF (part_edge_count(1) == 2) THEN
        !The particle locate the edge between CG(0.0) and DG(1.0).
        !The particle locate the edge between DG(0.5) and DG(0.5).

        IF (HE(5, part_edge_count(2)) == HE(5, part_edge_count(3))) THEN
          !The particle locate the edge between CG(0.0) and DG(1.0).(3CG and 1DG, the particle locate the corner node)
          n_element_old = HE(5,part_edge_count(2))

        ELSEIF (HE(5, part_edge_count(2)) /= HE(5, part_edge_count(3))) THEN
          !The particle locate the edge between DG(0.5) and DG(0.5).

          n_element_old = HE(5, part_edge_count(2))   !The 0.5 particle distribute DG element

          hx_partition = HP(1, HT(2, n_element_old)) - HP(1, HT(1, n_element_old))
          hy_partition = HP(2, HT(4, n_element_old)) - HP(2, HT(1, n_element_old))

          dx = (part_x - HP(1, HT(1, n_element_old))) / hx_partition
          dy = (part_y - HP(2, HT(1, n_element_old))) / hy_partition
          xcellmdx = 1.0 - dx
          ycellmdy = 1.0 - dy

          !LY add for axis-symmetric coordinates, 2022-3-22
          IF (delta == 0) THEN
            W1 = xcellmdx * ycellmdy
            W2 = dx       * ycellmdy
            W3 = dx       * dy
            W4 = xcellmdx * dy
          ELSEIF (delta == 1) THEN

            IF (ABS(HP(2,HT(1,n_element_old))-dymin) < SmallValue) THEN
              !There denote the element locate on the bottom edge, it equal to 'j==1' on Program 'hall_zwz'
              BETA = 0.75
            ELSE
              BETA = 1.0
            END IF

            R1 = dymin + HP(2, HT(1, n_element_old))
            R2 = dymin + HP(2, HT(4, n_element_old))
            R = part_y
            den = R2*R2 - R1*R1

            W1 = xcellmdx * (R2*R2-R*R) / den * BETA  !local 1
            W2 = dx       * (R2*R2-R*R) / den * BETA  !local 2
            W3 = dx       * (R*R-R1*R1) / den         !local 3
            W4 = xcellmdx * (R*R-R1*R1) / den         !local 4
          END IF

          rho_p(1, n_element_old, isp+1) = rho_p(1, n_element_old, isp+1) + W1 * 0.5
          rho_p(2, n_element_old, isp+1) = rho_p(2, n_element_old, isp+1) + W2 * 0.5
          rho_p(3, n_element_old, isp+1) = rho_p(3, n_element_old, isp+1) + W3 * 0.5
          rho_p(4, n_element_old, isp+1) = rho_p(4, n_element_old, isp+1) + W4 * 0.5
          
          If (ParticleGlobal(isp)%UnequalWeightFlag) Then
            W1 = W1 * ParticleGlobal(isp)%PO(i_part)%WQ
            W2 = W2 * ParticleGlobal(isp)%PO(i_part)%WQ
            W3 = W3 * ParticleGlobal(isp)%PO(i_part)%WQ
            W4 = W4 * ParticleGlobal(isp)%PO(i_part)%WQ
          End If
            
          rho_s(1, n_element_old, isp+1) = rho_s(1, n_element_old, isp+1) + W1 * 0.5
          rho_s(2, n_element_old, isp+1) = rho_s(2, n_element_old, isp+1) + W2 * 0.5
          rho_s(3, n_element_old, isp+1) = rho_s(3, n_element_old, isp+1) + W3 * 0.5
          rho_s(4, n_element_old, isp+1) = rho_s(4, n_element_old, isp+1) + W4 * 0.5
          
          Ek_tot(1, n_element_old, isp+1) = Ek_tot(1, n_element_old, isp+1) + W1 * Energy * 0.5
          Ek_tot(2, n_element_old, isp+1) = Ek_tot(2, n_element_old, isp+1) + W2 * Energy * 0.5
          Ek_tot(3, n_element_old, isp+1) = Ek_tot(3, n_element_old, isp+1) + W3 * Energy * 0.5
          Ek_tot(4, n_element_old, isp+1) = Ek_tot(4, n_element_old, isp+1) + W4 * Energy * 0.5

          n_element_bro = HE(5, part_edge_count(3))   !The 0.5 particle distribute DG element

          hx_partition = HP(1, HT(2, n_element_bro)) - HP(1, HT(1, n_element_bro))
          hy_partition = HP(2, HT(4, n_element_bro)) - HP(2, HT(1, n_element_bro))

          dx = (part_x - HP(1, HT(1, n_element_bro))) / hx_partition
          dy = (part_y - HP(2, HT(1, n_element_bro))) / hy_partition
          xcellmdx = 1.0 - dx
          ycellmdy = 1.0 - dy

          !LY add for axis-symmetric coordinates, 2022-3-22
          IF (delta == 0) THEN
            W1 = xcellmdx * ycellmdy
            W2 = dx       * ycellmdy
            W3 = dx       * dy
            W4 = xcellmdx * dy
          ELSEIF (delta == 1) THEN

            IF (ABS(HP(2,HT(1,n_element_bro))-dymin) < SmallValue) THEN
              !There denote the element locate on the bottom edge, it equal to 'j==1' on Program 'hall_zwz'
              BETA = 0.75
            ELSE
              BETA = 1.0
            END IF

            R1 = dymin + HP(2, HT(1, n_element_bro))
            R2 = dymin + HP(2, HT(4, n_element_bro))
            R = part_y
            den = R2*R2 - R1*R1

            W1 = xcellmdx * (R2*R2-R*R) / den * BETA  !local 1
            W2 = dx       * (R2*R2-R*R) / den * BETA  !local 2
            W3 = dx       * (R*R-R1*R1) / den         !local 3
            W4 = xcellmdx * (R*R-R1*R1) / den         !local 4
          END IF

          rho_p(1, n_element_bro, isp+1) = rho_p(1, n_element_bro, isp+1) + W1 * 0.5
          rho_p(2, n_element_bro, isp+1) = rho_p(2, n_element_bro, isp+1) + W2 * 0.5
          rho_p(3, n_element_bro, isp+1) = rho_p(3, n_element_bro, isp+1) + W3 * 0.5
          rho_p(4, n_element_bro, isp+1) = rho_p(4, n_element_bro, isp+1) + W4 * 0.5
          
          If (ParticleGlobal(isp)%UnequalWeightFlag) Then
            W1 = W1 * ParticleGlobal(isp)%PO(i_part)%WQ
            W2 = W2 * ParticleGlobal(isp)%PO(i_part)%WQ
            W3 = W3 * ParticleGlobal(isp)%PO(i_part)%WQ
            W4 = W4 * ParticleGlobal(isp)%PO(i_part)%WQ
          Endif
            
          rho_s(1, n_element_bro, isp+1) = rho_s(1, n_element_bro, isp+1) + W1 * 0.5
          rho_s(2, n_element_bro, isp+1) = rho_s(2, n_element_bro, isp+1) + W2 * 0.5
          rho_s(3, n_element_bro, isp+1) = rho_s(3, n_element_bro, isp+1) + W3 * 0.5
          rho_s(4, n_element_bro, isp+1) = rho_s(4, n_element_bro, isp+1) + W4 * 0.5
          
          Ek_tot(1, n_element_bro, isp+1) = Ek_tot(1, n_element_bro, isp+1) + W1 * Energy * 0.5
          Ek_tot(2, n_element_bro, isp+1) = Ek_tot(2, n_element_bro, isp+1) + W2 * Energy * 0.5
          Ek_tot(3, n_element_bro, isp+1) = Ek_tot(3, n_element_bro, isp+1) + W3 * Energy * 0.5
          Ek_tot(4, n_element_bro, isp+1) = Ek_tot(4, n_element_bro, isp+1) + W4 * Energy * 0.5

          CYCLE   !Jump the largest 'Do-End Do' circle, start interpolate the next particle's charge.
        END IF


      ELSEIF (part_edge_count(1) == 4) THEN
        !The particle locate the edge between CG(0.0) and DG(1.0). (condition 1)
        !The particle locate the edge between CG(0.5) and DG(0.5). (condition 2)

        edge_boundary = 0
        DO j = 2, 5
          IF (HE(6, part_edge_count(j)) == 0) THEN
            edge_boundary = edge_boundary + 1
          END IF
        END DO

        IF (edge_boundary == 2) THEN
          n_element_old = HE(5, part_edge_count(2))   !There is to pretect array from outing of bound.

        ELSE

          edge_bro_count = 0
          DO j = 2, 5
            IF (HT(5, HE(5, part_edge_count(j)))>0 .AND. HT(5, HE(6, part_edge_count(j)))==0) THEN
              edge_bro_element(1+edge_bro_count) = HE(6, part_edge_count(j))
              edge_bro_count = edge_bro_count + 1
            END IF
          END DO

          IF (edge_bro_count==2 .AND. edge_bro_element(1)==edge_bro_element(2)) THEN
            !The particle locate the edge between CG(0.5) and DG(0.5). (condition 2: 2CG and 2DG--opposite, CG element is same)

            n_element_old = HE(5, part_edge_count(2))   !The 0.5 particle distribute DG element

            hx_partition = HP(1, HT(2, n_element_old)) - HP(1, HT(1, n_element_old))
            hy_partition = HP(2, HT(4, n_element_old)) - HP(2, HT(1, n_element_old))

            dx = (part_x - HP(1, HT(1, n_element_old))) / hx_partition
            dy = (part_y - HP(2, HT(1, n_element_old))) / hy_partition
            xcellmdx = 1.0 - dx
            ycellmdy = 1.0 - dy

            !LY add for axis-symmetric coordinates, 2022-3-22
            IF (delta == 0) THEN
              W1 = xcellmdx * ycellmdy
              W2 = dx       * ycellmdy
              W3 = dx       * dy
              W4 = xcellmdx * dy
            ELSEIF (delta == 1) THEN

              IF (ABS(HP(2,HT(1,n_element_old))-dymin) < SmallValue) THEN
                !There denote the element locate on the bottom edge, it equal to 'j==1' on Program 'hall_zwz'
                BETA = 0.75
              ELSE
                BETA = 1.0
              END IF

              R1 = dymin + HP(2, HT(1, n_element_old))
              R2 = dymin + HP(2, HT(4, n_element_old))
              R = part_y
              den = R2*R2 - R1*R1

              W1 = xcellmdx * (R2*R2-R*R) / den * BETA  !local 1
              W2 = dx       * (R2*R2-R*R) / den * BETA  !local 2
              W3 = dx       * (R*R-R1*R1) / den         !local 3
              W4 = xcellmdx * (R*R-R1*R1) / den         !local 4
            END IF

            rho_p(1, n_element_old, isp+1) = rho_p(1, n_element_old, isp+1) + W1 * 0.5
            rho_p(2, n_element_old, isp+1) = rho_p(2, n_element_old, isp+1) + W2 * 0.5
            rho_p(3, n_element_old, isp+1) = rho_p(3, n_element_old, isp+1) + W3 * 0.5
            rho_p(4, n_element_old, isp+1) = rho_p(4, n_element_old, isp+1) + W4 * 0.5
            
            If (ParticleGlobal(isp)%UnequalWeightFlag) Then
              W1 = W1 * ParticleGlobal(isp)%PO(i_part)%WQ
              W2 = W2 * ParticleGlobal(isp)%PO(i_part)%WQ
              W3 = W3 * ParticleGlobal(isp)%PO(i_part)%WQ
              W4 = W4 * ParticleGlobal(isp)%PO(i_part)%WQ
            End If
            
            rho_s(1, n_element_old, isp+1) = rho_s(1, n_element_old, isp+1) + W1 * 0.5
            rho_s(2, n_element_old, isp+1) = rho_s(2, n_element_old, isp+1) + W2 * 0.5
            rho_s(3, n_element_old, isp+1) = rho_s(3, n_element_old, isp+1) + W3 * 0.5
            rho_s(4, n_element_old, isp+1) = rho_s(4, n_element_old, isp+1) + W4 * 0.5
            
            Ek_tot(1, n_element_old, isp+1) = Ek_tot(1, n_element_old, isp+1) + W1 * Energy * 0.5
            Ek_tot(2, n_element_old, isp+1) = Ek_tot(2, n_element_old, isp+1) + W2 * Energy * 0.5
            Ek_tot(3, n_element_old, isp+1) = Ek_tot(3, n_element_old, isp+1) + W3 * Energy * 0.5
            Ek_tot(4, n_element_old, isp+1) = Ek_tot(4, n_element_old, isp+1) + W4 * Energy * 0.5

            n_element_bro = edge_bro_element(1)   !The 0.5 particle distribute CG element

            hx_partition = HP(1, HT(2, n_element_bro)) - HP(1, HT(1, n_element_bro))
            hy_partition = HP(2, HT(4, n_element_bro)) - HP(2, HT(1, n_element_bro))

            dx = (part_x - HP(1, HT(1, n_element_bro))) / hx_partition
            dy = (part_y - HP(2, HT(1, n_element_bro))) / hy_partition
            xcellmdx = 1.0 - dx
            ycellmdy = 1.0 - dy

            !LY add for axis-symmetric coordinates, 2022-3-22
            IF (delta == 0) THEN
              W1 = xcellmdx * ycellmdy
              W2 = dx       * ycellmdy
              W3 = dx       * dy
              W4 = xcellmdx * dy
            ELSEIF (delta == 1) THEN

              IF (ABS(HP(2,HT(1,n_element_bro))-dymin) < SmallValue) THEN
                !There denote the element locate on the bottom edge, it equal to 'j==1' on Program 'hall_zwz'
                BETA = 0.75
              ELSE
                BETA = 1.0
              END IF

              R1 = dymin + HP(2, HT(1, n_element_bro))
              R2 = dymin + HP(2, HT(4, n_element_bro))
              R = part_y
              den = R2*R2 - R1*R1

              W1 = xcellmdx * (R2*R2-R*R) / den * BETA  !local 1
              W2 = dx       * (R2*R2-R*R) / den * BETA  !local 2
              W3 = dx       * (R*R-R1*R1) / den         !local 3
              W4 = xcellmdx * (R*R-R1*R1) / den         !local 4
            END IF

            rho_p(1, n_element_bro, isp+1) = rho_p(1, n_element_bro, isp+1) + W1 * 0.5
            rho_p(2, n_element_bro, isp+1) = rho_p(2, n_element_bro, isp+1) + W2 * 0.5
            rho_p(3, n_element_bro, isp+1) = rho_p(3, n_element_bro, isp+1) + W3 * 0.5
            rho_p(4, n_element_bro, isp+1) = rho_p(4, n_element_bro, isp+1) + W4 * 0.5
            
            If (ParticleGlobal(isp)%UnequalWeightFlag) Then
              W1 = W1 * ParticleGlobal(isp)%PO(i_part)%WQ
              W2 = W2 * ParticleGlobal(isp)%PO(i_part)%WQ
              W3 = W3 * ParticleGlobal(isp)%PO(i_part)%WQ
              W4 = W4 * ParticleGlobal(isp)%PO(i_part)%WQ
            End If
            
            rho_s(1, n_element_bro, isp+1) = rho_s(1, n_element_bro, isp+1) + W1 * 0.5
            rho_s(2, n_element_bro, isp+1) = rho_s(2, n_element_bro, isp+1) + W2 * 0.5
            rho_s(3, n_element_bro, isp+1) = rho_s(3, n_element_bro, isp+1) + W3 * 0.5
            rho_s(4, n_element_bro, isp+1) = rho_s(4, n_element_bro, isp+1) + W4 * 0.5
            
            Ek_tot(1, n_element_bro, isp+1) = Ek_tot(1, n_element_bro, isp+1) + W1 * Energy * 0.5
            Ek_tot(2, n_element_bro, isp+1) = Ek_tot(2, n_element_bro, isp+1) + W2 * Energy * 0.5
            Ek_tot(3, n_element_bro, isp+1) = Ek_tot(3, n_element_bro, isp+1) + W3 * Energy * 0.5
            Ek_tot(4, n_element_bro, isp+1) = Ek_tot(4, n_element_bro, isp+1) + W4 * Energy * 0.5

            CYCLE   !Jump the largest 'Do-End Do' circle, start interpolate the next particle's charge.

          ELSEIF (edge_bro_count==2 .AND. edge_bro_element(1)/=edge_bro_element(2)) THEN
            !The particle locate the edge between CG(0.0) and DG(1.0). (condition 1: 2CG and 2DG--opposite, CG element is different)
            n_element_old = HE(5, part_edge_count(2))

          ELSEIF (edge_bro_count==4) THEN
            !The particle locate the edge between CG(0.0) and DG(1.0). (condition 1: 2CG and 2DG--diagonal, CG element is different)
            n_element_old = HE(5, part_edge_count(2))
          END IF
        END IF


      ELSEIF (part_edge_count(1) == 6) THEN
        !The particle locate the edge between CG(1.0) and DG(0.0).(1CG and 3DG, the particle locate the corner node)��
        !The particle locate the edge between DG(0.5) and DG(0.5).(2DG and 1DG, the particle locate the hand node)��

        edge_bro_count = 0
        edge_temp_count = 0
        edge_self_element = 0
        edge_self_element_min = MINVAL(HT(5, HE(5,part_edge_count(2:7))))

        DO j = 2, 7
          IF (HT(5, HE(5,part_edge_count(j)))>0 .AND. HT(5, HE(6,part_edge_count(j)))==0) THEN
            edge_bro_element(1+edge_bro_count) = HE(6, part_edge_count(j))
            edge_bro_count = edge_bro_count + 1
          ELSEIF (HT(5, HE(5,part_edge_count(j)))>0 .AND. HT(5, HE(6,part_edge_count(j)))>0) THEN
            edge_temp_count = edge_temp_count + 1
            IF (edge_self_element_min == HT(5, HE(5,part_edge_count(j)))) THEN
              edge_self_element = part_edge_count(j)
            END IF
          END IF
        END DO

        IF (edge_bro_count==2 .AND. edge_bro_element(1)==edge_bro_element(2)) THEN
          !The particle locate the edge between CG(1.0) and DG(0.0).(1CG and 3DG, the particle locate the corner node)
          n_element_old = edge_bro_element(1)

        ELSEIF (edge_temp_count == 6) THEN
          !The particle locate the edge between DG(0.5) and DG(0.5).(2DG and 1DG, the particle locate the hand node)

          n_element_old = HE(5, edge_self_element)

          hx_partition = HP(1, HT(2, n_element_old)) - HP(1, HT(1, n_element_old))
          hy_partition = HP(2, HT(4, n_element_old)) - HP(2, HT(1, n_element_old))

          dx = (part_x - HP(1, HT(1, n_element_old))) / hx_partition
          dy = (part_y - HP(2, HT(1, n_element_old))) / hy_partition
          xcellmdx = 1.0 - dx
          ycellmdy = 1.0 - dy

          !LY add for axis-symmetric coordinates, 2022-3-22
          IF (delta == 0) THEN
            W1 = xcellmdx * ycellmdy
            W2 = dx       * ycellmdy
            W3 = dx       * dy
            W4 = xcellmdx * dy
          ELSEIF (delta == 1) THEN

            IF (ABS(HP(2,HT(1,n_element_old))-dymin) < SmallValue) THEN
              !There denote the element locate on the bottom edge, it equal to 'j==1' on Program 'hall_zwz'
              BETA = 0.75
            ELSE
              BETA = 1.0
            END IF

            R1 = dymin + HP(2, HT(1, n_element_old))
            R2 = dymin + HP(2, HT(4, n_element_old))
            R = part_y
            den = R2*R2 - R1*R1

            W1 = xcellmdx * (R2*R2-R*R) / den * BETA  !local 1
            W2 = dx       * (R2*R2-R*R) / den * BETA  !local 2
            W3 = dx       * (R*R-R1*R1) / den         !local 3
            W4 = xcellmdx * (R*R-R1*R1) / den         !local 4
          END IF

          rho_p(1, n_element_old, isp+1) = rho_p(1, n_element_old, isp+1) + W1 * 0.5
          rho_p(2, n_element_old, isp+1) = rho_p(2, n_element_old, isp+1) + W2 * 0.5
          rho_p(3, n_element_old, isp+1) = rho_p(3, n_element_old, isp+1) + W3 * 0.5
          rho_p(4, n_element_old, isp+1) = rho_p(4, n_element_old, isp+1) + W4 * 0.5
            
          If (ParticleGlobal(isp)%UnequalWeightFlag) Then
            W1 = W1 * ParticleGlobal(isp)%PO(i_part)%WQ
            W2 = W2 * ParticleGlobal(isp)%PO(i_part)%WQ
            W3 = W3 * ParticleGlobal(isp)%PO(i_part)%WQ
            W4 = W4 * ParticleGlobal(isp)%PO(i_part)%WQ
          End If
          
          rho_s(1, n_element_old, isp+1) = rho_s(1, n_element_old, isp+1) + W1 * 0.5
          rho_s(2, n_element_old, isp+1) = rho_s(2, n_element_old, isp+1) + W2 * 0.5
          rho_s(3, n_element_old, isp+1) = rho_s(3, n_element_old, isp+1) + W3 * 0.5
          rho_s(4, n_element_old, isp+1) = rho_s(4, n_element_old, isp+1) + W4 * 0.5
          
          Ek_tot(1, n_element_old, isp+1) = Ek_tot(1, n_element_old, isp+1) + W1 * Energy * 0.5
          Ek_tot(2, n_element_old, isp+1) = Ek_tot(2, n_element_old, isp+1) + W2 * Energy * 0.5
          Ek_tot(3, n_element_old, isp+1) = Ek_tot(3, n_element_old, isp+1) + W3 * Energy * 0.5
          Ek_tot(4, n_element_old, isp+1) = Ek_tot(4, n_element_old, isp+1) + W4 * Energy * 0.5

          IF (HE(6, edge_self_element)==0) THEN     !protect array out boundary.
            WRITE(6,*) 'ERROR DG BROTHER ELEMENT2'
            STOP
          END IF

          n_element_bro = HE(6, edge_self_element)   !The 0.5 particle distribute DG element

          hx_partition = HP(1, HT(2, n_element_bro)) - HP(1, HT(1, n_element_bro))
          hy_partition = HP(2, HT(4, n_element_bro)) - HP(2, HT(1, n_element_bro))

          dx = (part_x - HP(1, HT(1, n_element_bro))) / hx_partition
          dy = (part_y - HP(2, HT(1, n_element_bro))) / hy_partition
          xcellmdx = 1.0 - dx
          ycellmdy = 1.0 - dy

          !LY add for axis-symmetric coordinates, 2022-3-22
          IF (delta == 0) THEN
            W1 = xcellmdx * ycellmdy
            W2 = dx       * ycellmdy
            W3 = dx       * dy
            W4 = xcellmdx * dy
          ELSEIF (delta == 1) THEN

            IF (ABS(HP(2,HT(1,n_element_bro))-dymin) < SmallValue) THEN
              !There denote the element locate on the bottom edge, it equal to 'j==1' on Program 'hall_zwz'
              BETA = 0.75
            ELSE
              BETA = 1.0
            END IF

            R1 = dymin + HP(2, HT(1, n_element_bro))
            R2 = dymin + HP(2, HT(4, n_element_bro))
            R = part_y
            den = R2*R2 - R1*R1

            W1 = xcellmdx * (R2*R2-R*R) / den * BETA  !local 1
            W2 = dx       * (R2*R2-R*R) / den * BETA  !local 2
            W3 = dx       * (R*R-R1*R1) / den         !local 3
            W4 = xcellmdx * (R*R-R1*R1) / den         !local 4
          END IF

          rho_p(1, n_element_bro, isp+1) = rho_p(1, n_element_bro, isp+1) + W1 * 0.5
          rho_p(2, n_element_bro, isp+1) = rho_p(2, n_element_bro, isp+1) + W2 * 0.5
          rho_p(3, n_element_bro, isp+1) = rho_p(3, n_element_bro, isp+1) + W3 * 0.5
          rho_p(4, n_element_bro, isp+1) = rho_p(4, n_element_bro, isp+1) + W4 * 0.5
          
          If (ParticleGlobal(isp)%UnequalWeightFlag) Then
            W1 = W1 * ParticleGlobal(isp)%PO(i_part)%WQ
            W2 = W2 * ParticleGlobal(isp)%PO(i_part)%WQ
            W3 = W3 * ParticleGlobal(isp)%PO(i_part)%WQ
            W4 = W4 * ParticleGlobal(isp)%PO(i_part)%WQ
          End If
          
          rho_s(1, n_element_bro, isp+1) = rho_s(1, n_element_bro, isp+1) + W1 * 0.5
          rho_s(2, n_element_bro, isp+1) = rho_s(2, n_element_bro, isp+1) + W2 * 0.5
          rho_s(3, n_element_bro, isp+1) = rho_s(3, n_element_bro, isp+1) + W3 * 0.5
          rho_s(4, n_element_bro, isp+1) = rho_s(4, n_element_bro, isp+1) + W4 * 0.5
          
          Ek_tot(1, n_element_bro, isp+1) = Ek_tot(1, n_element_bro, isp+1) + W1 * Energy * 0.5
          Ek_tot(2, n_element_bro, isp+1) = Ek_tot(2, n_element_bro, isp+1) + W2 * Energy * 0.5
          Ek_tot(3, n_element_bro, isp+1) = Ek_tot(3, n_element_bro, isp+1) + W3 * Energy * 0.5
          Ek_tot(4, n_element_bro, isp+1) = Ek_tot(4, n_element_bro, isp+1) + W4 * Energy * 0.5

          CYCLE   !Jump the largest 'Do-End Do' circle, start interpolate the next particle's charge.
        END IF


      ELSEIF (part_edge_count(1) == 8) THEN
        !The particle locate the edge between DG(1.0) and DG(0.0).(4DG, the particle locate the corner node)��
        n_element_old = HE(5, part_edge_count(2))

      ELSEIF (part_edge_count(1) == 0) THEN
        !The particle locate on the edge of CG element without DG neighbor.
        n_element_old = CellMesh(n_element_old)%Finalindex

      ELSE
        WRITE(6,*) 'part_edge_count(1) = ', part_edge_count(1)
        WRITE(6,*) 'Part_index = ', i_part
        STOP
      END IF
      !========================particle locate on CG and DG intersect edge: CG, DG coupling===========================
    END IF
    
    hx_partition = HP(1, HT(2, n_element_old)) - HP(1, HT(1, n_element_old))
    hy_partition = HP(2, HT(4, n_element_old)) - HP(2, HT(1, n_element_old))

    dx = (part_x - HP(1, HT(1, n_element_old))) / hx_partition
    dy = (part_y - HP(2, HT(1, n_element_old))) / hy_partition

    xcellmdx = 1.0 - dx
    ycellmdy = 1.0 - dy

    !LY add for axis-symmetric coordinates, 2022-3-22
    IF (delta == 0) THEN
      W1 = xcellmdx * ycellmdy
      W2 = dx       * ycellmdy
      W3 = dx       * dy
      W4 = xcellmdx * dy
    ELSEIF (delta == 1) THEN

      IF (ABS(HP(2,HT(1,n_element_old))-dymin) < SmallValue) THEN
        !There denote the element locate on the bottom edge, it equal to 'j==1' on Program 'hall_zwz'
        BETA = 0.75
      ELSE
        BETA = 1.0
      END IF

      R1 = dymin + HP(2, HT(1, n_element_old))
      R2 = dymin + HP(2, HT(4, n_element_old))
      R = part_y
      den = R2*R2 - R1*R1

      W1 = xcellmdx * (R2*R2-R*R) / den * BETA  !local 1
      W2 = dx       * (R2*R2-R*R) / den * BETA  !local 2
      W3 = dx       * (R*R-R1*R1) / den         !local 3
      W4 = xcellmdx * (R*R-R1*R1) / den         !local 4
    END IF

    rho_p(1, n_element_old, isp+1) = rho_p(1, n_element_old, isp+1) + W1
    rho_p(2, n_element_old, isp+1) = rho_p(2, n_element_old, isp+1) + W2
    rho_p(3, n_element_old, isp+1) = rho_p(3, n_element_old, isp+1) + W3
    rho_p(4, n_element_old, isp+1) = rho_p(4, n_element_old, isp+1) + W4
    
    If (ParticleGlobal(isp)%UnequalWeightFlag) Then
      W1 = W1 * ParticleGlobal(isp)%PO(i_part)%WQ
      W2 = W2 * ParticleGlobal(isp)%PO(i_part)%WQ
      W3 = W3 * ParticleGlobal(isp)%PO(i_part)%WQ
      W4 = W4 * ParticleGlobal(isp)%PO(i_part)%WQ
    End If
            
    rho_s(1, n_element_old, isp+1) = rho_s(1, n_element_old, isp+1) + W1
    rho_s(2, n_element_old, isp+1) = rho_s(2, n_element_old, isp+1) + W2
    rho_s(3, n_element_old, isp+1) = rho_s(3, n_element_old, isp+1) + W3
    rho_s(4, n_element_old, isp+1) = rho_s(4, n_element_old, isp+1) + W4
    
    Ek_tot(1, n_element_old, isp+1) = Ek_tot(1, n_element_old, isp+1) + W1 * Energy
    Ek_tot(2, n_element_old, isp+1) = Ek_tot(2, n_element_old, isp+1) + W2 * Energy
    Ek_tot(3, n_element_old, isp+1) = Ek_tot(3, n_element_old, isp+1) + W3 * Energy
    Ek_tot(4, n_element_old, isp+1) = Ek_tot(4, n_element_old, isp+1) + W4 * Energy
    
  End Do
End Do

!==========Part 2 : Sum mesh nodes charge weight from basis nodes==========
Do isp = 0, ControlFlowGlobal%Ns
  Do i = 1, Size(HT,2)
    rho_p_sum(1, HT(1, i), isp+1) = rho_p_sum(1, HT(1, i), isp+1) + rho_p(1, i, isp+1)
    rho_p_sum(1, HT(2, i), isp+1) = rho_p_sum(1, HT(2, i), isp+1) + rho_p(2, i, isp+1)
    rho_p_sum(1, HT(3, i), isp+1) = rho_p_sum(1, HT(3, i), isp+1) + rho_p(3, i, isp+1)
    rho_p_sum(1, HT(4, i), isp+1) = rho_p_sum(1, HT(4, i), isp+1) + rho_p(4, i, isp+1)
    
    rho_s_sum(1, HT(1, i), isp+1) = rho_s_sum(1, HT(1, i), isp+1) + rho_s(1, i, isp+1)
    rho_s_sum(1, HT(2, i), isp+1) = rho_s_sum(1, HT(2, i), isp+1) + rho_s(2, i, isp+1)
    rho_s_sum(1, HT(3, i), isp+1) = rho_s_sum(1, HT(3, i), isp+1) + rho_s(3, i, isp+1)
    rho_s_sum(1, HT(4, i), isp+1) = rho_s_sum(1, HT(4, i), isp+1) + rho_s(4, i, isp+1)
    
    Ek_tot_sum(1, HT(1, i), isp+1) = Ek_tot_sum(1, HT(1, i), isp+1) + Ek_tot(1, i, isp+1)
    Ek_tot_sum(1, HT(2, i), isp+1) = Ek_tot_sum(1, HT(2, i), isp+1) + Ek_tot(2, i, isp+1) 
    Ek_tot_sum(1, HT(3, i), isp+1) = Ek_tot_sum(1, HT(3, i), isp+1) + Ek_tot(3, i, isp+1)
    Ek_tot_sum(1, HT(4, i), isp+1) = Ek_tot_sum(1, HT(4, i), isp+1) + Ek_tot(4, i, isp+1)
  End Do
End Do

!==========Part 3 : Calculate the charge density of mesh nodes==========
Do isp = 1, ispe_tot
  Do i = 1, Size(HP,2)
    
    !num = 0
    !Do j = 3, 6
    !  If (INT(P_average(j,i)) /= 0) Then
    !    num = num + 1
    !  End If
    !End Do
    
    !Do k = 1, num
      rho_p_mesh_sum(i, 1, isp) = rho_p_mesh_sum(i, 1, isp) + rho_p_sum(1, i, isp)
      rho_sum(i, 1, isp) = rho_sum(i, 1, isp) + rho_s_sum(1, i, isp)
      Ek_tot_mesh_sum(i, 1, isp) = Ek_tot_mesh_sum(i, 1, isp) + Ek_tot_sum(1, i, isp)
    !End Do
    
    If (Ek_tot_mesh_sum(i,1,isp)/=0 .AND. rho_sum(i,1,isp)/=0) Then
      Ek_s_mesh_sum(i, 1, isp) = Ek_tot_mesh_sum(i, 1, isp) / rho_sum(i, 1, isp)
    Elseif (Ek_tot_mesh_sum(i,1,isp)/=0 .AND. rho_sum(i,1,isp)==0) Then
      Print*, 'something wrong'
      Print*, 'error locate on GetParChg_2D.f90--Ek_tot_mesh_sum(i,1,isp)/=0 .AND. rho_sum(i,1,isp)==0'
      Stop
    End If
    
    rho_n_mesh_sum(i, 1, isp) = rho_sum(i, 1, isp)
    
    If (delta == 0) Then  !2D Cartesian coordinates.
      rho_sum(i, 1, isp) = ParticleGlobal(isp-1)%Weight * qs(isp) * rho_sum(i, 1, isp) / L_ref / L_ref / n_ref / Cell_Volume_zwz(2, i)
      Ek_tot_mesh_sum(i, 1, isp) = ParticleGlobal(isp-1)%Weight * Ek_tot_mesh_sum(i, 1, isp) / L_ref / L_ref / n_ref / Cell_Volume_zwz(2, i)
    Elseif (delta == 1) Then  !2D axis-symmetric cooridnates.
      If (ParticleGlobal(isp-1)%UnequalWeightFlag) Then
        rho_sum(i, 1, isp) = qs(isp) * rho_sum(i, 1, isp) / L_ref / L_ref / L_ref / n_ref / Cell_Volume_zwz(2, i)
        Ek_tot_mesh_sum(i, 1, isp) = Ek_tot_mesh_sum(i, 1, isp) / L_ref / L_ref / L_ref / n_ref / Cell_Volume_zwz(2, i)
      Else
        rho_sum(i, 1, isp) = ParticleGlobal(isp-1)%Weight * qs(isp) * rho_sum(i, 1, isp) / L_ref / L_ref / L_ref / n_ref / Cell_Volume_zwz(2, i)
        Ek_tot_mesh_sum(i, 1, isp) = ParticleGlobal(isp-1)%Weight * Ek_tot_mesh_sum(i, 1, isp) / L_ref / L_ref / L_ref / n_ref / Cell_Volume_zwz(2, i)
      End If
    End If
    
    rho(i, 1) = rho(i, 1) + rho_sum(i, 1, isp)
    
  End Do
End Do

!==========Part 4 : Result Visualization==========
!Do i = 1, Size(P_average,2)
!  
!  num = 0
!  Do j = 3, 6
!    If (INT(P_average(j,i)) /= 0) Then
!      num = num +1
!    End If
!  End Do
!  
!  rho_test(INT(P_average(3:2+num, i)), 1) = rho(i, 1)
!  
!  Do temp_count = 1, ispe_tot
!    rho_p_test(INT(P_average(3:2+num, i)), 1, temp_count) = rho_p_mesh_sum(i, 1, temp_count)
!    rho_s_test(INT(P_average(3:2+num, i)), 1, temp_count) = rho_sum(i, 1, temp_count)
!    Ek_tot_test(INT(P_average(3:2+num, i)), 1, temp_count) = Ek_tot_mesh_sum(i, 1, temp_count)
!    Ek_s(INT(P_average(3:2+num, i)), 1, temp_count) = Ek_s_mesh_sum(i, 1, temp_count)
!    rho_n(INT(P_average(3:2+num, i)), 1, temp_count) = rho_n_mesh_sum(i, 1, temp_count)
!  End Do
!  
!End Do

!We need again define the store space of rho, rho_s, rho_p
Deallocate(rho_s, rho_p)
Deallocate(Ek_tot)

Allocate(Ek_tot(Size(HP,2), 1, ispe_tot))
Ek_tot(1:Size(HP,2), 1, 1:ispe_tot) = 0.0

Allocate(rho_s(Size(HP,2),1,ispe_tot), rho_p(Size(HP,2),1,ispe_tot))
!rho(1:Size(HP,2), 1) = 0.0
rho_s(1:Size(HP,2), 1, 1:ispe_tot) = 0.0
rho_p(1:Size(HP,2), 1, 1:ispe_tot) = 0.0

Do i = 1, Size(HP,2)
  !rho(i, 1) = rho_test(i, 1)
  rho_s(i, 1, 1:ispe_tot) = rho_sum(i, 1, 1:ispe_tot)
  rho_p(i, 1, 1:ispe_tot) = rho_p_mesh_sum(i, 1, 1:ispe_tot)
  Ek_tot(i, 1, 1:ispe_tot) = Ek_tot_mesh_sum(i, 1, 1:ispe_tot)
End Do
!=========LY modification for Multi-Layer-Grid, 2022-7-25=========



!!=========ZWZ brother Old Code=========
!! Zero out the density (including the guard cells)
!rho_s	= Zero
!rho_p = Zero
!rho		= Zero
!Ek_tot = Zero
!
!n_node = 0
!
!!$ mb.ZWZ for using JW's particle data structure =======================\\
!Do isp=0,ControlFlowGlobal%Ns
!    Do i_part = 1, ParticleGlobal(isp)%NPar
!        rxp = (ParticleGlobal(isp)%PO(i_part)%X - Vert_o(1))*hxi(1)
!	    i = rxp
!	    dx = rxp-i
!
!	    ryp = (ParticleGlobal(isp)%PO(i_part)%Y - Vert_o(2))*hxi(2)
!	    j = ryp
!	    dy = ryp-j
!	
!	    xcellmdx = 1. -dx
!	    ycellmdy = 1. -dy
!        IF( delta == 0 ) THEN
!            P1 = xcellmdx*ycellmdy
!            P2 = dx*ycellmdy
!            P3 = xcellmdx*dy
!            P4 = dx*dy
!        ELSEIF( delta == 1 ) THEN
!            IF (j == 1) THEN
!                !BETA = 0.5
!                BETA = 0.75
!            ELSE 
!                BETA = 1.
!            ENDIF
!        
!            R1 = dymin + float(j - 1)*hx(2)
!            R2 = dymin + float(j)*hx(2)
!            R = ParticleGlobal(isp)%PO(i_part)%Y 
!            den = R2*R2-R1*R1
!        
!            P1 = xcellmdx*(R2*R2-R*R)/den * BETA 
!            P2 = dx      *(R2*R2-R*R)/den * BETA 
!            P3 = xcellmdx*(R*R-R1*R1)/den 
!            P4 = dx      *(R*R-R1*R1)/den 
!        ENDIF
!        
!        rho_p(i  ,j  ,isp+1) = rho_p(i  ,j  ,isp+1) + P1
!        rho_p(i+1,j  ,isp+1) = rho_p(i+1,j  ,isp+1) + P2
!        rho_p(i  ,j+1,isp+1) = rho_p(i  ,j+1,isp+1) + P3
!        rho_p(i+1,j+1,isp+1) = rho_p(i+1,j+1,isp+1) + P4
!        
!        If (ParticleGlobal(isp)%UnequalWeightFlag) Then
!            P1 = P1  * ParticleGlobal(isp)%PO(i_part)%WQ
!            P2 = P2  * ParticleGlobal(isp)%PO(i_part)%WQ
!            P3 = P3  * ParticleGlobal(isp)%PO(i_part)%WQ
!            P4 = P4  * ParticleGlobal(isp)%PO(i_part)%WQ
!        Endif
!        
!      rho_s(i  ,j  ,isp+1) = rho_s(i  ,j  ,isp+1) + P1
!	    rho_s(i+1,j  ,isp+1) = rho_s(i+1,j  ,isp+1) + P2 
!	    rho_s(i  ,j+1,isp+1) = rho_s(i  ,j+1,isp+1) + P3 
!	    rho_s(i+1,j+1,isp+1) = rho_s(i+1,j+1,isp+1) + P4 
!    
!     !   V2 = ParticleGlobal(isp)%PO(i_part)%Vx**2 + ParticleGlobal(isp)%PO(i_part)%Vy**2 + ParticleGlobal(isp)%PO(i_part)%Vz**2
!     !   Ek_tot(i  ,j  ,isp+1) = Ek_tot(i  ,j  ,isp+1) + P1 * V2 
!	    !Ek_tot(i+1,j  ,isp+1) = Ek_tot(i+1,j  ,isp+1) + P2 * V2 
!	    !Ek_tot(i  ,j+1,isp+1) = Ek_tot(i  ,j+1,isp+1) + P3 * V2 
!	    !Ek_tot(i+1,j+1,isp+1) = Ek_tot(i+1,j+1,isp+1) + P4 * V2 
!        
!      Energy = ParticleGlobal(isp)%PO(i_part)%Energy(ParticleGlobal(isp)%Mass,ParticleGlobal(isp)%VFactor)/JtoeV
!      Ek_tot(i  ,j  ,isp+1) = Ek_tot(i  ,j  ,isp+1) + P1 * Energy
!	    Ek_tot(i+1,j  ,isp+1) = Ek_tot(i+1,j  ,isp+1) + P2 * Energy
!	    Ek_tot(i  ,j+1,isp+1) = Ek_tot(i  ,j+1,isp+1) + P3 * Energy
!	    Ek_tot(i+1,j+1,isp+1) = Ek_tot(i+1,j+1,isp+1) + P4 * Energy
!    End Do
!End Do
!!$ mb.ZWZ for using JW's particle data structure =====================//
!
!! Particle loop
!!DO i_part = 1, ntot
!!    
!!	isp = INT(part(i_part,7))
!!	
!!	!IF (part(i_part,1)==f_left_wall(1)) part(i_part,1)=part(i_part,1)+0.0001
!!	!IF (part(i_part,1)==f_right_wall(1)) part(i_part,1)=part(i_part,1)-0.0001
!!	!IF (part(i_part,2)==f_left_wall(2)) part(i_part,2)=part(i_part,2)+0.0001
!!	!IF (part(i_part,2)==f_right_wall(2)) part(i_part,2)=part(i_part,2)-0.0001
!!
!!	rxp = (part(i_part,1) - Vert_o(1))*hxi(1)
!!	i = rxp
!!	dx = rxp-i
!!
!!	ryp = (part(i_part,2) - Vert_o(2))*hxi(2)
!!	j = ryp
!!	dy = ryp-j
!!	
!!	n = (j+(i-1)*ny)
!!    n_element = (j+(i-1)*(ny-1))
!!    
!!    x_part = part(i_part,1)
!!    y_part = part(i_part,2)
!!
!!	xcellmdx = 1. -dx
!!	ycellmdy = 1. -dy
!!    
!!    !$ =============== mb.ZWZ 2021/7/11 =================== \\
!!    IF( delta == 0 ) THEN
!!        P1 = xcellmdx*ycellmdy
!!        P2 = dx*ycellmdy
!!        P3 = xcellmdx*dy
!!        P4 = dx*dy
!!    ELSEIF( delta == 1 ) THEN
!!        !YFACTOR = part(i_part,2)/YINI(i_part)
!!        !YFACTOR = part(i_part,2)/R_ref
!!        YFACTOR = 1.
!!        !Rweight(i_part) = 1.
!!        IF (j == 1) THEN
!!            BETA = 0.5
!!        ELSE 
!!            BETA = 1.
!!        ENDIF
!!        
!!        R1 = f_left_wall(2) + float(j - 1)*hx(2)
!!        R2 = f_left_wall(2) + float(j)*hx(2)
!!        R = part(i_part,2)
!!        den=R2*R2-R1*R1
!!        !
!!        !P1 = YFACTOR * xcellmdx*(R2*R2-R*R)/den * BETA * Rweight(i_part)
!!        !P2 = YFACTOR * dx      *(R2*R2-R*R)/den * BETA * Rweight(i_part)
!!        !P3 = YFACTOR * xcellmdx*(R*R-R1*R1)/den * BETA * Rweight(i_part)
!!        !P4 = YFACTOR * dx      *(R*R-R1*R1)/den * BETA * Rweight(i_part)
!!        
!!        P1 = YFACTOR * xcellmdx*(R2*R2-R*R)/den * BETA 
!!        P2 = YFACTOR * dx      *(R2*R2-R*R)/den * BETA 
!!        P3 = YFACTOR * xcellmdx*(R*R-R1*R1)/den 
!!        P4 = YFACTOR * dx      *(R*R-R1*R1)/den 
!!        
!!        !if (j>0) THEN
!!            !print*,'j=',j
!!            !print*,P1,P2,P3,P4
!!            !pause
!!        !ENDIF
!!    ENDIF
!!    
!!    
!!    rho_s(i  ,j  ,isp) = rho_s(i  ,j  ,isp) + P1
!!	rho_s(i+1,j  ,isp) = rho_s(i+1,j  ,isp) + P2
!!	rho_s(i  ,j+1,isp) = rho_s(i  ,j+1,isp) + P3
!!	rho_s(i+1,j+1,isp) = rho_s(i+1,j+1,isp) + P4
!!    
!!    !$ =============== mb.ZWZ 2021/7/11 ================= //
!!    
!!    !$ =============== mb.ZWZ 2021/7/11 ================= \\
!!    V2 = part(i_part,4)**2 + part(i_part,5)**2 + part(i_part,6)**2
!!    Ek_tot(i  ,j  ,isp) = Ek_tot(i  ,j  ,isp) + P1 * V2
!!	Ek_tot(i+1,j  ,isp) = Ek_tot(i+1,j  ,isp) + P2 * V2
!!	Ek_tot(i  ,j+1,isp) = Ek_tot(i  ,j+1,isp) + P3 * V2
!!	Ek_tot(i+1,j+1,isp) = Ek_tot(i+1,j+1,isp) + P4 * V2
!!    !$ =============== mb.ZWZ 2021/7/11 ================= //
!!    
!! !   rho_s(i  ,j  ,isp) = rho_s(i  ,j  ,isp) + xcellmdx*ycellmdy
!!	!rho_s(i+1,j  ,isp) = rho_s(i+1,j  ,isp) + dx*ycellmdy
!!	!rho_s(i  ,j+1,isp) = rho_s(i  ,j+1,isp) + xcellmdx*dy
!!	!rho_s(i+1,j+1,isp) = rho_s(i+1,j+1,isp) + dx*dy
!!END DO 
!
!!$ == descripted by ZWZ for using different cell volume ==\\
!
!!! 处理边界
!!rho_s(1,1,:) = 4 * rho_s(1,1,:)
!!rho_s(nx,1,:) = 4 * rho_s(nx,1,:)
!!rho_s(1,ny,:) = 4 * rho_s(1,ny,:)
!!rho_s(nx,ny,:) = 4 * rho_s(nx,ny,:)
!!
!!rho_s(1,2:ny-1,:) = 2 * rho_s(1,2:ny-1,:)
!!rho_s(2:nx-1,1,:) = 2 * rho_s(2:nx-1,1,:)
!!rho_s(2:nx-1,ny,:) = 2 * rho_s(2:nx-1,ny,:)
!!rho_s(nx,2:ny-1,:) = 2 * rho_s(nx,2:ny-1,:)
!!Cell_volume_bjw = hx(1) * hx(2) !1.0
!
!!$ == descripted by ZWZ for using different cell volume ==//
!
!
!! Finally normalize the density and find the total charge density
!DO isp = 1, ispe_tot
!
!		DO jj = 1, ny
!			DO ii = 1, nx
!
!                n = (jj+(ii-1)*ny)
!
!                IF (Ek_tot(ii,jj,isp)/=0 .AND. rho_s(ii,jj,isp)/=0) THEN
!                    Ek_s(ii,jj,isp) = Ek_tot(ii,jj,isp)/ rho_s(ii,jj,isp)       !$ ab.ZWZ 2022/2/9
!                ELSE IF (Ek_tot(ii,jj,isp)/=0 .AND. rho_s(ii,jj,isp)==0) THEN
!                    print*,'someting wrong'
!                    stop
!                END IF
!                
!				!rho_s(ii,jj,isp) = affp(isp)*qs(isp)*rho_s(ii,jj,isp)/Cell_volume(n)
!                !rho_s(ii,jj,isp) = affp_bjw(isp)*qs(isp)*rho_s(ii,jj,isp)/L_ref/L_ref/n_ref/Cell_volume_bjw
!                !$ =========================== mb.ZWZ 2021/7/23 =============================== \\
!                rho_n(ii,jj,isp) = rho_s(ii,jj,isp)
!                !IF (delta == 0) THEN
!                !    rho_s(ii,jj,isp) = affp_bjw(isp)*qs(isp)*rho_s(ii,jj,isp)/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj) 
!                !    Ek_tot(ii,jj,isp) = affp_bjw(isp)*Ek_tot(ii,jj,isp)/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj) 
!                !ELSEIF (delta == 1) THEN
!                !    rho_s(ii,jj,isp) = affp_bjw(isp)*qs(isp)*rho_s(ii,jj,isp)/L_ref/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj)  
!                !    Ek_tot(ii,jj,isp) = affp_bjw(isp)*Ek_tot(ii,jj,isp)/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj) 
!                !ENDIF
!               
!                IF (delta == 0) THEN
!                    rho_s(ii,jj,isp) = ParticleGlobal(isp-1)%Weight*qs(isp)*rho_s(ii,jj,isp)/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj) 
!                    Ek_tot(ii,jj,isp) = ParticleGlobal(isp-1)%Weight*Ek_tot(ii,jj,isp)/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj) 
!                ELSEIF (delta == 1) THEN
!                    If (ParticleGlobal(isp-1)%UnequalWeightFlag) Then
!                        rho_s(ii,jj,isp) = qs(isp)*rho_s(ii,jj,isp)/L_ref/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj)  
!                        Ek_tot(ii,jj,isp) = Ek_tot(ii,jj,isp)/L_ref/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj) 
!                    Else
!                        rho_s(ii,jj,isp) = ParticleGlobal(isp-1)%Weight*qs(isp)*rho_s(ii,jj,isp)/L_ref/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj)  
!                        Ek_tot(ii,jj,isp) = ParticleGlobal(isp-1)%Weight*Ek_tot(ii,jj,isp)/L_ref/L_ref/L_ref/n_ref/Cell_volume_zwz (ii,jj) 
!                    Endif
!                ENDIF
!                !$ =========================== mb.ZWZ 2021/7/23 =============================== //
!                                   
!				rho(ii,jj) = rho(ii,jj) + rho_s(ii,jj,isp)
!              
!
!			END DO
!		END DO
!END DO
!!=========ZWZ brother Old Code=========









END

