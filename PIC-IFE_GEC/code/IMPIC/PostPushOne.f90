Subroutine PostPushOne(PO,isp,TimeMove, IvelFlag, IposFlag, delta, OriPosi, detaV)
    
    Use Domain_2D
    Use Particle_2D, Only: qm
    Use Field_2D
    Use IMPIC_Data_2D
    Use ModuleParticleOne
    !============LY add for SIDG-PIC couple, 2022-7-25============
    Use IFE_Data
    !============LY add for SIDG-PIC couple, 2022-7-25============
    Implicit none
    
    !Class(ParticleOne), intent(inout) :: PO
    Type(ParticleOne), intent(inout) :: PO
    Integer(4), intent(in) :: isp
    Real(8), intent(in) :: TimeMove
    Integer(4), intent(inout) :: IvelFlag, IposFlag
    Integer(4), intent(in) :: delta
    Real(8), intent(out) :: OriPosi(3)
    Real(8), intent(out) :: detaV(3)
    
    Real(8) :: xp, yp
    Real(8) :: dx, dy, xcellmdx, ycellmdy
    Real(8) :: Rp, R1, R2, den
    Real(8) :: P1, P2, P3, P4
    Integer(4) :: i, j
    
    Real(8) :: f
    Real(8) :: SVelocity(1:3)=0.d0
    Real(8) :: A_x, A_y, A_z
    Real(8), Save :: deta_vx, deta_vy, deta_vz
    Real(8) :: X, Y, Z, Theta
    Real(8) :: zefield, refield, tefield
    Real(8) :: xefield, yefield
    Real(8) :: zbfield, rbfield, tbfield
    Real(8) :: xbfield, ybfield
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    Real(8) :: out_object_beta, in_object_beta
    Integer :: traverse_flag
    Integer :: n_element_old, n_element_bro
    Real(8) :: part_ele_left, part_ele_right, part_ele_bottom, part_ele_top
    Real(8) :: bottom_top_centre, left_right_centre, delta_1, delta_2
    Integer :: edge_bro_element(4), part_edge_count(9)
    Integer :: edge_count
    Real(8) :: hx_partition, hy_partition
    Real(8) :: W1, W2, W3, W4
    Real(8) :: ex_part, ey_part
    Integer :: trial_basis_type
    Real(8) :: vertices(2,4)
    Integer :: index_ele, piece_flag, n
    Real(8) :: beta1, beta2
    Integer :: information_vector_1(18)
    Real(8) :: information_vector_2(8)
    Real(8) :: Ex_basis, Ey_basis
    Real(8) :: Efield_old(3), Efield_bro(3), Efield(3)
    Real(8) :: Bfield_old(3), Bfield_bro(3)
    Integer :: edge_boundary, edge_bro_count, edge_temp_count, edge_self_element, edge_self_element_min
    Integer :: number_ele_column
    Integer :: n_element_bro_initial
    Real(8), Parameter :: PI_LY	= 3.14159265358979D0
    Integer, Dimension(:), Allocatable :: EdgeArray
    Real(8), Dimension(:,:), Allocatable :: HE_particle
    
    Allocate(HE_particle(5, Size(HE,2)))
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    
    
    OriPosi = (/PO%X,PO%Y,PO%Z/)
    
    If(IvelFlag == 1) Then
      !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
      !==========Part 0 : Initialization==========
      edge_bro_element(1:4) = 0
      number_ele_column = ny - 1
      !out_object_beta = Global_Beta(1)
      !in_object_beta = Global_Beta(2)
      If (Allocated(EdgeArray)) Then
        Deallocate(EdgeArray)
      End If
      !==========Part 0 : Initialization==========
      
      !==========Part 1 : Positioning Particle and Moving==========
      traverse_flag = 0
      
      xp = (PO%X - Vert_o(1))*hxi(1)
      i = INT(xp)

      yp = (PO%Y - Vert_o(2))*hxi(2)
      j = INT(yp)
      
      If (j == ny) Then
        j = j - 1
      End If
      
      If (i == nx) Then
        i = i -1
      End If
      
      n_element_old = (j + (i-1)*(ny-1))  !Initial element on Tier1: no refinement.
      
      !=========LY modification for Speeding Particle Positioning, 2022-6-22=========
      Do While (CellMesh(n_element_old)%isSplitted == 1)
        !There denote the initial element is refinement element, it will be refine.
        !The following four variables denote the distance between particle and DG element boundary.
        part_ele_left = PO%X - CellMesh(n_element_old)%Boundary(1)    !Distance between particle and element's left boundary.
        part_ele_right = PO%X - CellMesh(n_element_old)%Boundary(2)   !Distance between particle and element's right boundary.
        part_ele_bottom = PO%Y - CellMesh(n_element_old)%Boundary(3)  !Distance between particle and element's bottom boundary.
        part_ele_top = PO%Y - CellMesh(n_element_old)%Boundary(4)     !Distance between particle and element's top boundary.

        !The following two variables 'delta_1' and 'delta_2' denote the distance between particle and DG element centre.
        left_right_centre = (CellMesh(n_element_old)%Boundary(1)+CellMesh(n_element_old)%Boundary(2)) / 2.0
        bottom_top_centre = (CellMesh(n_element_old)%Boundary(3)+CellMesh(n_element_old)%Boundary(4)) / 2.0
        delta_1 = PO%Y - bottom_top_centre  !Distance between particle and element's centre of bottom and top.
        delta_2 = PO%X - left_right_centre  !Distance between particle and element's centre of left and right.

        If (ABS(part_ele_left)<SmallValue .OR. ABS(part_ele_right)<SmallValue .OR. &
            ABS(part_ele_bottom)<SmallValue .OR. ABS(part_ele_top)<SmallValue) Then

          If (ABS(part_ele_left)<SmallValue .AND. ABS(part_ele_top)<SmallValue) Then
            !Condition A: part_ele_left = 0 and part_ele_top = 0
            n_element_old = CellMesh(n_element_old)%Child(2)  !Left-Top child element
          Elseif (ABS(part_ele_right)<SmallValue .AND. ABS(part_ele_top)<SmallValue) Then
            !Condition B: part_ele_right = 0 and part_ele_top = 0
            n_element_old = CellMesh(n_element_old)%Child(4)  !Right-Top child element
          Elseif (ABS(part_ele_left)<SmallValue .AND. ABS(part_ele_bottom)<SmallValue) Then
            !Condition C: part_ele_left = 0 and part_ele_bottom = 0
            n_element_old = CellMesh(n_element_old)%Child(1)  !Left-Bottom child element
          Elseif (ABS(part_ele_right)<SmallValue .AND. ABS(part_ele_bottom)<SmallValue) Then
            !Condition D: part_ele_right = 0 and part_ele_bottom = 0
            n_element_old = CellMesh(n_element_old)%Child(3)  !Right-Bottom child element
          Else
            !There denote the particle locate on the edge of element, but not corner.
            If (ABS(part_ele_left)<SmallValue) Then
              !Particle locate on the left edge but not endpoints.
              If (ABS(CellMesh(n_element_old)%Boundary(1) - dxmin)<SmallValue) Then
                !Particle locate on the left boundary edge but not endpoints.
                !We choose the edge set of self-initial element.
                traverse_flag = 1
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                Exit  !Jump out of the present 'DO WHILE' cycle.
              Else
                !Particle locate on the left inner edge but not endpoints.
                !We choose the edge set of self-initial element and bro-initial element.
                traverse_flag = 1
                n_element_bro_initial = CellMesh(n_element_old)%FinalParent - number_ele_column
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)+EdgeParent(1,n_element_bro_initial)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                j = 1
                Do i = EdgeParent(1,CellMesh(n_element_old)%FinalParent)+1, Size(EdgeArray)
                  EdgeArray(i) = EdgeParent(1+j,n_element_bro_initial)
                  j = j + 1
                End Do
                Exit  !Jump out of the present 'DO WHILE' cycle.
              End If

            Elseif (ABS(part_ele_right)<SmallValue) Then
              !Particle locate on the right edge but not endpoints.
              If (ABS(CellMesh(n_element_old)%Boundary(2) - dxmax)<SmallValue) Then
                !Particle locate on the right boundary edge but not endpoints.
                !We choose the edge set of self-initial element.
                traverse_flag = 1
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                Exit
              Else
                !Particle locate on the right inner edge but not endpoints.
                !We choose the edge set of self-initial element and bro-initial element.
                traverse_flag = 1
                n_element_bro_initial = CellMesh(n_element_old)%FinalParent + number_ele_column
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)+EdgeParent(1,n_element_bro_initial)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                j = 1
                Do i = EdgeParent(1,CellMesh(n_element_old)%FinalParent)+1, Size(EdgeArray)
                  EdgeArray(i) = EdgeParent(1+j,n_element_bro_initial)
                  j = j + 1
                End Do
                Exit  !Jump out of the present 'DO WHILE' cycle.
              End If

            Elseif (ABS(part_ele_bottom)<SmallValue) Then
              !Particle locate on the bottom edge but not endpoints.
              If (ABS(CellMesh(n_element_old)%Boundary(3) - dymin)<SmallValue) Then
                !Particle locate on the bottom boundary edge but not endpoints.
                !We choose the edge set of self-initial element.
                traverse_flag = 1
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                Exit
              Else
                !Particle locate on the bottom inner edge but not endpoints.
                !We choose the edge set of self-initial element and bro-initial element.
                traverse_flag = 1
                n_element_bro_initial = CellMesh(n_element_old)%FinalParent - 1
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)+EdgeParent(1,n_element_bro_initial)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                j = 1
                Do i = EdgeParent(1,CellMesh(n_element_old)%FinalParent)+1, Size(EdgeArray)
                  EdgeArray(i) = EdgeParent(1+j,n_element_bro_initial)
                  j = j + 1
                End Do
                Exit  !Jump out of the present 'DO WHILE' cycle.
              End If

            Elseif (ABS(part_ele_top)<SmallValue) Then
              !Particle locate on the top edge but not endpoints.
              If (ABS(CellMesh(n_element_old)%Boundary(4) - dymax)<SmallValue) Then
                !Particle locate on the top boundary edge but not endpoints.
                !We choose the edge set of self-initial element.
                traverse_flag = 1
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                Exit
              Else
                !Particle locate on the top inner edge but not endpoints.
                !We choose the edge set of self-initial element and bro-initial element.
                traverse_flag = 1
                n_element_bro_initial = CellMesh(n_element_old)%FinalParent + 1
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)+EdgeParent(1,n_element_bro_initial)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                j = 1
                Do i = EdgeParent(1,CellMesh(n_element_old)%FinalParent)+1, Size(EdgeArray)
                  EdgeArray(i) = EdgeParent(1+j,n_element_bro_initial)
                  j = j + 1
                End Do
                Exit  !Jump out of the present 'DO WHILE' cycle.
              End If
            End If
          End If

        Elseif (ABS(delta_1)<SmallValue .OR. ABS(delta_2)<SmallValue) Then

          If (ABS(delta_1)<SmallValue .AND. ABS(delta_2)<SmallValue) Then
            !Condition G: delta_1 = 0 and delta_2 = 0
            n_element_old = CellMesh(n_element_old)%Child(1)  !In here, It is ok whatever we choose any child element.
          Else
            !Particle locate on the left-right-centre or bottom-top-centre line but not endpoints and intersection.
            !We choose the edge set of self-initial element.
            traverse_flag = 1
            Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)))
            EdgeArray = 0
            Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
              EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
            End Do
            Exit  !Jump out of the present 'DO WHILE' cycle.
          End If

        Else
          !There denote the particle locate inside the DG element.
          !We need to discuss the particle locate inside which one child element.
          If (delta_1 > SmallValue) Then   !Particle locate inside left-top child element or right-top child element.
            If (delta_2 > SmallValue) Then   !Particle locate inside right-top child element.
              n_element_old = CellMesh(n_element_old)%Child(4)
            Elseif (delta_2 < -SmallValue) Then   !Particle locate inside left-top child element.
              n_element_old = CellMesh(n_element_old)%Child(2)
            Elseif (ABS(delta_2) < SmallValue) Then
              Write(6,*) 'The delta2 = PO%X - left_right_centre is 0, but traverse_flag = 0, ERROR--PostPush'
              Stop
            End If
          Elseif (delta_1 < -SmallValue) Then   !Particle locate inside left-bottom child element or right-bottom child element.
            If (delta_2 > SmallValue) Then   !Particle locate inside right-bottom child element.
              n_element_old = CellMesh(n_element_old)%Child(3)
            Elseif (delta_2 < -SmallValue) Then   !Particle locate inside left-bottom child element.
              n_element_old = CellMesh(n_element_old)%Child(1)
            Elseif (ABS(delta_2) < SmallValue) Then
              Write(6,*) 'The delta2 = PO%X - left_right_centre is 0, but traverse_flag = 0, ERROR--PostPush'
              Stop
            End If
          Elseif (ABS(delta_1) < SmallValue) Then
            Write(6,*) 'The delta1 = PO%Y - bottom_top_centre is 0, but traverse_flag = 0, ERROR--PostPush'
            Stop
          End If

        End If
      End Do
      
      If (CellMesh(n_element_old)%isSplitted==0 .AND. CellMesh(n_element_old)%Finalindex/=0) Then
        !There denote the element is no-refined element.

        !The following four variables denote the distance between particle and CG element boundary.
        part_ele_left = PO%X - CellMesh(n_element_old)%Boundary(1)    !Distance between particle and element's left boundary.
        part_ele_right = PO%X - CellMesh(n_element_old)%Boundary(2)   !Distance between particle and element's right boundary.
        part_ele_bottom = PO%Y - CellMesh(n_element_old)%Boundary(3)  !Distance between particle and element's bottom boundary.
        part_ele_top = PO%Y - CellMesh(n_element_old)%Boundary(4)     !Distance between particle and element's top boundary.

        If (ABS(part_ele_left)<SmallValue .OR. ABS(part_ele_right)<SmallValue .OR. &
            ABS(part_ele_bottom)<SmallValue .OR. ABS(part_ele_top)<SmallValue) Then

          If (ABS(part_ele_left)<SmallValue .AND. ABS(part_ele_top)<SmallValue) Then
            !Condition A: part_ele_left = 0 and part_ele_top = 0
            n_element_old = CellMesh(n_element_old)%Finalindex
          Elseif (ABS(part_ele_right)<SmallValue .AND. ABS(part_ele_top)<SmallValue) Then
            !Condition B: part_ele_right = 0 and part_ele_top = 0
            n_element_old = CellMesh(n_element_old)%Finalindex
          Elseif (ABS(part_ele_left)<SmallValue .AND. ABS(part_ele_bottom)<SmallValue) Then
            !Condition C: part_ele_left = 0 and part_ele_bottom = 0
            n_element_old = CellMesh(n_element_old)%Finalindex
          Elseif (ABS(part_ele_right)<SmallValue .AND. ABS(part_ele_bottom)<SmallValue) Then
            !Condition D: part_ele_right = 0 and part_ele_bottom = 0
            n_element_old = CellMesh(n_element_old)%Finalindex
          Else
            !There denote the particle locate on the edge of element, but not corner.
            If (ABS(part_ele_left)<SmallValue) Then
              !Particle locate on the left edge but no corner.
              If (ABS(CellMesh(n_element_old)%Boundary(1)-dxmin)<SmallValue) Then
                !Particle locate on the left boundary edge but not endpoints.
                !We choose the self-element Finalindex.
                n_element_old = CellMesh(n_element_old)%Finalindex
              Else
                !Particle locate on the left inner edge but not endpoints.
                !We choose the edge set of self-initial element and bro-initial element.
                traverse_flag = 1
                n_element_bro_initial = CellMesh(n_element_old)%FinalParent - number_ele_column
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)+EdgeParent(1,n_element_bro_initial)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                j = 1
                Do i = EdgeParent(1,CellMesh(n_element_old)%FinalParent)+1, Size(EdgeArray)
                  EdgeArray(i) = EdgeParent(1+j,n_element_bro_initial)
                  j = j + 1
                End Do
              End If

            Elseif (ABS(part_ele_right)<SmallValue) Then
              !Particle locate on the right edge but no corner.
              If (ABS(CellMesh(n_element_old)%Boundary(2)-dxmax)<SmallValue) Then
                !Particle locate on the right boundary edge but not endpoints.
                !We choose the self-element Finalindex.
                n_element_old = CellMesh(n_element_old)%Finalindex
              Else
                !Particle locate on the right inner edge but not endpoints.
                !We choose the edge set of self-initial element and bro-initial element.
                traverse_flag = 1
                n_element_bro_initial = CellMesh(n_element_old)%FinalParent + number_ele_column
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)+EdgeParent(1,n_element_bro_initial)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                j = 1
                Do i = EdgeParent(1,CellMesh(n_element_old)%FinalParent)+1, Size(EdgeArray)
                  EdgeArray(i) = EdgeParent(1+j,n_element_bro_initial)
                  j = j + 1
                End Do
              End If

            Elseif (ABS(part_ele_bottom)<SmallValue) Then
              !Particle locate on the bottom edge but no corner.
              If (ABS(CellMesh(n_element_old)%Boundary(3)-dymin)<SmallValue) Then
                !Particle locate on the bottom boundary edge but not endpoints.
                !We choose the self-element Finalindex.
                n_element_old = CellMesh(n_element_old)%Finalindex
              Else
                !Particle locate on the bottom inner edge but not endpoints.
                !We choose the edge set of self-initial element and bro-initial element.
                traverse_flag = 1
                n_element_bro_initial = CellMesh(n_element_old)%FinalParent - 1
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)+EdgeParent(1,n_element_bro_initial)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                j = 1
                Do i = EdgeParent(1,CellMesh(n_element_old)%FinalParent)+1, Size(EdgeArray)
                  EdgeArray(i) = EdgeParent(1+j,n_element_bro_initial)
                  j = j + 1
                End Do
              End If

            Elseif (ABS(part_ele_top)<SmallValue) Then
              !Particle locate on the top edge but no corner.
              If (ABS(CellMesh(n_element_old)%Boundary(4)-dymax)<SmallValue) Then
                !Particle locate on the top boundary edge but not endpoints.
                !We choose the self-element Finalindex.
                n_element_old = CellMesh(n_element_old)%Finalindex
              Else
                !Particle locate on the top inner edge but not endpoints.
                !We choose the edge set of self-initial element and bro-initial element.
                traverse_flag = 1
                n_element_bro_initial = CellMesh(n_element_old)%FinalParent + 1
                Allocate(EdgeArray(EdgeParent(1,CellMesh(n_element_old)%FinalParent)+EdgeParent(1,n_element_bro_initial)))
                EdgeArray = 0
                Do i = 1, EdgeParent(1,CellMesh(n_element_old)%FinalParent)
                  EdgeArray(i) = EdgeParent(1+i,CellMesh(n_element_old)%FinalParent)
                End Do
                j = 1
                Do i = EdgeParent(1,CellMesh(n_element_old)%FinalParent)+1, Size(EdgeArray)
                  EdgeArray(i) = EdgeParent(1+j,n_element_bro_initial)
                  j = j + 1
                End Do
              End If

            End If
          End If

        Else
          n_element_old = CellMesh(n_element_old)%Finalindex
        End If
      End If
      !=========LY modification for Speeding Particle Positioning, 2022-6-22=========
      
      If (traverse_flag == 1) Then
        !There denote the particle locate on edge. We need to traverse between particle and every edge.

        part_edge_count(1:9) = 0   !The variable 'part_edge_count' is record intersection number between particle and edge.
        edge_count = 2    !The variable 'edge_count' is temp counter.

        !=========LY modification for Speeding Particle Positioning, 2022-6-22=========
        If (Allocated(HE_particle)) Then
          Deallocate(HE_particle)
        End If
        Allocate(HE_particle(5,Size(EdgeArray)))
        HE_particle = 0.0

        Do j = 1, Size(EdgeArray)
          HE_particle(1,j) = DSQRT((PO%X-HP(1,HE(1,EdgeArray(j))))**2 + &
                                   (PO%Y-HP(2,HE(1,EdgeArray(j))))**2)
          HE_particle(2,j) = DSQRT((PO%X-HP(1,HE(2,EdgeArray(j))))**2 + &
                                   (PO%Y-HP(2,HE(2,EdgeArray(j))))**2)
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
        If (part_edge_count(1) == 1) Then
          !The particle locate the edge between CG(0.5) and DG(0.5).√

          If (HE(6, part_edge_count(2)) == 0) Then  !There denote the DG edge is exterior boundary edge. We need protect program from array.

            n_element_old = HE(5, part_edge_count(2))
            traverse_flag = 0

          Elseif (HE(6, part_edge_count(2)) /= 0) Then  !There denote the DG edge is interior edge.

            If (HT(5, HE(5,part_edge_count(2)))/=0 .AND. HT(5, HE(6,part_edge_count(2)))==0) Then

              !======The old element part======
              n_element_old = HE(5, part_edge_count(2)) !The 0.5 particle distribute DG element.

              hx_partition = HP(1, HT(2,n_element_old)) - HP(1, HT(1,n_element_old))
              hy_partition = HP(2, HT(4,n_element_old)) - HP(2, HT(1,n_element_old))

              dx = (PO%X - HP(1, HT(1,n_element_old))) / hx_partition
              dy = (PO%Y - HP(2, HT(1,n_element_old))) / hy_partition
              xcellmdx = 1.0 - dx
              ycellmdy = 1.0 - dy

              If (delta_global == 0) Then  !2D Cartesian coordinates.
                W1 = xcellmdx * ycellmdy
                W2 = dx       * ycellmdy
                W3 = dx       * dy
                W4 = xcellmdx * dy
              Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
                R1 = dymin + HP(2,HT(1,n_element_old))
                R2 = dymin + HP(2,HT(4,n_element_old))
                Rp = PO%Y
                den = R2*R2 - R1*R1

                !Here, we not consider the unequal weights for the moment, 2022-4-8.
                W1 = xcellmdx * (R2*R2-Rp*Rp) / den
                W2 = dx       * (R2*R2-Rp*Rp) / den
                W3 = dx       * (Rp*Rp-R1*R1) / den
                W4 = xcellmdx * (Rp*Rp-R1*R1) / den
              End If

              !Firstly, we need to calculate the Efield of old element.
              If (element_index(n_element_old) > 0) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
                ex_part = 0.0
                ey_part = 0.0
                trial_basis_type = 1
                vertices = HP(1:2, HT(1:4, n_element_old))
                index_ele = element_index(n_element_old)
                beta1 = information_2(1, index_ele)
                beta2 = information_2(2, index_ele)
                information_vector_1 = information_1(1:18, index_ele)
                information_vector_2 = information_2(1:8, index_ele)

                If (ABS(beta1 - out_object_beta) < SmallValue) Then
                  piece_flag = 1
                Elseif (ABS(beta2 - out_object_beta) < SmallValue) Then
                  piece_flag = 2
                End If

                Do n = 1, 4
                  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                  piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_old)
                  ex_part = ex_part + Phi(HT(n, n_element_old), 1) * Ex_basis

                  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                  piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_old)
                  ey_part = ey_part + Phi(HT(n, n_element_old), 1) * Ey_basis
                End Do

                Efield_old(1) = -ex_part*0.5
                Efield_old(2) = -ey_part*0.5

              Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
                Efield_old(1) = (W1 * efx(HT(1,n_element_old),1) + W2 * efx(HT(2,n_element_old),1) + &
                                 W3 * efx(HT(3,n_element_old),1) + W4 * efx(HT(4,n_element_old),1)) * 0.5
                Efield_old(2) = (W1 * efy(HT(1,n_element_old),1) + W2 * efy(HT(2,n_element_old),1) + &
                                 W3 * efy(HT(3,n_element_old),1) + W4 * efy(HT(4,n_element_old),1)) * 0.5
              End If
              !=========LY modification for 2D3V PIC model, 2022-5-27=========
              Efield_old(3) = (W1 * efz(HT(1,n_element_old),1) + W2 * efz(HT(2,n_element_old),1) + &
                               W3 * efz(HT(3,n_element_old),1) + W4 * efz(HT(4,n_element_old),1)) * 0.5
              !=========LY modification for 2D3V PIC model, 2022-5-27=========

              !Secondly, we need to calculate the Bfield of old element.
              Bfield_old(1) = (W1 * bfx(HT(1,n_element_old),1) + W2 * bfx(HT(2,n_element_old),1) + &
                               W3 * bfx(HT(3,n_element_old),1) + W4 * bfx(HT(4,n_element_old),1)) * 0.5
              Bfield_old(2) = (W1 * bfy(HT(1,n_element_old),1) + W2 * bfy(HT(2,n_element_old),1) + &
                               W3 * bfy(HT(3,n_element_old),1) + W4 * bfy(HT(4,n_element_old),1)) * 0.5
              Bfield_old(3) = (W1 * bfz(HT(1,n_element_old),1) + W2 * bfz(HT(2,n_element_old),1) + &
                               W3 * bfz(HT(3,n_element_old),1) + W4 * bfz(HT(4,n_element_old),1)) * 0.5
              !======The old element part======

              !======The bro element part======
              n_element_bro = HE(6, part_edge_count(2)) !The 0.5 particle distribute CG element.

              hx_partition = HP(1, HT(2,n_element_bro)) - HP(1, HT(1,n_element_bro))
              hy_partition = HP(2, HT(4,n_element_bro)) - HP(2, HT(1,n_element_bro))

              dx = (PO%X - HP(1, HT(1,n_element_bro))) / hx_partition
              dy = (PO%Y - HP(2, HT(1,n_element_bro))) / hy_partition
              xcellmdx = 1.0 - dx
              ycellmdy = 1.0 - dy

              If (delta_global == 0) Then  !2D Cartesian coordinates.
                W1 = xcellmdx * ycellmdy
                W2 = dx       * ycellmdy
                W3 = dx       * dy
                W4 = xcellmdx * dy
              Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
                R1 = dymin + HP(2,HT(1,n_element_bro))
                R2 = dymin + HP(2,HT(4,n_element_bro))
                Rp = PO%Y
                den = R2*R2 - R1*R1

                !Here, we not consider the unequal weights for the moment, 2022-4-8.
                W1 = xcellmdx * (R2*R2-Rp*Rp) / den
                W2 = dx       * (R2*R2-Rp*Rp) / den
                W3 = dx       * (Rp*Rp-R1*R1) / den
                W4 = xcellmdx * (Rp*Rp-R1*R1) / den
              End If

              !Firstly, we need to calculate the Efield of bro element.
              If (element_index(n_element_bro) > 0) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
                ex_part = 0.0
                ey_part = 0.0
                trial_basis_type = 1
                vertices = HP(1:2, HT(1:4, n_element_bro))
                index_ele = element_index(n_element_bro)
                beta1 = information_2(1, index_ele)
                beta2 = information_2(2, index_ele)
                information_vector_1 = information_1(1:18, index_ele)
                information_vector_2 = information_2(1:8, index_ele)

                If (Abs(beta1 - out_object_beta) < SmallValue) Then
                  piece_flag = 1
                Elseif (Abs(beta2 - out_object_beta) < SmallValue) Then
                  piece_flag = 2
                End If

                Do n = 1, 4
                  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                  piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_bro)
                  ex_part = ex_part + Phi(HT(n, n_element_bro), 1) * Ex_basis

                  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                  piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_bro)
                  ey_part = ey_part + Phi(HT(n, n_element_bro), 1) * Ey_basis
                End Do

                Efield_bro(1) = -ex_part*0.5
                Efield_bro(2) = -ey_part*0.5

              Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
                Efield_bro(1) = (W1 * efx(HT(1,n_element_bro),1) + W2 * efx(HT(2,n_element_bro),1) + &
                                 W3 * efx(HT(3,n_element_bro),1) + W4 * efx(HT(4,n_element_bro),1)) * 0.5
                Efield_bro(2) = (W1 * efy(HT(1,n_element_bro),1) + W2 * efy(HT(2,n_element_bro),1) + &
                                 W3 * efy(HT(3,n_element_bro),1) + W4 * efy(HT(4,n_element_bro),1)) * 0.5
              End If
              !=========LY modification for 2D3V PIC model, 2022-5-27=========
              Efield_bro(3) = (W1 * efz(HT(1,n_element_bro),1) + W2 * efz(HT(2,n_element_bro),1) + &
                               W3 * efz(HT(3,n_element_bro),1) + W4 * efz(HT(4,n_element_bro),1)) * 0.5
              !=========LY modification for 2D3V PIC model, 2022-5-27=========

              !Secondly, we need to calculate the Bfield of bro element.
              Bfield_bro(1) = (W1 * bfx(HT(1,n_element_bro),1) + W2 * bfx(HT(2,n_element_bro),1) + &
                               W3 * bfx(HT(3,n_element_bro),1) + W4 * bfx(HT(4,n_element_bro),1)) * 0.5
              Bfield_bro(2) = (W1 * bfy(HT(1,n_element_bro),1) + W2 * bfy(HT(2,n_element_bro),1) + &
                               W3 * bfy(HT(3,n_element_bro),1) + W4 * bfy(HT(4,n_element_bro),1)) * 0.5
              Bfield_bro(3) = (W1 * bfz(HT(1,n_element_bro),1) + W2 * bfz(HT(2,n_element_bro),1) + &
                               W3 * bfz(HT(3,n_element_bro),1) + W4 * bfz(HT(4,n_element_bro),1)) * 0.5
              !======The bro element part======

              If (delta_global == 0) Then !2D Cartesian coordinates.
                xefield = Efield_old(1) + Efield_bro(1)
                yefield = Efield_old(2) + Efield_bro(2)
                zefield = Efield_old(3) + Efield_bro(3)

                A_x = xefield * qm(isp+1)
                A_y = yefield * qm(isp+1)
                A_z = zefield * qm(isp+1)
              Elseif (delta_global == 1) Then !2D axis-symmetric coordinates.
                !$ axisymmetric e field
                zefield = Efield_old(1) + Efield_bro(1)
                refield = Efield_old(2) + Efield_bro(2)
                tefield = Efield_old(3) + Efield_bro(3)
                !$ Convert axisymmetric efield to Cartesian efield
                Theta = PO%Z
                xefield = refield*DCOS(Theta) - tefield*DSIN(Theta)
                yefield = refield*DSIN(Theta) + tefield*DCOS(Theta)

                A_x = zefield * qm(isp+1)
                A_y = xefield * qm(isp+1)
                A_z = yefield * qm(isp+1)
              End If

              A_n_plus_1(1) = A_x
              A_n_plus_1(2) = A_y
              A_n_plus_1(3) = A_z

              !Thirdly, we need to calculate the Efield force and Bfield force of the particle.
              If (Bfiled_index) Then  !with magnetic field.
                If (delta_global == 0) Then  !2D Cartesian coordinates.
                  xbfield = Bfield_old(1) + Bfield_bro(1)
                  ybfield = Bfield_old(2) + Bfield_bro(2)
                  zbfield = Bfield_old(3) + Bfield_bro(3)

                  Omega(1) = xbfield * qm(isp+1) * 0.5 * TimeMove
                  Omega(2) = ybfield * qm(isp+1) * 0.5 * TimeMove
                  Omega(3) = zbfield * qm(isp+1) * 0.5 * TimeMove
                Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
                  !$ axisymmetric b field
                  zbfield = Bfield_old(1) + Bfield_bro(1)
                  rbfield = Bfield_old(2) + Bfield_bro(2)
                  tbfield = Bfield_old(3) + Bfield_bro(3)
                  !$ Convert axisymmetric bfield to Cartesian bfield
                  Theta = PO%Z
                  xbfield = rbfield*DCOS(Theta) - tbfield*DSIN(Theta)
                  ybfield = rbfield*DSIN(Theta) + tbfield*DCOS(Theta)
                  
                  Omega(1) = zbfield * qm(isp+1) * 0.5 * TimeMove
                  Omega(2) = xbfield * qm(isp+1) * 0.5 * TimeMove
                  Omega(3) = ybfield * qm(isp+1) * 0.5 * TimeMove
                End If

                OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
                TransB(1,1)=1.d0+Omega(1)*Omega(1)
                TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
                TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
                TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
                TransB(2,2)=1.d0+Omega(2)*Omega(2)
                TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
                TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
                TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
                TransB(3,3)=1.d0+Omega(3)*Omega(3)
                TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT

                SVelocity(1) = A_n_plus_1(1) * 0.5 * TimeMove
                SVelocity(2) = A_n_plus_1(2) * 0.5 * TimeMove
                SVelocity(3) = A_n_plus_1(3) * 0.5 * TimeMove

                deta_vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
                deta_vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
                deta_vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)

                PO%Vx = PO%Vx + deta_vx
                PO%Vy = PO%Vy + deta_vy
                PO%Vz = PO%Vz + deta_vz

              Else  !without magnetic field.
                deta_vx = 0.5*TimeMove*A_n_plus_1(1)
                deta_vy = 0.5*TimeMove*A_n_plus_1(2)
                deta_vz = 0.5*TimeMove*A_n_plus_1(3)

                PO%Vx = PO%Vx + deta_vx
                PO%Vy = PO%Vy + deta_vy
                PO%Vz = PO%Vz + deta_vz
              End If

              A_bar_n_minus_1(1) = PO%Ax
              A_bar_n_minus_1(2) = PO%Ay
              A_bar_n_minus_1(3) = PO%Az

              A_bar_n(1) = 0.5*(A_bar_n_minus_1(1) + A_n_plus_1(1))
              A_bar_n(2) = 0.5*(A_bar_n_minus_1(2) + A_n_plus_1(2))
              A_bar_n(3) = 0.5*(A_bar_n_minus_1(3) + A_n_plus_1(3))

              PO%Ax = A_bar_n(1)
              PO%Ay = A_bar_n(2)
              PO%Az = A_bar_n(3)

            Else
              Write(6,*) 'Program error, pelase check : part_edge_count(1) == 1!'
              Stop
            End If
          End If

        Elseif (part_edge_count(1) == 2) Then
          !The particle locate the edge between CG(0.0) and DG(1.0).√
          !The particle locate the edge between DG(0.5) and DG(0.5).√

          If (HE(5, part_edge_count(2)) == HE(5, part_edge_count(3))) Then
            !The particle locate the edge between CG(0.0) and DG(1.0).(3CG and 1DG, the particle locate the corner node)
            n_element_old = HE(5,part_edge_count(2))
            traverse_flag = 0

          Elseif (HE(5, part_edge_count(2)) /= HE(5, part_edge_count(3))) Then
            !The particle locate the edge between DG(0.5) and DG(0.5).

            !======The old element part======
            n_element_old = HE(5, part_edge_count(2)) !The 0.5 particle distribute DG element

            hx_partition = HP(1, HT(2,n_element_old)) - HP(1, HT(1,n_element_old))
            hy_partition = HP(2, HT(4,n_element_old)) - HP(2, HT(1,n_element_old))

            dx = (PO%X - HP(1, HT(1,n_element_old))) / hx_partition
            dy = (PO%Y - HP(2, HT(1,n_element_old))) / hy_partition
            xcellmdx = 1.0 - dx
            ycellmdy = 1.0 - dy

            If (delta_global == 0) Then  !2D Cartesian coordinates.
              W1 = xcellmdx * ycellmdy
              W2 = dx       * ycellmdy
              W3 = dx       * dy
              W4 = xcellmdx * dy
            Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
              R1 = dymin + HP(2,HT(1,n_element_old))
              R2 = dymin + HP(2,HT(4,n_element_old))
              Rp = PO%Y
              den = R2*R2 - R1*R1

              !Here, we not consider the unequal weights for the moment, 2022-4-8.
              W1 = xcellmdx * (R2*R2-Rp*Rp) / den
              W2 = dx       * (R2*R2-Rp*Rp) / den
              W3 = dx       * (Rp*Rp-R1*R1) / den
              W4 = xcellmdx * (Rp*Rp-R1*R1) / den
            End If

            !Firstly, we need to calculate the Efield of old element.
            If (element_index(n_element_old) > 0) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
              ex_part = 0.0
              ey_part = 0.0
              trial_basis_type = 1
              vertices = HP(1:2, HT(1:4, n_element_old))
              index_ele = element_index(n_element_old)
              beta1 = information_2(1, index_ele)
              beta2 = information_2(2, index_ele)
              information_vector_1 = information_1(1:18, index_ele)
              information_vector_2 = information_2(1:8, index_ele)

              If (Abs(beta1 - out_object_beta) < SmallValue) Then
                piece_flag = 1
              Elseif (Abs(beta2 - out_object_beta) < SmallValue) Then
                piece_flag = 2
              End If

              Do n = 1, 4
                Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_old)
                ex_part = ex_part + Phi(HT(n, n_element_old), 1) * Ex_basis

                Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_old)
                ey_part = ey_part + Phi(HT(n, n_element_old), 1) * Ey_basis
              End Do

              Efield_old(1) = -ex_part*0.5
              Efield_old(2) = -ey_part*0.5

            Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
              Efield_old(1) = (W1 * efx(HT(1,n_element_old),1) + W2 * efx(HT(2,n_element_old),1) + &
                               W3 * efx(HT(3,n_element_old),1) + W4 * efx(HT(4,n_element_old),1)) * 0.5
              Efield_old(2) = (W1 * efy(HT(1,n_element_old),1) + W2 * efy(HT(2,n_element_old),1) + &
                               W3 * efy(HT(3,n_element_old),1) + W4 * efy(HT(4,n_element_old),1)) * 0.5
            End If
            !=========LY modification for 2D3V PIC model, 2022-5-27=========
            Efield_old(3) = (W1 * efz(HT(1,n_element_old),1) + W2 * efz(HT(2,n_element_old),1) + &
                             W3 * efz(HT(3,n_element_old),1) + W4 * efz(HT(4,n_element_old),1)) * 0.5
            !=========LY modification for 2D3V PIC model, 2022-5-27=========

            !Secondly, we need to calculate the Bfield of old element.
            Bfield_old(1) = (W1 * bfx(HT(1,n_element_old),1) + W2 * bfx(HT(2,n_element_old),1) + &
                             W3 * bfx(HT(3,n_element_old),1) + W4 * bfx(HT(4,n_element_old),1)) * 0.5
            Bfield_old(2) = (W1 * bfy(HT(1,n_element_old),1) + W2 * bfy(HT(2,n_element_old),1) + &
                             W3 * bfy(HT(3,n_element_old),1) + W4 * bfy(HT(4,n_element_old),1)) * 0.5
            Bfield_old(3) = (W1 * bfz(HT(1,n_element_old),1) + W2 * bfz(HT(2,n_element_old),1) + &
                             W3 * bfz(HT(3,n_element_old),1) + W4 * bfz(HT(4,n_element_old),1)) * 0.5
            !======The old element part======

            !======The bro element part======
            n_element_bro = HE(5, part_edge_count(3)) !The 0.5 particle distribute DG element

            hx_partition = HP(1, HT(2,n_element_bro)) - HP(1, HT(1,n_element_bro))
            hy_partition = HP(2, HT(4,n_element_bro)) - HP(2, HT(1,n_element_bro))

            dx = (PO%X - HP(1, HT(1,n_element_bro))) / hx_partition
            dy = (PO%Y - HP(2, HT(1,n_element_bro))) / hy_partition
            xcellmdx = 1.0 - dx
            ycellmdy = 1.0 - dy

            If (delta_global == 0) Then  !2D Cartesian coordinates.
              W1 = xcellmdx * ycellmdy
              W2 = dx       * ycellmdy
              W3 = dx       * dy
              W4 = xcellmdx * dy
            Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
              R1 = dymin + HP(2,HT(1,n_element_bro))
              R2 = dymin + HP(2,HT(4,n_element_bro))
              Rp = PO%Y
              den = R2*R2 - R1*R1

              !Here, we not consider the unequal weights for the moment, 2022-4-8.
              W1 = xcellmdx * (R2*R2-Rp*Rp) / den
              W2 = dx       * (R2*R2-Rp*Rp) / den
              W3 = dx       * (Rp*Rp-R1*R1) / den
              W4 = xcellmdx * (Rp*Rp-R1*R1) / den
            End If

            !Firstly, we need to calculate the Efield of bro element.
            If (element_index(n_element_bro) > 0) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
              ex_part = 0.0
              ey_part = 0.0
              trial_basis_type = 1
              vertices = HP(1:2, HT(1:4, n_element_bro))
              index_ele = element_index(n_element_bro)
              beta1 = information_2(1, index_ele)
              beta2 = information_2(2, index_ele)
              information_vector_1 = information_1(1:18, index_ele)
              information_vector_2 = information_2(1:8, index_ele)

              If (Abs(beta1 - out_object_beta) < SmallValue) Then
                piece_flag = 1
              Elseif (Abs(beta2 - out_object_beta) < SmallValue) Then
                piece_flag = 2
              End If

              Do n = 1, 4
                Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_bro)
                ex_part = ex_part + Phi(HT(n, n_element_bro), 1) * Ex_basis

                Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_bro)
                ey_part = ey_part + Phi(HT(n, n_element_bro), 1) * Ey_basis
              End Do

              Efield_bro(1) = -ex_part*0.5
              Efield_bro(2) = -ey_part*0.5

            Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
              Efield_bro(1) = (W1 * efx(HT(1,n_element_bro),1) + W2 * efx(HT(2,n_element_bro),1) + &
                               W3 * efx(HT(3,n_element_bro),1) + W4 * efx(HT(4,n_element_bro),1)) * 0.5
              Efield_bro(2) = (W1 * efy(HT(1,n_element_bro),1) + W2 * efy(HT(2,n_element_bro),1) + &
                               W3 * efy(HT(3,n_element_bro),1) + W4 * efy(HT(4,n_element_bro),1)) * 0.5
            End If
            !=========LY modification for 2D3V PIC model, 2022-5-27=========
            Efield_bro(3) = (W1 * efz(HT(1,n_element_bro),1) + W2 * efz(HT(2,n_element_bro),1) + &
                             W3 * efz(HT(3,n_element_bro),1) + W4 * efz(HT(4,n_element_bro),1)) * 0.5
            !=========LY modification for 2D3V PIC model, 2022-5-27=========

            !Secondly, we need to calculate the Bfield of bro element.
            Bfield_bro(1) = (W1 * bfx(HT(1,n_element_bro),1) + W2 * bfx(HT(2,n_element_bro),1) + &
                             W3 * bfx(HT(3,n_element_bro),1) + W4 * bfx(HT(4,n_element_bro),1)) * 0.5
            Bfield_bro(2) = (W1 * bfy(HT(1,n_element_bro),1) + W2 * bfy(HT(2,n_element_bro),1) + &
                             W3 * bfy(HT(3,n_element_bro),1) + W4 * bfy(HT(4,n_element_bro),1)) * 0.5
            Bfield_bro(3) = (W1 * bfz(HT(1,n_element_bro),1) + W2 * bfz(HT(2,n_element_bro),1) + &
                             W3 * bfz(HT(3,n_element_bro),1) + W4 * bfz(HT(4,n_element_bro),1)) * 0.5
            !======The bro element part======

            If (delta_global == 0) Then !2D Cartesian coordinates.
              xefield = Efield_old(1) + Efield_bro(1)
              yefield = Efield_old(2) + Efield_bro(2)
              zefield = Efield_old(3) + Efield_bro(3)

              A_x = xefield * qm(isp+1)
              A_y = yefield * qm(isp+1)
              A_z = zefield * qm(isp+1)
            Elseif (delta_global == 1) Then !2D axis-symmetric coordinates.
              !$ axisymmetric e field
              zefield = Efield_old(1) + Efield_bro(1)
              refield = Efield_old(2) + Efield_bro(2)
              tefield = Efield_old(3) + Efield_bro(3)
              !$ Convert axisymmetric efield to Cartesian efield
              Theta = PO%Z
              xefield = refield*DCOS(Theta) - tefield*DSIN(Theta)
              yefield = refield*DSIN(Theta) + tefield*DCOS(Theta)

              A_x = zefield * qm(isp+1)
              A_y = xefield * qm(isp+1)
              A_z = yefield * qm(isp+1)
            End If

            A_n_plus_1(1) = A_x
            A_n_plus_1(2) = A_y
            A_n_plus_1(3) = A_z

            !Thirdly, we need to calculate the Efield force and Bfield force of the particle.
            If (Bfiled_index) Then  !with magnetic field.
              If (delta_global == 0) Then  !2D Cartesian coordinates.
                xbfield = Bfield_old(1) + Bfield_bro(1)
                ybfield = Bfield_old(2) + Bfield_bro(2)
                zbfield = Bfield_old(3) + Bfield_bro(3)

                Omega(1) = xbfield * qm(isp+1) * 0.5 * TimeMove
                Omega(2) = ybfield * qm(isp+1) * 0.5 * TimeMove
                Omega(3) = zbfield * qm(isp+1) * 0.5 * TimeMove
              Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
                !$ axisymmetric b field
                zbfield = Bfield_old(1) + Bfield_bro(1)
                rbfield = Bfield_old(2) + Bfield_bro(2)
                tbfield = Bfield_old(3) + Bfield_bro(3)
                !$ Convert axisymmetric bfield to Cartesian bfield
                Theta = PO%Z
                xbfield = rbfield*DCOS(Theta) - tbfield*DSIN(Theta)
                ybfield = rbfield*DSIN(Theta) + tbfield*DCOS(Theta)

                Omega(1) = zbfield * qm(isp+1) * 0.5 * TimeMove
                Omega(2) = xbfield * qm(isp+1) * 0.5 * TimeMove
                Omega(3) = ybfield * qm(isp+1) * 0.5 * TimeMove
              End If

              OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
              TransB(1,1)=1.d0+Omega(1)*Omega(1)
              TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
              TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
              TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
              TransB(2,2)=1.d0+Omega(2)*Omega(2)
              TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
              TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
              TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
              TransB(3,3)=1.d0+Omega(3)*Omega(3)
              TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT

              SVelocity(1) = A_n_plus_1(1) * 0.5 * TimeMove
              SVelocity(2) = A_n_plus_1(2) * 0.5 * TimeMove
              SVelocity(3) = A_n_plus_1(3) * 0.5 * TimeMove

              deta_vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
              deta_vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
              deta_vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)

              PO%Vx = PO%Vx + deta_vx
              PO%Vy = PO%Vy + deta_vy
              PO%Vz = PO%Vz + deta_vz

            Else  !without magnetic field.
              deta_vx = 0.5*TimeMove*A_n_plus_1(1)
              deta_vy = 0.5*TimeMove*A_n_plus_1(2)
              deta_vz = 0.5*TimeMove*A_n_plus_1(3)

              PO%Vx = PO%Vx + deta_vx
              PO%Vy = PO%Vy + deta_vy
              PO%Vz = PO%Vz + deta_vz
            End If

            A_bar_n_minus_1(1) = PO%Ax
            A_bar_n_minus_1(2) = PO%Ay
            A_bar_n_minus_1(3) = PO%Az

            A_bar_n(1) = 0.5*(A_bar_n_minus_1(1) + A_n_plus_1(1))
            A_bar_n(2) = 0.5*(A_bar_n_minus_1(2) + A_n_plus_1(2))
            A_bar_n(3) = 0.5*(A_bar_n_minus_1(3) + A_n_plus_1(3))

            PO%Ax = A_bar_n(1)
            PO%Ay = A_bar_n(2)
            PO%Az = A_bar_n(3)

          End If

        Elseif (part_edge_count(1) == 4) Then
          !The particle locate the edge between CG(0.0) and DG(1.0). (condition 1)√
          !The particle locate the edge between CG(0.5) and DG(0.5). (condition 2)√

          edge_boundary = 0
          Do j = 2, 5
            If (HE(6, part_edge_count(j)) == 0) Then
              edge_boundary = edge_boundary + 1
            End If
          End Do

          If (edge_boundary == 2) Then
            n_element_old = HE(5, part_edge_count(2))   !There is to pretect array from outing of bound.
            traverse_flag = 0

          Else
            edge_bro_count = 0
            Do j = 2, 5
              If (HT(5, HE(5, part_edge_count(j)))>0 .AND. HT(5, HE(6, part_edge_count(j)))==0) Then
                edge_bro_element(1+edge_bro_count) = HE(6, part_edge_count(j))
                edge_bro_count = edge_bro_count + 1
              End If
            End Do

            If (edge_bro_count==2 .AND. edge_bro_element(1)==edge_bro_element(2)) Then
              !The particle locate the edge between CG(0.5) and DG(0.5). (condition 2: 2CG and 2DG--opposite, CG element is same)

              !======The old element part======
              n_element_old = HE(5, part_edge_count(2))   !The 0.5 particle distribute DG element

              hx_partition = HP(1, HT(2,n_element_old)) - HP(1, HT(1,n_element_old))
              hy_partition = HP(2, HT(4,n_element_old)) - HP(2, HT(1,n_element_old))

              dx = (PO%X - HP(1, HT(1,n_element_old))) / hx_partition
              dy = (PO%Y - HP(2, HT(1,n_element_old))) / hy_partition
              xcellmdx = 1.0 - dx
              ycellmdy = 1.0 - dy

              If (delta_global == 0) Then  !2D Cartesian coordinates.
                W1 = xcellmdx * ycellmdy
                W2 = dx       * ycellmdy
                W3 = dx       * dy
                W4 = xcellmdx * dy
              Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
                R1 = dymin + HP(2,HT(1,n_element_old))
                R2 = dymin + HP(2,HT(4,n_element_old))
                Rp = PO%Y
                den = R2*R2 - R1*R1

                !Here, we not consider the unequal weights for the moment, 2022-4-8.
                W1 = xcellmdx * (R2*R2-Rp*Rp) / den
                W2 = dx       * (R2*R2-Rp*Rp) / den
                W3 = dx       * (Rp*Rp-R1*R1) / den
                W4 = xcellmdx * (Rp*Rp-R1*R1) / den
              End If

              !Firstly, we need to calculate the Efield of old element.
              If (element_index(n_element_old) > 0) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
                ex_part = 0.0
                ey_part = 0.0
                trial_basis_type = 1
                vertices = HP(1:2, HT(1:4, n_element_old))
                index_ele = element_index(n_element_old)
                beta1 = information_2(1, index_ele)
                beta2 = information_2(2, index_ele)
                information_vector_1 = information_1(1:18, index_ele)
                information_vector_2 = information_2(1:8, index_ele)

                If (ABS(beta1 - out_object_beta) < SmallValue) Then
                  piece_flag = 1
                Elseif (ABS(beta2 - out_object_beta) < SmallValue) Then
                  piece_flag = 2
                End If

                Do n = 1, 4
                  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                  piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_old)
                  ex_part = ex_part + Phi(HT(n, n_element_old), 1) * Ex_basis

                  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                  piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_old)
                  ey_part = ey_part + Phi(HT(n, n_element_old), 1) * Ey_basis
                End Do

                Efield_old(1) = -ex_part*0.5
                Efield_old(2) = -ey_part*0.5

              Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
                Efield_old(1) = (W1 * efx(HT(1,n_element_old),1) + W2 * efx(HT(2,n_element_old),1) + &
                                 W3 * efx(HT(3,n_element_old),1) + W4 * efx(HT(4,n_element_old),1)) * 0.5
                Efield_old(2) = (W1 * efy(HT(1,n_element_old),1) + W2 * efy(HT(2,n_element_old),1) + &
                                 W3 * efy(HT(3,n_element_old),1) + W4 * efy(HT(4,n_element_old),1)) * 0.5
              End If
              !=========LY modification for 2D3V PIC model, 2022-5-27=========
              Efield_old(3) = (W1 * efz(HT(1,n_element_old),1) + W2 * efz(HT(2,n_element_old),1) + &
                               W3 * efz(HT(3,n_element_old),1) + W4 * efz(HT(4,n_element_old),1)) * 0.5
              !=========LY modification for 2D3V PIC model, 2022-5-27=========

              !Secondly, we need to calculate the Bfield of old element.
              Bfield_old(1) = (W1 * bfx(HT(1,n_element_old),1) + W2 * bfx(HT(2,n_element_old),1) + &
                               W3 * bfx(HT(3,n_element_old),1) + W4 * bfx(HT(4,n_element_old),1)) * 0.5
              Bfield_old(2) = (W1 * bfy(HT(1,n_element_old),1) + W2 * bfy(HT(2,n_element_old),1) + &
                               W3 * bfy(HT(3,n_element_old),1) + W4 * bfy(HT(4,n_element_old),1)) * 0.5
              Bfield_old(3) = (W1 * bfz(HT(1,n_element_old),1) + W2 * bfz(HT(2,n_element_old),1) + &
                               W3 * bfz(HT(3,n_element_old),1) + W4 * bfz(HT(4,n_element_old),1)) * 0.5
              !======The old element part======

              !======The bro element part======
              n_element_bro = edge_bro_element(1)   !The 0.5 particle distribute CG element

              hx_partition = HP(1, HT(2,n_element_bro)) - HP(1, HT(1,n_element_bro))
              hy_partition = HP(2, HT(4,n_element_bro)) - HP(2, HT(1,n_element_bro))

              dx = (PO%X - HP(1, HT(1,n_element_bro))) / hx_partition
              dy = (PO%Y - HP(2, HT(1,n_element_bro))) / hy_partition
              xcellmdx = 1.0 - dx
              ycellmdy = 1.0 - dy

              If (delta_global == 0) Then
                W1 = xcellmdx * ycellmdy
                W2 = dx       * ycellmdy
                W3 = dx       * dy
                W4 = xcellmdx * dy
              Elseif (delta_global == 1) Then
                R1 = dymin + HP(2,HT(1,n_element_bro))
                R2 = dymin + HP(2,HT(4,n_element_bro))
                Rp = PO%Y
                den = R2*R2 - R1*R1

                !Here, we not consider the unequal weights for the moment, 2022-4-8.
                W1 = xcellmdx * (R2*R2-Rp*Rp) / den
                W2 = dx       * (R2*R2-Rp*Rp) / den
                W3 = dx       * (Rp*Rp-R1*R1) / den
                W4 = xcellmdx * (Rp*Rp-R1*R1) / den
              End If

              !Firstly, we need to calculate the Efield of bro element.
              If (element_index(n_element_bro) > 0) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
                ex_part = 0.0
                ey_part = 0.0
                trial_basis_type = 1
                vertices = HP(1:2, HT(1:4, n_element_bro))
                index_ele = element_index(n_element_bro)
                beta1 = information_2(1, index_ele)
                beta2 = information_2(2, index_ele)
                information_vector_1 = information_1(1:18, index_ele)
                information_vector_2 = information_2(1:8, index_ele)

                If (Abs(beta1 - out_object_beta) < SmallValue) Then
                  piece_flag = 1
                Elseif (Abs(beta2 - out_object_beta) < SmallValue) Then
                  piece_flag = 2
                End If

                Do n = 1, 4
                  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                  piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_bro)
                  ex_part = ex_part + Phi(HT(n, n_element_bro), 1) * Ex_basis

                  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                  piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_bro)
                  ey_part = ey_part + Phi(HT(n, n_element_bro), 1) * Ey_basis
                End Do

                Efield_bro(1) = -ex_part*0.5
                Efield_bro(2) = -ey_part*0.5

              Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
                Efield_bro(1) = (W1 * efx(HT(1,n_element_bro),1) + W2 * efx(HT(2,n_element_bro),1) + &
                                 W3 * efx(HT(3,n_element_bro),1) + W4 * efx(HT(4,n_element_bro),1)) * 0.5
                Efield_bro(2) = (W1 * efy(HT(1,n_element_bro),1) + W2 * efy(HT(2,n_element_bro),1) + &
                                 W3 * efy(HT(3,n_element_bro),1) + W4 * efy(HT(4,n_element_bro),1)) * 0.5
              End If
              !=========LY modification for 2D3V PIC model, 2022-5-27=========
              Efield_bro(3) = (W1 * efz(HT(1,n_element_bro),1) + W2 * efz(HT(2,n_element_bro),1) + &
                               W3 * efz(HT(3,n_element_bro),1) + W4 * efz(HT(4,n_element_bro),1)) * 0.5
              !=========LY modification for 2D3V PIC model, 2022-5-27=========

              !Secondly, we need to calculate the Bfield of bro element.
              Bfield_bro(1) = (W1 * bfx(HT(1,n_element_bro),1) + W2 * bfx(HT(2,n_element_bro),1) + &
                               W3 * bfx(HT(3,n_element_bro),1) + W4 * bfx(HT(4,n_element_bro),1)) * 0.5
              Bfield_bro(2) = (W1 * bfy(HT(1,n_element_bro),1) + W2 * bfy(HT(2,n_element_bro),1) + &
                               W3 * bfy(HT(3,n_element_bro),1) + W4 * bfy(HT(4,n_element_bro),1)) * 0.5
              Bfield_bro(3) = (W1 * bfz(HT(1,n_element_bro),1) + W2 * bfz(HT(2,n_element_bro),1) + &
                               W3 * bfz(HT(3,n_element_bro),1) + W4 * bfz(HT(4,n_element_bro),1)) * 0.5
              !======The bro element part======

              If (delta_global == 0) Then !2D Cartesian coordinates.
                xefield = Efield_old(1) + Efield_bro(1)
                yefield = Efield_old(2) + Efield_bro(2)
                zefield = Efield_old(3) + Efield_bro(3)

                A_x = xefield * qm(isp+1)
                A_y = yefield * qm(isp+1)
                A_z = zefield * qm(isp+1)
              Elseif (delta_global == 1) Then !2D axis-symmetric coordinates.
                !$ axisymmetric e field
                zefield = Efield_old(1) + Efield_bro(1)
                refield = Efield_old(2) + Efield_bro(2)
                tefield = Efield_old(3) + Efield_bro(3)
                !$ Convert axisymmetric efield to Cartesian efield
                Theta = PO%Z
                xefield = refield*DCOS(Theta) - tefield*DSIN(Theta)
                yefield = refield*DSIN(Theta) + tefield*DCOS(Theta)

                A_x = zefield * qm(isp+1)
                A_y = xefield * qm(isp+1)
                A_z = yefield * qm(isp+1)
              End If

              A_n_plus_1(1) = A_x
              A_n_plus_1(2) = A_y
              A_n_plus_1(3) = A_z

              !Thirdly, we need to calculate the Efield force and Bfield force of the particle.
              If (Bfiled_index) Then  !with magnetic field.
                If (delta_global == 0) Then  !2D Cartesian coordinates.
                  xbfield = Bfield_old(1) + Bfield_bro(1)
                  ybfield = Bfield_old(2) + Bfield_bro(2)
                  zbfield = Bfield_old(3) + Bfield_bro(3)

                  Omega(1) = xbfield * qm(isp+1) * 0.5 * TimeMove
                  Omega(2) = ybfield * qm(isp+1) * 0.5 * TimeMove
                  Omega(3) = zbfield * qm(isp+1) * 0.5 * TimeMove
                Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
                  !$ axisymmetric b field
                  zbfield = Bfield_old(1) + Bfield_bro(1)
                  rbfield = Bfield_old(2) + Bfield_bro(2)
                  tbfield = Bfield_old(3) + Bfield_bro(3)
                  !$ Convert axisymmetric bfield to Cartesian bfield
                  Theta = PO%Z
                  xbfield = rbfield*DCOS(Theta) - tbfield*DSIN(Theta)
                  ybfield = rbfield*DSIN(Theta) + tbfield*DCOS(Theta)

                  Omega(1) = zbfield * qm(isp+1) * 0.5 * TimeMove
                  Omega(2) = xbfield * qm(isp+1) * 0.5 * TimeMove
                  Omega(3) = ybfield * qm(isp+1) * 0.5 * TimeMove
                End If

                OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
                TransB(1,1)=1.d0+Omega(1)*Omega(1)
                TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
                TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
                TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
                TransB(2,2)=1.d0+Omega(2)*Omega(2)
                TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
                TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
                TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
                TransB(3,3)=1.d0+Omega(3)*Omega(3)
                TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT

                SVelocity(1) = A_n_plus_1(1) * 0.5 * TimeMove
                SVelocity(2) = A_n_plus_1(2) * 0.5 * TimeMove
                SVelocity(3) = A_n_plus_1(3) * 0.5 * TimeMove

                deta_vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
                deta_vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
                deta_vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)

                PO%Vx = PO%Vx + deta_vx
                PO%Vy = PO%Vy + deta_vy
                PO%Vz = PO%Vz + deta_vz

              Else  !without magnetic field.
                deta_vx = 0.5*TimeMove*A_n_plus_1(1)
                deta_vy = 0.5*TimeMove*A_n_plus_1(2)
                deta_vz = 0.5*TimeMove*A_n_plus_1(3)

                PO%Vx = PO%Vx + deta_vx
                PO%Vy = PO%Vy + deta_vy
                PO%Vz = PO%Vz + deta_vz
              End If

              A_bar_n_minus_1(1) = PO%Ax
              A_bar_n_minus_1(2) = PO%Ay
              A_bar_n_minus_1(3) = PO%Az

              A_bar_n(1) = 0.5*(A_bar_n_minus_1(1) + A_n_plus_1(1))
              A_bar_n(2) = 0.5*(A_bar_n_minus_1(2) + A_n_plus_1(2))
              A_bar_n(3) = 0.5*(A_bar_n_minus_1(3) + A_n_plus_1(3))

              PO%Ax = A_bar_n(1)
              PO%Ay = A_bar_n(2)
              PO%Az = A_bar_n(3)

            Elseif (edge_bro_count==2 .AND. edge_bro_element(1)/=edge_bro_element(2)) Then
              !The particle locate the edge between CG(0.0) and DG(1.0). (condition 1: 2CG and 2DG--opposite, CG element is different)
              n_element_old = HE(5, part_edge_count(2))
              traverse_flag = 0

            Elseif (edge_bro_count==4) Then
              !The particle locate the edge between CG(0.0) and DG(1.0). (condition 1: 2CG and 2DG--diagonal, CG element is different)
              n_element_old = HE(5, part_edge_count(2))
              traverse_flag = 0

            End If
          End If

        Elseif (part_edge_count(1) == 6) Then
          !The particle locate the edge between CG(1.0) and DG(0.0).(1CG and 3DG, the particle locate the corner node)√
          !The particle locate the edge between DG(0.5) and DG(0.5).(2DG and 1DG, the particle locate the hand node)√

          edge_bro_count = 0
          edge_temp_count = 0
          edge_self_element = 0
          edge_self_element_min = Minval(HT(5, HE(5,part_edge_count(2:7))))

          Do j = 2, 7
            If (HT(5, HE(5,part_edge_count(j)))>0 .AND. HT(5, HE(6,part_edge_count(j)))==0) Then
              edge_bro_element(1+edge_bro_count) = HE(6, part_edge_count(j))
              edge_bro_count = edge_bro_count + 1
            Elseif (HT(5, HE(5,part_edge_count(j)))>0 .AND. HT(5, HE(6,part_edge_count(j)))>0) Then
              edge_temp_count = edge_temp_count + 1
              If (edge_self_element_min == HT(5, HE(5,part_edge_count(j)))) Then
                edge_self_element = part_edge_count(j)
              End If
            End If
          End Do

          If (edge_bro_count==2 .AND. edge_bro_element(1)==edge_bro_element(2)) Then
            !The particle locate the edge between CG(1.0) and DG(0.0).(1CG and 3DG, the particle locate the corner node)
            n_element_old = edge_bro_element(1)
            traverse_flag = 0

          Elseif (edge_temp_count == 6) Then
            !The particle locate the edge between DG(0.5) and DG(0.5).(2DG and 1DG, the particle locate the hand node)

            !======The old element part======
            n_element_old = HE(5, edge_self_element)

            hx_partition = HP(1, HT(2,n_element_old)) - HP(1, HT(1,n_element_old))
            hy_partition = HP(2, HT(4,n_element_old)) - HP(2, HT(1,n_element_old))

            dx = (PO%X - HP(1, HT(1,n_element_old))) / hx_partition
            dy = (PO%Y - HP(2, HT(1,n_element_old))) / hy_partition
            xcellmdx = 1.0 - dx
            ycellmdy = 1.0 - dy

            If (delta_global == 0) Then  !2D Cartesian coordinates.
              W1 = xcellmdx * ycellmdy
              W2 = dx       * ycellmdy
              W3 = dx       * dy
              W4 = xcellmdx * dy
            Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
              R1 = dymin + HP(2,HT(1,n_element_old))
              R2 = dymin + HP(2,HT(4,n_element_old))
              Rp = PO%Y
              den = R2*R2 - R1*R1

              !Here, we not consider the unequal weights for the moment, 2022-4-8.
              W1 = xcellmdx * (R2*R2-Rp*Rp) / den
              W2 = dx       * (R2*R2-Rp*Rp) / den
              W3 = dx       * (Rp*Rp-R1*R1) / den
              W4 = xcellmdx * (Rp*Rp-R1*R1) / den
            End If

            !Firstly, we need to calculate the Efield of old element.
            If (element_index(n_element_old) > 0) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
              ex_part = 0.0
              ey_part = 0.0
              trial_basis_type = 1
              vertices = HP(1:2, HT(1:4, n_element_old))
              index_ele = element_index(n_element_old)
              beta1 = information_2(1, index_ele)
              beta2 = information_2(2, index_ele)
              information_vector_1 = information_1(1:18, index_ele)
              information_vector_2 = information_2(1:8, index_ele)

              If (Abs(beta1 - out_object_beta) < SmallValue) Then
                piece_flag = 1
              Elseif (Abs(beta2 - out_object_beta) < SmallValue) Then
                piece_flag = 2
              End If

              Do n = 1, 4
                Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_old)
                ex_part = ex_part + Phi(HT(n, n_element_old), 1) * Ex_basis

                Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_old)
                ey_part = ey_part + Phi(HT(n, n_element_old), 1) * Ey_basis
              End Do

              Efield_old(1) = -ex_part*0.5
              Efield_old(2) = -ey_part*0.5

            Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
              Efield_old(1) = (W1 * efx(HT(1,n_element_old),1) + W2 * efx(HT(2,n_element_old),1) + &
                               W3 * efx(HT(3,n_element_old),1) + W4 * efx(HT(4,n_element_old),1)) * 0.5
              Efield_old(2) = (W1 * efy(HT(1,n_element_old),1) + W2 * efy(HT(2,n_element_old),1) + &
                               W3 * efy(HT(3,n_element_old),1) + W4 * efy(HT(4,n_element_old),1)) * 0.5
            End If
            !=========LY modification for 2D3V PIC model, 2022-5-27=========
            Efield_old(3) = (W1 * efz(HT(1,n_element_old),1) + W2 * efz(HT(2,n_element_old),1) + &
                             W3 * efz(HT(3,n_element_old),1) + W4 * efz(HT(4,n_element_old),1)) * 0.5
            !=========LY modification for 2D3V PIC model, 2022-5-27=========

            !Secondly, we need to calculate the Bfield of old element.
            Bfield_old(1) = (W1 * bfx(HT(1,n_element_old),1) + W2 * bfx(HT(2,n_element_old),1) + &
                             W3 * bfx(HT(3,n_element_old),1) + W4 * bfx(HT(4,n_element_old),1)) * 0.5
            Bfield_old(2) = (W1 * bfy(HT(1,n_element_old),1) + W2 * bfy(HT(2,n_element_old),1) + &
                             W3 * bfy(HT(3,n_element_old),1) + W4 * bfy(HT(4,n_element_old),1)) * 0.5
            Bfield_old(3) = (W1 * bfz(HT(1,n_element_old),1) + W2 * bfz(HT(2,n_element_old),1) + &
                             W3 * bfz(HT(3,n_element_old),1) + W4 * bfz(HT(4,n_element_old),1)) * 0.5
            !======The old element part======

            !======The bro element part======
            n_element_bro = HE(6, edge_self_element)   !The 0.5 particle distribute DG element

            hx_partition = HP(1, HT(2,n_element_bro)) - HP(1, HT(1,n_element_bro))
            hy_partition = HP(2, HT(4,n_element_bro)) - HP(2, HT(1,n_element_bro))

            dx = (PO%X - HP(1, HT(1,n_element_bro))) / hx_partition
            dy = (PO%Y - HP(2, HT(1,n_element_bro))) / hy_partition
            xcellmdx = 1.0 - dx
            ycellmdy = 1.0 - dy

            If (delta_global == 0) Then  !2D Cartesian coordinates.
              W1 = xcellmdx * ycellmdy
              W2 = dx       * ycellmdy
              W3 = dx       * dy
              W4 = xcellmdx * dy
            Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
              R1 = dymin + HP(2,HT(1,n_element_bro))
              R2 = dymin + HP(2,HT(4,n_element_bro))
              Rp = PO%Y
              den = R2*R2 - R1*R1

              !Here, we not consider the unequal weights for the moment, 2022-4-8.
              W1 = xcellmdx * (R2*R2-Rp*Rp) / den
              W2 = dx       * (R2*R2-Rp*Rp) / den
              W3 = dx       * (Rp*Rp-R1*R1) / den
              W4 = xcellmdx * (Rp*Rp-R1*R1) / den
            End If

            !Firstly, we need to calculate the Efield of bro element.
            If (element_index(n_element_bro) > 0) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
              ex_part = 0.0
              ey_part = 0.0
              trial_basis_type = 1
              vertices = HP(1:2, HT(1:4, n_element_bro))
              index_ele = element_index(n_element_bro)
              beta1 = information_2(1, index_ele)
              beta2 = information_2(2, index_ele)
              information_vector_1 = information_1(1:18, index_ele)
              information_vector_2 = information_2(1:8, index_ele)

              If (Abs(beta1 - out_object_beta) < SmallValue) Then
                piece_flag = 1
              Elseif (Abs(beta2 - out_object_beta) < SmallValue) Then
                piece_flag = 2
              End If

              Do n = 1, 4
                Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_bro)
                ex_part = ex_part + Phi(HT(n, n_element_bro), 1) * Ex_basis

                Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                                piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_bro)
                ey_part = ey_part + Phi(HT(n, n_element_bro), 1) * Ey_basis
              End Do

              Efield_bro(1) = -ex_part*0.5
              Efield_bro(2) = -ey_part*0.5

            Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
              Efield_bro(1) = (W1 * efx(HT(1,n_element_bro),1) + W2 * efx(HT(2,n_element_bro),1) + &
                               W3 * efx(HT(3,n_element_bro),1) + W4 * efx(HT(4,n_element_bro),1)) * 0.5
              Efield_bro(2) = (W1 * efy(HT(1,n_element_bro),1) + W2 * efy(HT(2,n_element_bro),1) + &
                               W3 * efy(HT(3,n_element_bro),1) + W4 * efy(HT(4,n_element_bro),1)) * 0.5
            End If
            !=========LY modification for 2D3V PIC model, 2022-5-27=========
            Efield_bro(3) = (W1 * efz(HT(1,n_element_bro),1) + W2 * efz(HT(2,n_element_bro),1) + &
                             W3 * efz(HT(3,n_element_bro),1) + W4 * efz(HT(4,n_element_bro),1)) * 0.5
            !=========LY modification for 2D3V PIC model, 2022-5-27=========

            !Secondly, we need to calculate the Bfield of bro element.
            Bfield_bro(1) = (W1 * bfx(HT(1,n_element_bro),1) + W2 * bfx(HT(2,n_element_bro),1) + &
                             W3 * bfx(HT(3,n_element_bro),1) + W4 * bfx(HT(4,n_element_bro),1)) * 0.5
            Bfield_bro(2) = (W1 * bfy(HT(1,n_element_bro),1) + W2 * bfy(HT(2,n_element_bro),1) + &
                             W3 * bfy(HT(3,n_element_bro),1) + W4 * bfy(HT(4,n_element_bro),1)) * 0.5
            Bfield_bro(3) = (W1 * bfz(HT(1,n_element_bro),1) + W2 * bfz(HT(2,n_element_bro),1) + &
                             W3 * bfz(HT(3,n_element_bro),1) + W4 * bfz(HT(4,n_element_bro),1)) * 0.5
            !======The bro element part======

            If (delta_global == 0) Then !2D Cartesian coordinates.
              xefield = Efield_old(1) + Efield_bro(1)
              yefield = Efield_old(2) + Efield_bro(2)
              zefield = Efield_old(3) + Efield_bro(3)

              A_x = xefield * qm(isp+1)
              A_y = yefield * qm(isp+1)
              A_z = zefield * qm(isp+1)
            Elseif (delta_global == 1) Then !2D axis-symmetric coordinates.
              !$ axisymmetric e field
              zefield = Efield_old(1) + Efield_bro(1)
              refield = Efield_old(2) + Efield_bro(2)
              tefield = Efield_old(3) + Efield_bro(3)
              !$ Convert axisymmetric efield to Cartesian efield
              Theta = PO%Z
              xefield = refield*DCOS(Theta) - tefield*DSIN(Theta)
              yefield = refield*DSIN(Theta) + tefield*DCOS(Theta)

              A_x = zefield * qm(isp+1)
              A_y = xefield * qm(isp+1)
              A_z = yefield * qm(isp+1)
            End If

            A_n_plus_1(1) = A_x
            A_n_plus_1(2) = A_y
            A_n_plus_1(3) = A_z

            !Thirdly, we need to calculate the Efield force and Bfield force of the particle.
            If (Bfiled_index) Then  !with magnetic field.
              If (delta_global == 0) Then  !2D Cartesian coordinates.
                xbfield = Bfield_old(1) + Bfield_bro(1)
                ybfield = Bfield_old(2) + Bfield_bro(2)
                zbfield = Bfield_old(3) + Bfield_bro(3)

                Omega(1) = xbfield * qm(isp+1) * 0.5 * TimeMove
                Omega(2) = ybfield * qm(isp+1) * 0.5 * TimeMove
                Omega(3) = zbfield * qm(isp+1) * 0.5 * TimeMove
              Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
                !$ axisymmetric b field
                zbfield = Bfield_old(1) + Bfield_bro(1)
                rbfield = Bfield_old(2) + Bfield_bro(2)
                tbfield = Bfield_old(3) + Bfield_bro(3)
                !$ Convert axisymmetric bfield to Cartesian bfield
                Theta = PO%Z
                xbfield = rbfield*DCOS(Theta) - tbfield*DSIN(Theta)
                ybfield = rbfield*DSIN(Theta) + tbfield*DCOS(Theta)

                Omega(1) = zbfield * qm(isp+1) * 0.5 * TimeMove
                Omega(2) = xbfield * qm(isp+1) * 0.5 * TimeMove
                Omega(3) = ybfield * qm(isp+1) * 0.5 * TimeMove
              End If

              OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
              TransB(1,1)=1.d0+Omega(1)*Omega(1)
              TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
              TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
              TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
              TransB(2,2)=1.d0+Omega(2)*Omega(2)
              TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
              TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
              TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
              TransB(3,3)=1.d0+Omega(3)*Omega(3)
              TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT

              SVelocity(1) = A_n_plus_1(1) * 0.5 * TimeMove
              SVelocity(2) = A_n_plus_1(2) * 0.5 * TimeMove
              SVelocity(3) = A_n_plus_1(3) * 0.5 * TimeMove

              deta_vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
              deta_vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
              deta_vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)

              PO%Vx = PO%Vx + deta_vx
              PO%Vy = PO%Vy + deta_vy
              PO%Vz = PO%Vz + deta_vz

            Else  !without magnetic field.
              deta_vx = 0.5*TimeMove*A_n_plus_1(1)
              deta_vy = 0.5*TimeMove*A_n_plus_1(2)
              deta_vz = 0.5*TimeMove*A_n_plus_1(3)

              PO%Vx = PO%Vx + deta_vx
              PO%Vy = PO%Vy + deta_vy
              PO%Vz = PO%Vz + deta_vz
            End If

            A_bar_n_minus_1(1) = PO%Ax
            A_bar_n_minus_1(2) = PO%Ay
            A_bar_n_minus_1(3) = PO%Az

            A_bar_n(1) = 0.5*(A_bar_n_minus_1(1) + A_n_plus_1(1))
            A_bar_n(2) = 0.5*(A_bar_n_minus_1(2) + A_n_plus_1(2))
            A_bar_n(3) = 0.5*(A_bar_n_minus_1(3) + A_n_plus_1(3))

            PO%Ax = A_bar_n(1)
            PO%Ay = A_bar_n(2)
            PO%Az = A_bar_n(3)

          End If

        Elseif (part_edge_count(1) == 8) Then
          !The particle locate the edge between DG(1.0) and DG(0.0).(4DG, the particle locate the corner node)√
          n_element_old = HE(5, part_edge_count(2))
          traverse_flag = 0

        Elseif (part_edge_count(1) == 0) Then
          !The particle locate on the edge of CG element without DG neighbor.
          n_element_old = CellMesh(n_element_old)%Finalindex
          traverse_flag = 0

        Else
          Write(6,*) 'part_edge_count(1) = ', part_edge_count(1)
          Write(6,*) 'PO%X = ', PO%X, 'PO%Y = ', PO%Y
          Stop
        End If
        !========================particle locate on CG and DG intersect edge: CG, DG coupling===========================
      End If
      
      If (traverse_flag == 0) Then
        !There denote the particle locate inside the element.
        hx_partition = HP(1, HT(2,n_element_old)) - HP(1, HT(1,n_element_old))
        hy_partition = HP(2, HT(4,n_element_old)) - HP(2, HT(1,n_element_old))

        dx = (PO%X - HP(1, HT(1,n_element_old))) / hx_partition
        dy = (PO%Y - HP(2, HT(1,n_element_old))) / hy_partition

        xcellmdx = 1.0 - dx
        ycellmdy = 1.0 - dy

        If (delta_global == 0) Then  !2D Cartesian coordinates.
          W1 = xcellmdx * ycellmdy
          W2 = dx       * ycellmdy
          W3 = dx       * dy
          W4 = xcellmdx * dy
        Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
          R1 = dymin + HP(2,HT(1,n_element_old))
          R2 = dymin + HP(2,HT(4,n_element_old))
          Rp = PO%Y
          den = R2*R2 - R1*R1

          !Here, we not consider the unequal weights for the moment, 2022-4-8.
          W1 = xcellmdx * (R2*R2-Rp*Rp) / den
          W2 = dx       * (R2*R2-Rp*Rp) / den
          W3 = dx       * (Rp*Rp-R1*R1) / den
          W4 = xcellmdx * (Rp*Rp-R1*R1) / den
        End If

        !Firstly, we need to calculate the Efield between interface element and non-interface element.
        If (element_index(n_element_old) >0 ) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
          ex_part = 0.0
          ey_part = 0.0
          trial_basis_type = 1
          vertices = HP(1:2, HT(1:4, n_element_old))
          index_ele = element_index(n_element_old)
          beta1 = information_2(1, index_ele)
          beta2 = information_2(2, index_ele)
          information_vector_1 = information_1(1:18, index_ele)
          information_vector_2 = information_2(1:8, index_ele)

          If (Abs(beta1 - out_object_beta) < SmallValue) Then
            piece_flag = 1
          Elseif (Abs(beta2 - out_object_beta) < SmallValue) Then
            piece_flag = 2
          End If

          Do n = 1, 4
            Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                            piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_old)
            ex_part = ex_part + Phi(HT(n, n_element_old), 1) * Ex_basis

            Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
                                            piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_old)
            ey_part = ey_part + Phi(HT(n, n_element_old), 1) * Ey_basis
          End Do

          Efield(1) = -ex_part
          Efield(2) = -ey_part

        Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
          Efield(1) = W1 * efx(HT(1,n_element_old),1) + W2 * efx(HT(2,n_element_old),1) + &
                      W3 * efx(HT(3,n_element_old),1) + W4 * efx(HT(4,n_element_old),1)
          Efield(2) = W1 * efy(HT(1,n_element_old),1) + W2 * efy(HT(2,n_element_old),1) + &
                      W3 * efy(HT(3,n_element_old),1) + W4 * efy(HT(4,n_element_old),1)
        End If
        !=========LY modification for 2D3V PIC model, 2022-5-27=========
        Efield(3) = W1 * efz(HT(1,n_element_old),1) + W2 * efz(HT(2,n_element_old),1) + &
                    W3 * efz(HT(3,n_element_old),1) + W4 * efz(HT(4,n_element_old),1)
        !=========LY modification for 2D3V PIC model, 2022-5-27=========

        If (delta_global == 0) Then !2D Cartesian coordinates.
          xefield = Efield(1)
          yefield = Efield(2)
          zefield = Efield(3)

          A_x = xefield * qm(isp+1)
          A_y = yefield * qm(isp+1)
          A_z = zefield * qm(isp+1)
        Elseif (delta_global == 1) Then !2D axis-symmetric coordinates.
          !$ axisymmetric e field
          zefield = Efield(1)
          refield = Efield(2)
          tefield = Efield(3)
          !$ Convert axisymmetric efield to Cartesian efield
          Theta = PO%Z
          xefield = refield*DCOS(Theta) - tefield*DSIN(Theta)
          yefield = refield*DSIN(Theta) + tefield*DCOS(Theta)

          A_x = zefield * qm(isp+1)
          A_y = xefield * qm(isp+1)
          A_z = yefield * qm(isp+1)
        End If

        A_n_plus_1(1) = A_x
        A_n_plus_1(2) = A_y
        A_n_plus_1(3) = A_z

        !Secondly, we need to calculate the Efield force and Bfield force.
        If (Bfiled_index) Then  !with magnetic field.
          If (delta_global == 0) Then  !2D Cartesian coordinates.

            xbfield = W1 * bfx(HT(1,n_element_old),1) + W2 * bfx(HT(2,n_element_old),1) + &
                      W3 * bfx(HT(3,n_element_old),1) + W4 * bfx(HT(4,n_element_old),1)
            ybfield = W1 * bfy(HT(1,n_element_old),1) + W2 * bfy(HT(2,n_element_old),1) + &
                      W3 * bfy(HT(3,n_element_old),1) + W4 * bfy(HT(4,n_element_old),1)
            zbfield = W1 * bfz(HT(1,n_element_old),1) + W2 * bfz(HT(2,n_element_old),1) + &
                      W3 * bfz(HT(3,n_element_old),1) + W4 * bfz(HT(4,n_element_old),1)

            Omega(1) = xbfield * qm(isp+1) * 0.5 * TimeMove
            Omega(2) = ybfield * qm(isp+1) * 0.5 * TimeMove
            Omega(3) = zbfield * qm(isp+1) * 0.5 * TimeMove

          Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
            !$ axisymmetric b field
            zbfield = W1 * bfx(HT(1,n_element_old),1) + W2 * bfx(HT(2,n_element_old),1) + &
                      W3 * bfx(HT(3,n_element_old),1) + W4 * bfx(HT(4,n_element_old),1)
            rbfield = W1 * bfy(HT(1,n_element_old),1) + W2 * bfy(HT(2,n_element_old),1) + &
                      W3 * bfy(HT(3,n_element_old),1) + W4 * bfy(HT(4,n_element_old),1)
            tbfield = W1 * bfz(HT(1,n_element_old),1) + W2 * bfz(HT(2,n_element_old),1) + &
                      W3 * bfz(HT(3,n_element_old),1) + W4 * bfz(HT(4,n_element_old),1)
            !$ Convert axisymmetric bfield to Cartesian bfield
            Theta = PO%Z
            xbfield = rbfield*DCOS(Theta) - tbfield*DSIN(Theta)
            ybfield = rbfield*DSIN(Theta) + tbfield*DCOS(Theta)
          
            Omega(1) = zbfield * qm(isp+1) * 0.5 * TimeMove
            Omega(2) = xbfield * qm(isp+1) * 0.5 * TimeMove
            Omega(3) = ybfield * qm(isp+1) * 0.5 * TimeMove
          End If

          OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
          TransB(1,1)=1.d0+Omega(1)*Omega(1)
          TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
          TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
          TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
          TransB(2,2)=1.d0+Omega(2)*Omega(2)
          TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
          TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
          TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
          TransB(3,3)=1.d0+Omega(3)*Omega(3)
          TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT

          SVelocity(1) = A_n_plus_1(1) * 0.5 * TimeMove
          SVelocity(2) = A_n_plus_1(2) * 0.5 * TimeMove
          SVelocity(3) = A_n_plus_1(3) * 0.5 * TimeMove

          deta_vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
          deta_vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
          deta_vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)

          PO%Vx = PO%Vx + deta_vx
          PO%Vy = PO%Vy + deta_vy
          PO%Vz = PO%Vz + deta_vz

        Else  !without magnetic field.

          deta_vx = 0.5 * TimeMove * A_n_plus_1(1)
          deta_vy = 0.5 * TimeMove * A_n_plus_1(2)
          deta_vz = 0.5 * TimeMove * A_n_plus_1(3)

          PO%Vx = PO%Vx + deta_vx
          PO%Vy = PO%Vy + deta_vy
          PO%Vz = PO%Vz + deta_vz
        End If

        A_bar_n_minus_1(1) = PO%Ax
        A_bar_n_minus_1(2) = PO%Ay
        A_bar_n_minus_1(3) = PO%Az
        
        A_bar_n(1) = 0.5 * (A_bar_n_minus_1(1) + A_n_plus_1(1))
        A_bar_n(2) = 0.5 * (A_bar_n_minus_1(2) + A_n_plus_1(2))
        A_bar_n(3) = 0.5 * (A_bar_n_minus_1(3) + A_n_plus_1(3))
        
        PO%Ax = A_bar_n(1)
        PO%Ay = A_bar_n(2)
        PO%Az = A_bar_n(3)
      End If
      !=========Part 1 : Positioning Particle and Solving velocity=========
      
      !xcellmdx = 1. -dx
      !ycellmdy = 1. -dy
      !
      !IF(delta_global == 0) THEN
      !  P1 = xcellmdx*ycellmdy
      !  P2 = dx*ycellmdy
      !  P3 = xcellmdx*dy
      !  P4 = dx*dy
      !ELSEIF(delta_global == 1)THEN
      !  R1=dymin + float(j - 1)*hx(2)
      !  R2=dymin + float(j)*hx(2)
      !  R = PO%Y
      !  den=R2*R2-R1*R1
      !
      !  P1 = xcellmdx*(R2*R2-R*R)/den
      !  P2 = dx*(R2*R2-R*R)/den
      !  P3 = xcellmdx*(R*R-R1*R1)/den
      !  P4 = dx*(R*R-R1*R1)/den
      !ENDIF
      !
      !IF (delta_global == 0) THEN
      !  xefield=efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4
      !  yefield=efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4
      !  zefield=efz(i,j)*P1 + efz(i+1,j)*P2 + efz(i,j+1)*P3 + efz(i+1,j+1)*P4
      !
      !  A_x = xefield * qm(isp+1)
      !  A_y = yefield * qm(isp+1)
      !  A_z = zefield * qm(isp+1)
      !ELSEIF (delta_global == 1) THEN
      !  !$ axisymmetric e field
      !  zefield = efx(i,j)*P1 + efx(i+1,j)*P2 + efx(i,j+1)*P3 + efx(i+1,j+1)*P4
      !  refield = efy(i,j)*P1 + efy(i+1,j)*P2 + efy(i,j+1)*P3 + efy(i+1,j+1)*P4
      !  tefield = efz(i,j)*P1 + efz(i+1,j)*P2 + efz(i,j+1)*P3 + efz(i+1,j+1)*P4
      !  !$ Convert axisymmetric efield to Cartesian efield
      !  Theta = PO%Z
      !  xefield = refield*DCOS(Theta) - tefield*DSIN(Theta)
      !  yefield = refield*DSIN(Theta) + tefield*DCOS(Theta)
      !
      !  A_x = zefield * qm(isp+1)
      !  A_y = xefield * qm(isp+1)
      !  A_z = yefield * qm(isp+1)
      !ENDIF
      !
      !A_n_plus_1(1) = A_x
      !A_n_plus_1(2) = A_y
      !A_n_plus_1(3) = A_z
      !
      !IF (Bfiled_index) THEN   !!!! 有磁场
      !
      !  IF (delta_global == 0) THEN
      !
      !    xbfield=bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4
      !    ybfield=bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4
      !    zbfield=bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4
      !
      !    Omega(1) = xbfield * qm(isp+1) * 0.5 * TimeMove
      !    Omega(2) = ybfield * qm(isp+1) * 0.5 * TimeMove
      !    Omega(3) = zbfield * qm(isp+1) * 0.5 * TimeMove
      !  ELSE IF (delta_global == 1) THEN
      !    !$ axisymmetric e field
      !    zbfield = bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4
      !    rbfield = bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4
      !    tbfield = bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4
      !    !$ Convert axisymmetric efield to Cartesian efield
      !    Theta = PO%Z
      !    xbfield = rbfield*DCOS(Theta) - tbfield*DSIN(Theta)
      !    ybfield = rbfield*DSIN(Theta) + tbfield*DCOS(Theta)
      !
      !    Omega(1) = zbfield * qm(isp+1) * 0.5 * TimeMove
      !    Omega(2) = xbfield * qm(isp+1) * 0.5 * TimeMove
      !    Omega(3) = ybfield * qm(isp+1) * 0.5 * TimeMove
      !  ENDIF
      !
      !  OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
      !  TransB(1,1)=1.d0+Omega(1)*Omega(1)
      !  TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
      !  TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
      !  TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
      !  TransB(2,2)=1.d0+Omega(2)*Omega(2)
      !  TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
      !  TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
      !  TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
      !  TransB(3,3)=1.d0+Omega(3)*Omega(3)
      !  TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT
      !
      !  SVelocity(1) = A_n_plus_1(1) * 0.5 * TimeMove
      !  SVelocity(2) = A_n_plus_1(2) * 0.5 * TimeMove
      !  SVelocity(3) = A_n_plus_1(3) * 0.5 * TimeMove
      !
      !  deta_vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
      !  deta_vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
      !  deta_vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)
      !
      !  PO%Vx = PO%Vx + deta_vx
      !  PO%Vy = PO%Vy + deta_vy
      !  PO%Vz = PO%Vz + deta_vz
      !
      !ELSE  !!!! 无磁场
      !
      !  deta_vx = 0.5*TimeMove*A_n_plus_1(1)
      !  deta_vy = 0.5*TimeMove*A_n_plus_1(2)
      !  deta_vz = 0.5*TimeMove*A_n_plus_1(3)
      !
      !  PO%Vx = PO%Vx + deta_vx
      !  PO%Vy = PO%Vy + deta_vy
      !  PO%Vz = PO%Vz + deta_vz
      !
      !END IF
      !
      !A_bar_n_minus_1(1) = PO%Ax
      !A_bar_n_minus_1(2) = PO%Ay
      !A_bar_n_minus_1(3) = PO%Az
      !
      !A_bar_n(1) = 0.5*(A_bar_n_minus_1(1) + A_n_plus_1(1))
      !A_bar_n(2) = 0.5*(A_bar_n_minus_1(2) + A_n_plus_1(2))
      !A_bar_n(3) = 0.5*(A_bar_n_minus_1(3) + A_n_plus_1(3))
      !
      !PO%Ax = A_bar_n(1)
      !PO%Ay = A_bar_n(2)
      !PO%Az = A_bar_n(3)

      IvelFlag = 0
    End If
    
    detaV(1) = deta_vx
    detaV(2) = deta_vy
    detaV(3) = deta_vz

    IF(delta == 0)THEN
      PO%X = PO%X + TimeMove*deta_vx
      PO%Y = PO%Y + TimeMove*deta_vy
      PO%Z = PO%Z + TimeMove*deta_vz
    ELSEIF(delta == 1)THEN
      X=PO%Y*DCOS(PO%Z)
      Y=PO%Y*DSIN(PO%Z)
      Z=PO%X

      Z = Z + TimeMove*deta_vx
      X = X + TimeMove*deta_vy
      IF(X == 0.) THEN
        X=hx(2)*1.0E-5
      ENDIF
      Y = Y + TimeMove*deta_vz

      PO%X = Z
      PO%Y = DSQRT(X*X+Y*Y)
      PO%Z = DATAN(Y/X)
      IF (X <= 0.0) THEN
        PO%Z=PO%Z+PI_LY
      ENDIF
    ENDIF

    IposFlag = 0
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    Deallocate(HE_particle)
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    
End Subroutine PostPushOne


