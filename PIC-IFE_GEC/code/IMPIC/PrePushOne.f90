SUBROUTINE PrePushOne(PO,isp,TimeMove, IvelFlag, IposFlag, delta, OriPosi)
  
    USE Domain_2D
    USE Particle_2D, Only: qm
    USE Field_2D
    USE IMPIC_Data_2D
    Use ModuleParticleOne
    !============LY add for SIDG-PIC couple, 2022-7-25============
    Use IFE_Data
    !============LY add for SIDG-PIC couple, 2022-7-25============
    IMPLICIT NONE
    
    !Class(ParticleOne), intent(inout) :: PO
    Type(ParticleOne), intent(inout) :: PO
    Integer(4), intent(in) :: isp
    Real(8), intent(in) :: TimeMove
    Integer(4), intent(inout) :: IvelFlag, IposFlag
    Integer(4), intent(in) :: delta
    Real(8), intent(out) :: OriPosi(3)
    
    Real(8) :: xp, yp
    Real(8) :: dx, dy, xcellmdx, ycellmdy
    Real(8) :: Rp, R1, R2, den
    Real(8) :: P1, P2, P3, P4
    Integer(4) :: i, j
    
    Real(8) :: f
    Real(8) :: SVelocity(1:3)=0.d0
    Real(8) :: deta_vx, deta_vy, deta_vz
    REAL(8) :: X, Y, Z, Theta
    REAL(8) :: zefield, refield, tefield
    REAL(8) :: xefield, yefield
    REAL(8) :: zbfield, rbfield, tbfield
    REAL(8) :: xbfield, ybfield
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    Integer :: traverse_flag
    Integer :: n_element_old, n_element_bro
    Real(8) :: part_ele_left, part_ele_right, part_ele_bottom, part_ele_top
    Real(8) :: bottom_top_centre, left_right_centre, delta_1, delta_2
    Integer :: edge_count
    Integer :: edge_bro_element(4), part_edge_count(9)
    Real(8) :: hx_partition, hy_partition
    Real(8) :: W1, W2, W3, W4
    Real(8) :: Bfield_old(3), Bfield_bro(3)
    Integer :: edge_boundary, edge_bro_count, edge_temp_count, edge_self_element, edge_self_element_min
    Integer :: number_ele_column
    Integer :: n_element_bro_initial
    Real(8), Parameter :: PI_LY	= 3.14159265358979D0
    Real(8), Dimension(:,:), Allocatable :: HE_particle
    Integer, Dimension(:), Allocatable :: EdgeArray
    
    Allocate(HE_particle(5, Size(HE,2)))
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    

    OriPosi = (/PO%X,PO%Y,PO%Z/)

    IF(IvelFlag == 1) THEN
      
      !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
      !==========Part 0 : Initialization==========
      edge_bro_element(1:4) = 0
      number_ele_column = ny - 1
      If (Allocated(EdgeArray)) Then
        Deallocate(EdgeArray)
      End If
      !==========Part 0 : Initialization==========
      
      !==========Part 1 : Positioning Particle and Moving==========
      IF (Bfiled_index) THEN   !with magnetic field.
        
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
                Write(6,*) 'The delta2 = PO%X - left_right_centre is 0, but traverse_flag = 0, ERROR--PrePsuh'
                Stop
              End If
            Elseif (delta_1 < -SmallValue) Then   !Particle locate inside left-bottom child element or right-bottom child element.
              If (delta_2 > SmallValue) Then   !Particle locate inside right-bottom child element.
                n_element_old = CellMesh(n_element_old)%Child(3)
              Elseif (delta_2 < -SmallValue) Then   !Particle locate inside left-bottom child element.
                n_element_old = CellMesh(n_element_old)%Child(1)
              Elseif (ABS(delta_2) < SmallValue) Then
                Write(6,*) 'The delta2 = PO%X - left_right_centre is 0, but traverse_flag = 0, ERROR--PrePsuh'
                Stop
              End If
            Elseif (ABS(delta_1) < SmallValue) Then
              Write(6,*) 'The delta1 = PO%Y - bottom_top_centre is 0, but traverse_flag = 0, ERROR--PrePsuh'
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

                !Firstly, we need to calculate the Bfield of old element.
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

                !Secondly, we need to calculate the Bfield of bro element.
                Bfield_bro(1) = (W1 * bfx(HT(1,n_element_bro),1) + W2 * bfx(HT(2,n_element_bro),1) + &
                                 W3 * bfx(HT(3,n_element_bro),1) + W4 * bfx(HT(4,n_element_bro),1)) * 0.5
                Bfield_bro(2) = (W1 * bfy(HT(1,n_element_bro),1) + W2 * bfy(HT(2,n_element_bro),1) + &
                                 W3 * bfy(HT(3,n_element_bro),1) + W4 * bfy(HT(4,n_element_bro),1)) * 0.5
                Bfield_bro(3) = (W1 * bfz(HT(1,n_element_bro),1) + W2 * bfz(HT(2,n_element_bro),1) + &
                                 W3 * bfz(HT(3,n_element_bro),1) + W4 * bfz(HT(4,n_element_bro),1)) * 0.5
                !======The bro element part======

                !Thirdly, we need to calculate the Bfield force of the particle.
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

                SVelocity(1) = PO%Vx + 0.5 * PO%Ax * TimeMove + &
                               PO%Vy * Omega(3) - PO%Vz * Omega(2)
                SVelocity(2) = PO%Vy + 0.5 * PO%Ay * TimeMove + &
                               PO%Vz * Omega(1) - PO%Vx * Omega(3)
                SVelocity(3) = PO%Vz + 0.5 * PO%Az * TimeMove + &
                               PO%Vx * Omega(2) - PO%Vy * Omega(1)

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

                PO%Vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
                PO%Vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
                PO%Vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)

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

              !Firstly, we need to calculate the Bfield of old element.
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

              !Secondly, we need to calculate the Bfield of bro element.
              Bfield_bro(1) = (W1 * bfx(HT(1,n_element_bro),1) + W2 * bfx(HT(2,n_element_bro),1) + &
                               W3 * bfx(HT(3,n_element_bro),1) + W4 * bfx(HT(4,n_element_bro),1)) * 0.5
              Bfield_bro(2) = (W1 * bfy(HT(1,n_element_bro),1) + W2 * bfy(HT(2,n_element_bro),1) + &
                               W3 * bfy(HT(3,n_element_bro),1) + W4 * bfy(HT(4,n_element_bro),1)) * 0.5
              Bfield_bro(3) = (W1 * bfz(HT(1,n_element_bro),1) + W2 * bfz(HT(2,n_element_bro),1) + &
                               W3 * bfz(HT(3,n_element_bro),1) + W4 * bfz(HT(4,n_element_bro),1)) * 0.5
              !======The bro element part======

              !Thirdly, we need to calculate the Bfield force of the particle.
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

              SVelocity(1) = PO%Vx + 0.5 * PO%Ax * TimeMove + &
                             PO%Vy * Omega(3) - PO%Vz * Omega(2)
              SVelocity(2) = PO%Vy + 0.5 * PO%Ay * TimeMove + &
                             PO%Vz * Omega(1) - PO%Vx * Omega(3)
              SVelocity(3) = PO%Vz + 0.5 * PO%Az * TimeMove + &
                             PO%Vx * Omega(2) - PO%Vy * Omega(1)

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

              PO%Vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
              PO%Vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
              PO%Vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3) 
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

                !Firstly, we need to calculate the Bfield of old element.
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

                !Secondly, we need to calculate the Bfield of bro element.
                Bfield_bro(1) = (W1 * bfx(HT(1,n_element_bro),1) + W2 * bfx(HT(2,n_element_bro),1) + &
                                 W3 * bfx(HT(3,n_element_bro),1) + W4 * bfx(HT(4,n_element_bro),1)) * 0.5
                Bfield_bro(2) = (W1 * bfy(HT(1,n_element_bro),1) + W2 * bfy(HT(2,n_element_bro),1) + &
                                 W3 * bfy(HT(3,n_element_bro),1) + W4 * bfy(HT(4,n_element_bro),1)) * 0.5
                Bfield_bro(3) = (W1 * bfz(HT(1,n_element_bro),1) + W2 * bfz(HT(2,n_element_bro),1) + &
                                 W3 * bfz(HT(3,n_element_bro),1) + W4 * bfz(HT(4,n_element_bro),1)) * 0.5
                !======The bro element part======

                !Thirdly, we need to calculate the Efield force and Bfield force of the particle.
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

                SVelocity(1) = PO%Vx + 0.5 * PO%Ax * TimeMove + &
                               PO%Vy * Omega(3) - PO%Vz * Omega(2)
                SVelocity(2) = PO%Vy + 0.5 * PO%Ay * TimeMove + &
                               PO%Vz * Omega(1) - PO%Vx * Omega(3)
                SVelocity(3) = PO%Vz + 0.5 * PO%Az * TimeMove + &
                               PO%Vx * Omega(2) - PO%Vy * Omega(1)

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

                PO%Vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
                PO%Vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
                PO%Vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)

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
            edge_self_element_min = MINVAL(HT(5, HE(5,part_edge_count(2:7))))

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

              !Firstly, we need to calculate the Bfield of old element.
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

              !Secondly, we need to calculate the Bfield of bro element.
              Bfield_bro(1) = (W1 * bfx(HT(1,n_element_bro),1) + W2 * bfx(HT(2,n_element_bro),1) + &
                               W3 * bfx(HT(3,n_element_bro),1) + W4 * bfx(HT(4,n_element_bro),1)) * 0.5
              Bfield_bro(2) = (W1 * bfy(HT(1,n_element_bro),1) + W2 * bfy(HT(2,n_element_bro),1) + &
                               W3 * bfy(HT(3,n_element_bro),1) + W4 * bfy(HT(4,n_element_bro),1)) * 0.5
              Bfield_bro(3) = (W1 * bfz(HT(1,n_element_bro),1) + W2 * bfz(HT(2,n_element_bro),1) + &
                               W3 * bfz(HT(3,n_element_bro),1) + W4 * bfz(HT(4,n_element_bro),1)) * 0.5
              !======The bro element part======

              !Thirdly, we need to calculate the Efield force and Bfield force of the particle.
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

              SVelocity(1) = PO%Vx + 0.5 * PO%Ax * TimeMove + &
                             PO%Vy * Omega(3) - PO%Vz * Omega(2)
              SVelocity(2) = PO%Vy + 0.5 * PO%Ay * TimeMove + &
                             PO%Vz * Omega(1) - PO%Vx * Omega(3)
              SVelocity(3) = PO%Vz + 0.5 * PO%Az * TimeMove + &
                             PO%Vx * Omega(2) - PO%Vy * Omega(1)

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

              PO%Vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
              PO%Vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
              PO%Vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)
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

          !Firstly, we need to calculate the Bfield.
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

          SVelocity(1) = PO%Vx + 0.5 * PO%Ax * TimeMove + &
                         PO%Vy * Omega(3) - PO%Vz * Omega(2)
          SVelocity(2) = PO%Vy + 0.5 * PO%Ay * TimeMove + &
                         PO%Vz * Omega(1) - PO%Vx * Omega(3)
          SVelocity(3) = PO%Vz + 0.5 * PO%Az * TimeMove + &
                         PO%Vx * Omega(2) - PO%Vy * Omega(1)

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

          PO%Vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
          PO%Vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
          PO%Vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)

        End If
        
      Else  !without magnetic field.

        PO%Vx = PO%Vx + 0.5*TimeMove*PO%Ax
        PO%Vy = PO%Vy + 0.5*TimeMove*PO%Ay
        PO%Vz = PO%Vz + 0.5*TimeMove*PO%Az

      End if
      !==========Part 1 : Positioning Particle and Moving==========
      !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
      
      !IF (Bfiled_index) THEN   !!!! 有磁场
      !  xp = (PO%X - Vert_o(1))*hxi(1)
      !  i = xp        !!! 粒子在第i个单元里
      !  dx = xp-i
      !
      !  yp = (PO%Y - Vert_o(2))*hxi(2)
      !  j = yp        !!! 粒子在第j个单元里
      !  dy = yp-j
      !
      !  xcellmdx = 1. -dx
      !  ycellmdy = 1. -dy
      !
      !  IF(delta == 0) THEN
      !    P1 = xcellmdx*ycellmdy
      !    P2 = dx      *ycellmdy
      !    P3 = xcellmdx*dy
      !    P4 = dx      *dy
      !  ELSEIF(delta == 1)THEN
      !    R1=dymin + float(j - 1)*hx(2)
      !    R2=dymin + float(j)*hx(2)
      !    Rp = PO%Y
      !    den=R2*R2-R1*R1
      !
      !    P1 = xcellmdx*(R2*R2-Rp*Rp)/den
      !    P2 = dx      *(R2*R2-Rp*Rp)/den
      !    P3 = xcellmdx*(Rp*Rp-R1*R1)/den
      !    P4 = dx      *(Rp*Rp-R1*R1)/den
      !  ENDIF
      !
      !  IF (delta == 0) THEN
      !
      !    xbfield=bfx(i,j)*P1 + bfx(i+1,j)*P2 + bfx(i,j+1)*P3 + bfx(i+1,j+1)*P4
      !    ybfield=bfy(i,j)*P1 + bfy(i+1,j)*P2 + bfy(i,j+1)*P3 + bfy(i+1,j+1)*P4
      !    zbfield=bfz(i,j)*P1 + bfz(i+1,j)*P2 + bfz(i,j+1)*P3 + bfz(i+1,j+1)*P4
      !
      !    Omega(1) = xbfield * qm(isp+1) * 0.5 * TimeMove
      !    Omega(2) = ybfield * qm(isp+1) * 0.5 * TimeMove
      !    Omega(3) = zbfield * qm(isp+1) * 0.5 * TimeMove
      !
      !  ELSE IF (delta == 1) THEN
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
      !  SVelocity(1) = PO%Vx + 0.5 * PO%Ax * TimeMove + &
      !    PO%Vy * Omega(3) - PO%Vz * Omega(2)
      !  SVelocity(2) = PO%Vy + 0.5 * PO%Ay * TimeMove + &
      !    PO%Vz * Omega(1) - PO%Vx * Omega(3)
      !  SVelocity(3) = PO%Vz + 0.5 * PO%Az * TimeMove + &
      !    PO%Vx * Omega(2) - PO%Vy * Omega(1)
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
      !  PO%Vx = TransB(1,1)*SVelocity(1)+TransB(1,2)*SVelocity(2)+TransB(1,3)*SVelocity(3)
      !  PO%Vy = TransB(2,1)*SVelocity(1)+TransB(2,2)*SVelocity(2)+TransB(2,3)*SVelocity(3)
      !  PO%Vz = TransB(3,1)*SVelocity(1)+TransB(3,2)*SVelocity(2)+TransB(3,3)*SVelocity(3)
      !
      !Else  !!!! 无磁场
      !
      !  PO%Vx = PO%Vx + 0.5*TimeMove*PO%Ax
      !  PO%Vy = PO%Vy + 0.5*TimeMove*PO%Ay
      !  PO%Vz = PO%Vz + 0.5*TimeMove*PO%Az
      !
      !End if
      IvelFlag = 0
    End if

    If (delta == 0) Then
      PO%X = PO%X + TimeMove*PO%Vx
      PO%Y = PO%Y + TimeMove*PO%Vy
    Elseif (delta == 1) Then
      !> transform positions to cartesian coordinates
      X=PO%Y * DCOS(PO%Z)
      Y=PO%Y * DSIN(PO%Z)
      Z=PO%X
      !> update cartesian positions
      X=X+TimeMove * PO%Vy
      IF(X == 0.) THEN
        X=hx(2)*1.0E-5
      ENDIF
      Y=Y+TimeMove * PO%Vz

      PO%X = PO%X + TimeMove*PO%Vx
      !> update polar positions
      PO%Y = DSQRT(X*X+Y*Y)
      PO%Z = DATAN(Y/X)
      !> place particle in proper quadrant
      If (X <= 0.0) Then
        PO%Z=PO%Z+PI_LY
      Endif
    EndIf
    IposFlag = 0
    
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    Deallocate(HE_particle)
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
End Subroutine PrePushOne