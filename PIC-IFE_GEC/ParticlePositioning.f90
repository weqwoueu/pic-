Subroutine ParticlePositioning(part_x, part_y, n_element_old, traverse_flag, EdgeArray)
Use IFE_Data
Use Domain_2D
IMPLICIT NONE

!=========LY modification for Multi-Layer-Grid, 2022-7-25=========

!========= Input and Output ================
Integer, Intent(InOut) :: n_element_old
!================================

!========= Input ================
Real(8), Intent(In) :: part_x, part_y
!================================

!========= Output ================
Integer, Intent(Out) :: traverse_flag
Integer, Allocatable, Intent(Out) :: EdgeArray(:)
!================================


Integer :: i, j, k, l
Integer :: n_element_bro
Integer :: n_element_bro_initial, number_ele_column
Real(8) :: part_ele_left, part_ele_right, part_ele_bottom, part_ele_top
Real(8) :: left_right_centre, bottom_top_centre, delta_1, delta_2
Integer :: part_edge_count(9)
Integer :: edge_count
Integer :: num
Real(8) :: hx_partition, hy_partition
Real(8) :: W1, W2, W3, W4
Integer :: edge_boundary
Integer :: edge_bro_count, edge_temp_count, edge_self_element, edge_self_element_min
Integer :: edge_bro_element(4)
Integer :: temp_count
Real(8), Dimension(:,:), Allocatable :: HE_particle
!=========LY modification for Multi-Layer-Grid, 2022-7-25=========

number_ele_column = ny - 1

Do While (CellMesh(n_element_old)%isSplitted == 1)
    !In SIDG/IDG method, if we don't refine any element, that is use initial grid,
    !so we set the 'isSplitted' flag of every element is 0 in Setup_IFE_Mesh_2D.f90
    !The Do-While circle only use grid refinement. There denote the initial element will be refined.

    part_ele_left = part_x - CellMesh(n_element_old)%Boundary(1)    !Distance between particle and element's left boundary.
    part_ele_right = part_x - CellMesh(n_element_old)%Boundary(2)   !Distance between particle and element's right boundary.
    part_ele_bottom = part_y - CellMesh(n_element_old)%Boundary(3)  !Distance between particle and element's bottom boundary.
    part_ele_top = part_y - CellMesh(n_element_old)%Boundary(4)     !Distance between particle and element's top boundary.

    left_right_centre = (CellMesh(n_element_old)%Boundary(1)+CellMesh(n_element_old)%Boundary(2)) / 2.0
    bottom_top_centre = (CellMesh(n_element_old)%Boundary(3)+CellMesh(n_element_old)%Boundary(4)) / 2.0
    delta_1 = part_y - bottom_top_centre  !Distance between particle and element's centre of bottom and top.
    delta_2 = part_x - left_right_centre  !Distance between particle and element's centre of left and right.

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
            Exit
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
            Exit
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
            Exit
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
            Exit
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
            Exit
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
        Exit
    End If

    Else

    If (delta_1 > SmallValue) Then !Left-Top or Right-Top child element
        If (delta_2 > SmallValue) Then !Right-Top child element
        n_element_old = CellMesh(n_element_old)%Child(4)
        Elseif (delta_2 < -SmallValue) Then  !Left-Top child element
        n_element_old = CellMesh(n_element_old)%Child(2)
        Elseif (ABS(delta_2) < SmallValue) Then
        Write(6,*) 'The delta2 = part_x - left_right_centre is 0, but traverse_flag = 0, ERROR--GetParChg_2D'
        Stop
        End If
    Elseif (delta_1 < -SmallValue) Then  !Left-Bottom or Right-Bottom child element
        If (delta_2 > SmallValue) Then !Right-Bottom child element
        n_element_old = CellMesh(n_element_old)%Child(3)
        Elseif (delta_2 < -SmallValue) Then  !Left-Bottom child element
        n_element_old = CellMesh(n_element_old)%Child(1)
        Elseif (ABS(delta_2) < SmallValue) Then
        Write(6,*) 'The delta2 = part_x - left_right_centre is 0, but traverse_flag = 0, ERROR--GetParChg_2D'
        Stop
        End If
    Elseif (ABS(delta_1) < SmallValue) Then
        Write(6,*) 'The delta1 = part_y - bottom_top_centre is 0, but traverse_flag = 0, ERROR--GetParChg_2D'
        Stop
    End If

    End If
End Do
    
If (CellMesh(n_element_old)%isSplitted==0 .AND. CellMesh(n_element_old)%Finalindex/=0) Then
    !There denote the element is no-refined element.

    part_ele_left = part_x - CellMesh(n_element_old)%Boundary(1)    !Distance between particle and element's left boundary.
    part_ele_right = part_x - CellMesh(n_element_old)%Boundary(2)   !Distance between particle and element's right boundary.
    part_ele_bottom = part_y - CellMesh(n_element_old)%Boundary(3)  !Distance between particle and element's bottom boundary.
    part_ele_top = part_y - CellMesh(n_element_old)%Boundary(4)     !Distance between particle and element's top boundary.

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


End Subroutine ParticlePositioning