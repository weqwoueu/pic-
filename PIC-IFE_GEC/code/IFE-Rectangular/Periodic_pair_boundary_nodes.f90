Subroutine Periodic_pair_boundary_nodes(nnx, nny, bc_index, p_basic, bc_point_1, bc_point_2, PairNodes)
  
  Implicit None
  
  Integer :: nnx, nny
  Integer, Dimension(:), Pointer :: bc_index, match_flag
  Real(8), Dimension(:,:), Pointer :: p_basic
  Real(8), Dimension(:,:), Pointer :: bc_point_1, bc_point_2
  Integer, Dimension(:,:), Pointer :: PairNodesTemp
  Integer, Dimension(:,:), Pointer :: PairNodes
  Integer, Dimension(:,:), Pointer ::	check_periodic_boundary
  Integer, Dimension(:,:), Pointer ::	ALQ_index
  Real(8), Parameter :: SmallValue = 1.0D-8
  Integer :: x_min, y_min, x_max, y_max, num_unknows, i, j, m, n, n_boundary
  Real(8) :: x1_value_min, x1_value_max, y1_value_min, y1_value_max, x2_value_min, x2_value_max, y2_value_min, y2_value_max
  Integer :: num_periodic_boundary, num, num_corresponding_boundary
  
  num_unknows = SIZE(p_basic,2) !Number of Basis nodes
  x_min = MINVAL(p_basic(1,:))  !left boundary
  x_max = MAXVAL(p_basic(1,:))  !right boundary
  y_min = MINVAL(p_basic(2,:))  !bottom boundary
  y_max = MAXVAL(p_basic(2,:))  !top boundary
  
  n_boundary = SIZE(bc_index)   !Number of boundary
  
  num_periodic_boundary = 0     !Number of periodic boundary
  
  Allocate(ALQ_index(5,n_boundary))
  ALQ_index = -200

  Do i = 1,n_boundary
    If (bc_index(i) == -1) Then
      num_periodic_boundary = num_periodic_boundary+1
      ALQ_index(1,i) = i
      !ALQ_index(2:5,num):the num line's begin and end point
      ALQ_index(2,i) = bc_point_1(1,i)  !Boundary start point x
      ALQ_index(3,i) = bc_point_1(2,i)  !Boundary start point y
      ALQ_index(4,i) = bc_point_2(1,i)  !Boundary end point x
      ALQ_index(5,i) = bc_point_2(2,i)  !Boundary end point y
    End If
  End Do
  
  !===================Check the boundary conditions===========================
  If (MOD(num_periodic_boundary,2) /= 0) Then
    Write(*,*) "Please check the periodic boudary conditions! "
    Write(*,*) "The periodic boundary conditions must be appear in pairs."
    Stop
  End If
  !===========================================================================
  
  Allocate(check_periodic_boundary(2,num_periodic_boundary/2))
  check_periodic_boundary = -200
  
  Allocate(PairNodesTemp(2, nnx*nny)) ! big enough
  PairNodesTemp = -200
  
  num_corresponding_boundary=0
  num=0
  
  Allocate(match_flag(num_unknows))
  match_flag = 0
  
  Do i = 1, n_boundary
    
    If (ABS(ALQ_index(2,i)-x_min)<SmallValue .AND. ABS(ALQ_index(4,i)-x_min)<SmallValue) Then !vertical periodic boundary

      y1_value_min = MIN(ALQ_index(3,i),ALQ_index(5,i))
      y1_value_max = MAX(ALQ_index(3,i),ALQ_index(5,i))
      num_corresponding_boundary = num_corresponding_boundary + 1
      check_periodic_boundary(1,num_corresponding_boundary) = i

      Do j = 1, n_boundary

        y2_value_min = MIN(ALQ_index(3,j),ALQ_index(5,j))
        y2_value_max = MAX(ALQ_index(3,j),ALQ_index(5,j))

        If (ABS(ALQ_index(2,j)-x_max)<SmallValue .AND. ABS(ALQ_index(4,j)-x_max)<SmallValue .AND. &   ! match another periodic boundary
            ABS(y1_value_min-y2_value_min)<SmallValue .AND. ABS(y1_value_max-y2_value_max)<SmallValue) Then

          check_periodic_boundary(2,num_corresponding_boundary)=j

          Do n = 1, num_unknows
            
            If (ABS(p_basic(1,n)-x_min)<SmallValue .AND. p_basic(2,n)-y1_value_min>SmallValue .AND. &
                p_basic(2,n)-y1_value_max<-SmallValue) Then

              num = num + 1
              PairNodesTemp(1,num) = n  !The 'unknow_nodes(1,num)' is store the index of node
              match_flag(n) = 1 !periodic boundary node flag

              Do m = 1, num_unknows
                If (match_flag(m) == 0) Then
                  If (ABS(p_basic(1,m)-x_max)<SmallValue .AND. ABS(p_basic(2,n)-p_basic(2,m))<SmallValue) Then

                    PairNodesTemp(2,num) = m  !The 'unknow_nodes(2,num)' is store the index of brother node of the node in another periodic boundary
                    match_flag(m) = 1 !periodic boundary node flag
                    Exit
                  End If
                End If
              End Do
            End If
          End Do
        End If
      End Do

    Elseif(ABS(ALQ_index(3,i)-y_min)<SmallValue .AND. ABS(ALQ_index(5,i)-y_min)<SmallValue) Then !horizontal periodic boundary

      x1_value_min=MIN(ALQ_index(2,i),ALQ_index(4,i))
      x1_value_max=MAX(ALQ_index(2,i),ALQ_index(4,i))
      num_corresponding_boundary=num_corresponding_boundary+1
      check_periodic_boundary(1,num_corresponding_boundary)=i

      Do j=1,SIZE(ALQ_index,2)

        x2_value_min=MIN(ALQ_index(2,j),ALQ_index(4,j))
        x2_value_max=MAX(ALQ_index(2,j),ALQ_index(4,j))

        If (ABS(ALQ_index(3,j)-y_max)<SmallValue .AND. ABS(ALQ_index(5,j)-y_max)<SmallValue .AND. &
            ABS(x1_value_min-x2_value_min)<SmallValue .AND. ABS(x1_value_max-x2_value_max)<SmallValue) Then

          check_periodic_boundary(2,num_corresponding_boundary)=j

          Do n=1,num_unknows
            If (ABS(p_basic(2,n)-y_min)<SmallValue .AND. p_basic(1,n)-x1_value_min>SmallValue .AND. &
                p_basic(1,n)-x1_value_max<-SmallValue) Then

              num=num+1
              PairNodesTemp(1,num)=n  !The 'unknow_nodes(1,num)' is store the index of node
              match_flag(n) = 1 !periodic boundary node flag

              Do m=1,num_unknows
                If (match_flag(m) == 0) Then
                  If (ABS(p_basic(2,m)-y_max)<SmallValue .AND. ABS(p_basic(1,n)-p_basic(1,m))<SmallValue) Then

                    PairNodesTemp(2,num)=m  !The 'unknow_nodes(2,num)' is store the index of brother node of the node in another periodic boundary
                    match_flag(m) = 1 !periodic boundary node flag
                    Exit
                  End If
                End If
              End Do
            End If
          End Do
        End If
      End Do
    End If
  End Do
  
  !===================Check the boundary conditions===========================
  !The unknow_nodes default value is -200. After code running, the value is great zero
  Do i=1,num
    If (PairNodesTemp(1,i)<0 .OR. PairNodesTemp(2,i)<0) Then
      Write(*,*) "Please check the periodic boundary conditions or the start and end points of the boundary!"
      Write(*,*) "PairNodesTemp error."
      Stop
    End If
  End Do
  Do i=1,num_corresponding_boundary
    If (check_periodic_boundary(2,num_corresponding_boundary)<0 .OR. &
        check_periodic_boundary(1,num_corresponding_boundary)<0) Then
      Write(*,*) "Please check the periodic boundary conditions or the start and end points of the boundary!"
      Write(*,*) "check_periodic_boundary."
      Stop
    End If
  End Do
  !===========================================================================

  If (num > 0) Then
    Allocate(PairNodes(2,num))
    Do i = 1, num
      PairNodes(:,i) = PairNodesTemp(:,i)
    End Do
  End If

  !PairNodes:
  !(1,:) is the index of node in left or bottom periodic boundary
  !(2,:) is the index of brother node of the (1,:) in right or top periodic boundary
  !Size(PairNodes,2) is the node number of left or bottom periodic boundary

  !Note: Left and Right boundary is a pairs periodic boundary.
  !      Bottom and Top boundary is a pairs periodic boundary.
  Deallocate(ALQ_index, check_periodic_boundary, PairNodesTemp)
  
End Subroutine Periodic_pair_boundary_nodes