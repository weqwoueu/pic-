SUBROUTINE Update_Edge_Objects_Intersection_2D(DGE, DGP, DGT, information_2, node_index, element_index, objects, &
                                                information_3, edge_index)
!----------------------- wsy add information3 and edge_index 2021/7/30 -------------------------------------
!��һ�еڶ���Ϊ�߶˵��ȫ������
!������Ϊ�����ڵ�Ԫ
!�����е�����Ϊ��������
!�����е������Ǳ߶˵��λ����Ϣ out : -2 in : -1  �������Ѿ�ͳһ -2Ϊout -1Ϊin

USE Object_Data_2D
USE IFE_MAIN_PARAM !LY REVISE, 2022-1-6, Float number error

IMPLICIT NONE

INTEGER, DIMENSION(:), POINTER    ::  edge_index, node_index, element_index
REAL(8), DIMENSION(:,:), POINTER 	::	information_2, information_3, information_3_temp
REAL(8), DIMENSION(:,:), POINTER 	::	DGP
REAL(8), DIMENSION(:), POINTER    ::  information_2_vector, start_point, end_point
INTEGER, DIMENSION(:,:), POINTER	::	DGT, DGE
INTEGER                           ::  number_of_edge, i, j, N_Objects
INTEGER                           ::  xc, yc, r, out_index, in_index, count
TYPE(ObjectType), DIMENSION(:), POINTER		::	objects

number_of_edge = size(DGE, 2)
ALLOCATE(information_3_temp(7,number_of_edge))
ALLOCATE(information_2_vector(8))
ALLOCATE(edge_index(number_of_edge))
AllOCATE(start_point(2), end_point(2))
N_Objects = size(objects)
count = 0

!=========LY modification for Multi-Objects-Bool calculation, 2022-6-9=========
out_index = -1
Do i = 1, number_of_edge
  If (element_index(DGE(5,i)) == out_index) Then !Outside no-interface element.
    
    edge_index(i) = out_index
    
  Elseif (element_index(DGE(5,i)) < out_index) Then  !Inside no-interface element.
    
    edge_index(i) = element_index(DGE(5,i))
    
  Elseif (element_index(DGE(5,i)) > 0) Then !Interface element.
    
    information_2_vector(:) = information_2(1:8, element_index(DGE(5,i)))
    
    If (node_index(DGE(1,i))==out_index .AND. node_index(DGE(2,i))==out_index) Then   !node1 and node2 locate on outside objects.
      
      edge_index(i) = out_index
      
    Elseif (node_index(DGE(1,i))<out_index .AND. node_index(DGE(2,i))<out_index) Then   !node1 and node2 locate on inside objects.
      
      If (node_index(DGE(1,i))/=node_index(DGE(2,i))) Then
        Write(*,*) '=========error: The edge belong to two different objects.========='
        Write(*,*) 'The edge_index = ', i
        Write(*,*) 'The node_index1 /= node_index2.'
        Write(*,*) 'node_index1 = ', node_index(DGE(1,i))
        Write(*,*) 'node_index2 = ', node_index(DGE(2,i))
        Write(*,*) '=========error: The edge belong to two different objects.========='
        Stop
      End If
      
      edge_index(i) = node_index(DGE(1,i))
      
    Elseif (node_index(DGE(1,i))==out_index .AND. node_index(DGE(2,i))<out_index) Then    !node1 locate on outside objects and node2 locate on inside objects.
      
      If (ABS(DGP(1,DGE(2,i))-information_2_vector(3))<SmallValue .AND. &
          ABS(DGP(2,DGE(2,i))-information_2_vector(4))<SmallValue) Then

        edge_index(i) = out_index

      Elseif (ABS(DGP(1,DGE(2,i))-information_2_vector(5))<SmallValue .AND. &
              ABS(DGP(2,DGE(2,i))-information_2_vector(6))<SmallValue) Then

        edge_index(i) = out_index

      Else
        count = count + 1
        edge_index(i) = count
        start_point = DGP(1:2, DGE(1, i))
        end_point = DGP(1:2, DGE(2, i))
        information_3_temp(1:2, i) = DGE(1:2, i)
        information_3_temp(3, i) = DGE(5, i)
        information_3_temp(6, i) = out_index
        information_3_temp(7, i) = node_index(DGE(2,i))

        !LY REVISE, 2022-1-6, Float number error
        !initial code is 'start_point(1) == end_point(1)'
        If (ABS(start_point(1)-end_point(1)) < SmallValue) Then        !x1 = x2, this is vertical line

          !initial code is 'start_point(1) == information_2_vector(3)'
          If (ABS(start_point(1)-information_2_vector(3)) < SmallValue) Then     !the intersection point is D
            information_3_temp(4, i) = information_2_vector(3)
            information_3_temp(5, i) = information_2_vector(4)

            !initial code is 'start_point(1) == information_2_vector(5)'
          Elseif (ABS(start_point(1)-information_2_vector(5)) < SmallValue) Then !the intersection point is E
            information_3_temp(4, i) = information_2_vector(5)
            information_3_temp(5, i) = information_2_vector(6)
          End If

          !LY REVISE, 2022-1-6, Float number error
          !initial code is 'start_point(2) == end_point(2)'
        Elseif (ABS(start_point(2)-end_point(2)) < SmallValue) Then    !y1 = y2, this is horizonal line

          !initial code is 'start_point(2) == information_2_vector(4)'
          If (ABS(start_point(2)-information_2_vector(4)) < SmallValue) Then     !the intersection point is D
            information_3_temp(4, i) = information_2_vector(3)
            information_3_temp(5, i) = information_2_vector(4)

            !initial code is 'start_point(2) == information_2_vector(6)'
          Elseif (ABS(start_point(2)-information_2_vector(6)) < SmallValue) Then !the intersection point is E
            information_3_temp(4, i) = information_2_vector(5)
            information_3_temp(5, i) = information_2_vector(6)
          End If

        Else
          Print *,'error in Update_Edge_Objects_Intersection_2D'
          Write(6, *) 'start_point', start_point, 'end_point', end_point
          Stop
        End If
      End If
      
    Elseif (node_index(DGE(1,i))<out_index .AND. node_index(DGE(2,i))==out_index) Then    !node1 locate on inside objects and node2 locate on outside objects.
      
      If (ABS(DGP(1,DGE(1,i))-information_2_vector(3))<SmallValue .AND. &
          ABS(DGP(2,DGE(1,i))-information_2_vector(4))<SmallValue) Then

        edge_index(i) = out_index

      Elseif (ABS(DGP(1,DGE(1,i))-information_2_vector(5))<SmallValue .AND. &
              ABS(DGP(2,DGE(1,i))-information_2_vector(6))<SmallValue) Then

        edge_index(i) = out_index

      Else

        count = count + 1
        edge_index(i) = count
        start_point = DGP(:, DGE(1, i))
        end_point = DGP(:, DGE(2, i))
        information_3_temp(1:2, i) = DGE(1:2, i)
        information_3_temp(3, i) = DGE(5, i)
        information_3_temp(6, i) = node_index(DGE(1,i))
        information_3_temp(7, i) = out_index

        !LY REVISE, 2022-1-6, Float number error
        !initial code is 'start_point(1) == end_point(1)'
        If (ABS(start_point(1)-end_point(1)) < SmallValue) Then        !x1 = x2, this is vertical line

          !initial code is 'start_point(1) == information_2_vector(3)'
          If (ABS(start_point(1)-information_2_vector(3)) < SmallValue) Then       !the intersection point is D
            information_3_temp(4, i) = information_2_vector(3)
            information_3_temp(5, i) = information_2_vector(4)

            !initial code is 'start_point(1) == information_2_vector(5)'
          Elseif (ABS(start_point(1)-information_2_vector(5)) < SmallValue) Then   !the intersection point is E
            information_3_temp(4, i) = information_2_vector(5)
            information_3_temp(5, i) = information_2_vector(6)
          End If

          !LY REVISE, 2022-1-6, Float number error
          !initial code is 'start_point(2) == end_point(2)'
        Elseif (ABS(start_point(2)-end_point(2)) < SmallValue) Then    !y1 = y2, this is horizonal line

          !initial code is 'start_point(2) == information_2_vector(4)'
          If (ABS(start_point(2)-information_2_vector(4)) < SmallValue) Then       !the intersection point is D
            information_3_temp(4, i) = information_2_vector(3)
            information_3_temp(5, i) = information_2_vector(4)

            !initial code is 'start_point(2) == information_2_vector(6)'
          Elseif (ABS(start_point(2)-information_2_vector(6)) < SmallValue) Then   !the intersection point is E
            information_3_temp(4, i) = information_2_vector(5)
            information_3_temp(5, i) = information_2_vector(6)
          End If

        Else
          Print *,'error in Update_Edge_Objects_Intersection_2D'
          Write(6, *) 'start_point', start_point, 'end_point', end_point
          Stop
        End If
      End If
    End If
  End If
End Do
!=========LY modification for Multi-Objects-Bool calculation, 2022-7-25=========

!=========Old Code=========
!do i = 1, number_of_edge
!    do j = 1, N_Objects
!        
!        if (objects(j)%Shape==1 .AND. objects(j)%Axis==0) then !Circle  %Axisʲô��˼
!        
!            in_index = objects(j)%Regions(1)    !in_index = -2
!            out_index = objects(j)%Regions(2)   !out_index = -1
!            xc = objects(j)%Locations(1, 1)     !circle_x
!            yc = objects(j)%Locations(1, 2)     !circle_y
!            r = objects(j)%Dimensions(1)        !circle_r
!            
!            if (element_index(DGE(5, i)) == out_index) then       !outside objects
!                
!                edge_index(i) = out_index
!                
!            elseif (element_index(DGE(5, i)) == in_index) then    !inside objects
!                
!                edge_index(i) = in_index
!                
!            elseif (element_index(DGE(5, i)) > 0) then            !interface element
!                
!                information_2_vector(:) = information_2(:, element_index( DGE(5, i) ))
!                
!                if (node_index(DGE(1, i)) == out_index .AND. node_index(DGE(2, i)) == out_index) then     !node1 and node2 locate on outside objects
!                    
!                    edge_index(i) = out_index
!                    
!                elseif (node_index(DGE(1, i)) == in_index .AND. node_index(DGE(2, i)) == in_index) then   !node1 and node2 locate on inside objects
!                    
!                    edge_index(i) = in_index
!                    
!                elseif (node_index(DGE(1, i)) == out_index .AND. node_index(DGE(2, i)) == in_index) then  !node1 locate on outside objects and node2 locate on inside objects
!                    
!                    count = count + 1
!                    edge_index(i) = count
!                    start_point = DGP(:, DGE(1, i))
!                    end_point = DGP(:, DGE(2, i))
!                    information_3_temp(1:2, i) = DGE(1:2, i)
!                    information_3_temp(3, i) = DGE(5, i)
!                    information_3_temp(6, i) = out_index
!                    information_3_temp(7, i) = in_index
!                    
!                    !LY REVISE, 2022-1-6, Float number error
!                    !initial code is 'start_point(1) == end_point(1)'
!                    if (ABS(start_point(1)-end_point(1)) < SmallValue) then        !x1 = x2, this is vertical line
!                        
!                        !initial code is 'start_point(1) == information_2_vector(3)'
!                        if (ABS(start_point(1)-information_2_vector(3)) < SmallValue) then     !the intersection point is D
!                            information_3_temp(4, i) = information_2_vector(3)
!                            information_3_temp(5, i) = information_2_vector(4)
!                            
!                        !initial code is 'start_point(1) == information_2_vector(5)'
!                        elseif (ABS(start_point(1)-information_2_vector(5)) < SmallValue) then !the intersection point is E
!                            information_3_temp(4, i) = information_2_vector(5)
!                            information_3_temp(5, i) = information_2_vector(6)
!                        endif
!  
!                    !LY REVISE, 2022-1-6, Float number error
!                    !initial code is 'start_point(2) == end_point(2)'
!                    elseif (ABS(start_point(2)-end_point(2)) < SmallValue) then    !y1 = y2, this is horizonal line
!                        
!                        !initial code is 'start_point(2) == information_2_vector(4)'
!                        if (ABS(start_point(2)-information_2_vector(4)) < SmallValue) then     !the intersection point is D
!                            information_3_temp(4, i) = information_2_vector(3)
!                            information_3_temp(5, i) = information_2_vector(4)
!                            
!                        !initial code is 'start_point(2) == information_2_vector(6)'
!                        elseif (ABS(start_point(2)-information_2_vector(6)) < SmallValue) then !the intersection point is E
!                            information_3_temp(4, i) = information_2_vector(5)
!                            information_3_temp(5, i) = information_2_vector(6)
!                        endif
!                        
!                    else
!                        print *,'error in Update_Edge_Objects_Intersection_2D'
!                        WRITE(6, *) 'start_point', start_point, 'end_point', end_point
!                        WRITE(6, *) 'element_index =', element_index(DGE(5, i))
!                        WRITE(6, *) 'element = ', DGE(5, i)
!                        WRITE(6, *) 'node_index_1 =', node_index(DGE(1, i))
!                        WRITE(6, *) 'node_index_2 =', node_index(DGE(2, i))
!                        WRITE(6, *) 'edge_index =', i
!                        STOP
!                    endif
!                
!                elseif (node_index(DGE(1, i)) == in_index .AND. node_index(DGE(2, i)) == out_index) then  !node1 locate on inside objects and node2 locate on outside objects
!                    
!                    count = count + 1
!                    edge_index(i) = count
!                    start_point = DGP(:, DGE(1, i))
!                    end_point = DGP(:, DGE(2, i))
!                    information_3_temp(1:2, i) = DGE(1:2, i)
!                    information_3_temp(3, i) = DGE(5, i)
!                    information_3_temp(6, i) = in_index
!                    information_3_temp(7, i) = out_index
!                    
!                    !LY REVISE, 2022-1-6, Float number error
!                    !initial code is 'start_point(1) == end_point(1)'
!                    if (ABS(start_point(1)-end_point(1)) < SmallValue) then        !x1 = x2, this is vertical line
!                        
!                        !initial code is 'start_point(1) == information_2_vector(3)'
!                        if (ABS(start_point(1)-information_2_vector(3)) < SmallValue) then       !the intersection point is D
!                            information_3_temp(4, i) = information_2_vector(3)
!                            information_3_temp(5, i) = information_2_vector(4)
!                            
!                        !initial code is 'start_point(1) == information_2_vector(5)'
!                        elseif (ABS(start_point(1)-information_2_vector(5)) < SmallValue) then   !the intersection point is E
!                            information_3_temp(4, i) = information_2_vector(5)
!                            information_3_temp(5, i) = information_2_vector(6)
!                        endif
!  
!                    !LY REVISE, 2022-1-6, Float number error
!                    !initial code is 'start_point(2) == end_point(2)'
!                    elseif (ABS(start_point(2)-end_point(2)) < SmallValue) then    !y1 = y2, this is horizonal line
!                        
!                        !initial code is 'start_point(2) == information_2_vector(4)'
!                        if (ABS(start_point(2)-information_2_vector(4)) < SmallValue) then       !the intersection point is D
!                            information_3_temp(4, i) = information_2_vector(3)
!                            information_3_temp(5, i) = information_2_vector(4)
!                            
!                        !initial code is 'start_point(2) == information_2_vector(6)'
!                        elseif (ABS(start_point(2)-information_2_vector(6)) < SmallValue) then   !the intersection point is E
!                            information_3_temp(4, i) = information_2_vector(5)
!                            information_3_temp(5, i) = information_2_vector(6)
!                        endif
!                        
!                    else
!                        print *,'error in Update_Edge_Objects_Intersection_2D'
!                        WRITE(6, *) 'start_point', start_point, 'end_point', end_point
!                        WRITE(6, *) 'element_index =', element_index(DGE(5, i))
!                        WRITE(6, *) 'element = ', DGE(5, i)
!                        WRITE(6, *) 'node_index_1 =', node_index(DGE(1, i))
!                        WRITE(6, *) 'node_index_2 =', node_index(DGE(2, i))
!                        WRITE(6, *) 'edge_index =', i
!                        STOP
!                    endif
!                    
!                endif
!                
!            else
!                print *,'element_index error in Update_Edge_Objects_Intersection_2D'
!            endif
!            
!            
!        elseif (objects(j)%Shape==3) then   !Box LY 2021-11-19
!        
!            in_index = objects(j)%Regions(1)    !in_index = -2
!            out_index = objects(j)%Regions(2)   !out_index = -1
!            
!            if (element_index(DGE(5, i)) == out_index) then       !outside objects
!                
!                edge_index(i) = out_index
!                
!            elseif (element_index(DGE(5, i)) == in_index) then    !inside objects
!                
!                edge_index(i) = in_index
!                
!            elseif (element_index(DGE(5, i)) > 0) then            !interface element
!                
!                information_2_vector(:) = information_2(:, element_index( DGE(5, i) ))
!                
!                if (node_index(DGE(1, i)) == out_index .AND. node_index(DGE(2, i)) == out_index) then     !node1 and node2 locate on outside objects
!                    
!                    edge_index(i) = out_index
!                    
!                elseif (node_index(DGE(1, i)) == in_index .AND. node_index(DGE(2, i)) == in_index) then   !node1 and node2 locate on inside objects
!                    
!                    edge_index(i) = in_index
!                    
!                elseif (node_index(DGE(1, i)) == out_index .AND. node_index(DGE(2, i)) == in_index) then  !node1 locate on outside objects and node2 locate on inside objects
!                    
!                    count = count + 1
!                    edge_index(i) = count
!                    start_point = DGP(:, DGE(1, i))
!                    end_point = DGP(:, DGE(2, i))
!                    information_3_temp(1:2, i) = DGE(1:2, i)
!                    information_3_temp(3, i) = DGE(5, i)
!                    information_3_temp(6, i) = out_index
!                    information_3_temp(7, i) = in_index
!                    
!                    !LY REVISE, 2022-1-6, Float number error
!                    !initial code is 'start_point(1) == end_point(1)'
!                    if (ABS(start_point(1)-end_point(1)) < SmallValue) then        !x1 = x2, this is vertical line
!                        
!                        !initial code is 'start_point(1) == information_2_vector(3)'
!                        if (ABS(start_point(1)-information_2_vector(3)) < SmallValue) then     !the intersection point is D
!                            information_3_temp(4, i) = information_2_vector(3)
!                            information_3_temp(5, i) = information_2_vector(4)
!                            
!                        !initial code is 'start_point(1) == information_2_vector(5)'
!                        elseif (ABS(start_point(1)-information_2_vector(5)) < SmallValue) then !the intersection point is E
!                            information_3_temp(4, i) = information_2_vector(5)
!                            information_3_temp(5, i) = information_2_vector(6)
!                        endif
!  
!                    !LY REVISE, 2022-1-6, Float number error
!                    !initial code is 'start_point(2) == end_point(2)'
!                    elseif (ABS(start_point(2)-end_point(2)) < SmallValue) then    !y1 = y2, this is horizonal line
!                        
!                        !initial code is 'start_point(2) == information_2_vector(4)'
!                        if (ABS(start_point(2)-information_2_vector(4)) < SmallValue) then     !the intersection point is D
!                            information_3_temp(4, i) = information_2_vector(3)
!                            information_3_temp(5, i) = information_2_vector(4)
!                            
!                        !initial code is 'start_point(2) == information_2_vector(6)'
!                        elseif (ABS(start_point(2)-information_2_vector(6)) < SmallValue) then !the intersection point is E
!                            information_3_temp(4, i) = information_2_vector(5)
!                            information_3_temp(5, i) = information_2_vector(6)
!                        endif
!                        
!                    else
!                        print *,'error in Update_Edge_Objects_Intersection_2D'
!                        WRITE(6, *) 'start_point', start_point, 'end_point', end_point
!                        STOP
!                    endif
!                    Exit
!                elseif (node_index(DGE(1, i)) == in_index .AND. node_index(DGE(2, i)) == out_index) then  !node1 locate on inside objects and node2 locate on outside objects
!                    
!                    count = count + 1
!                    edge_index(i) = count
!                    start_point = DGP(:, DGE(1, i))
!                    end_point = DGP(:, DGE(2, i))
!                    information_3_temp(1:2, i) = DGE(1:2, i)
!                    information_3_temp(3, i) = DGE(5, i)
!                    information_3_temp(6, i) = in_index
!                    information_3_temp(7, i) = out_index
!                    
!                    !LY REVISE, 2022-1-6, Float number error
!                    !initial code is 'start_point(1) == end_point(1)'
!                    if (ABS(start_point(1)-end_point(1)) < SmallValue) then        !x1 = x2, this is vertical line
!                        
!                        !initial code is 'start_point(1) == information_2_vector(3)'
!                        if (ABS(start_point(1)-information_2_vector(3)) < SmallValue) then       !the intersection point is D
!                            information_3_temp(4, i) = information_2_vector(3)
!                            information_3_temp(5, i) = information_2_vector(4)
!                            
!                        !initial code is 'start_point(1) == information_2_vector(5)'
!                        elseif (ABS(start_point(1)-information_2_vector(5)) < SmallValue) then   !the intersection point is E
!                            information_3_temp(4, i) = information_2_vector(5)
!                            information_3_temp(5, i) = information_2_vector(6)
!                        endif
!  
!                    !LY REVISE, 2022-1-6, Float number error
!                    !initial code is 'start_point(2) == end_point(2)'
!                    elseif (ABS(start_point(2)-end_point(2)) < SmallValue) then    !y1 = y2, this is horizonal line
!                        
!                        !initial code is 'start_point(2) == information_2_vector(4)'
!                        if (ABS(start_point(2)-information_2_vector(4)) < SmallValue) then       !the intersection point is D
!                            information_3_temp(4, i) = information_2_vector(3)
!                            information_3_temp(5, i) = information_2_vector(4)
!                            
!                        !initial code is 'start_point(2) == information_2_vector(6)'
!                        elseif (ABS(start_point(2)-information_2_vector(6)) < SmallValue) then   !the intersection point is E
!                            information_3_temp(4, i) = information_2_vector(5)
!                            information_3_temp(5, i) = information_2_vector(6)
!                        endif
!                        
!                    else
!                        print *,'error in Update_Edge_Objects_Intersection_2D'
!                        WRITE(6, *) 'start_point', start_point, 'end_point', end_point
!                        STOP
!                    endif
!                    
!                endif
!                
!            else
!                print *,'element_index error in Update_Edge_Objects_Intersection_2D'
!            endif
!            Exit
!        endif
!    
!    end do
!end do
!=========Old Code=========

ALLOCATE(information_3(7, count))
count = 1
do i = 1, number_of_edge
    if (edge_index(i) > 0) then 
        information_3(: , count) = information_3_temp(:, i)
        count = count + 1
    endif
end do
    
END