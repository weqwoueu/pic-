SUBROUTINE generate_Adaptive_DGT_IFE(DGP, DGT, objects, DGT_IFE_partition, area_flag)
  
!------------------LY add DGT_IFE for Adaptive Mesh 2021-9-25------------------
!Note: DGT_IFE can be seen as refinement flag.

USE Object_Data_2D
USE IFE_MAIN_PARAM
USE IFE_INTERFACE, ONLY: Adaptive_interface_function

IMPLICIT NONE

REAL(8), DIMENSION(:,:), POINTER :: DGP
INTEGER, DIMENSION(:,:), POINTER :: DGT, DGT_IFE_partition
TYPE(ObjectType), DIMENSION(:), INTENT(IN) :: objects
INTEGER :: area_flag   !LY REVISE, 2021-12-19

INTEGER :: row_number_of_T_partition, number_of_element, count, n, k, i, n_objects, num_of_area, e
REAL(8) :: r
REAL(8) :: vertices(2,4)
REAL(8) :: Box_low_corner_x, Box_low_corner_y, Box_hig_corner_x, Box_hig_corner_y
REAL(8) :: local_1_x, local_1_y, local_2_x, local_2_y, local_3_x, local_3_y, local_4_x, local_4_y
REAL(8) :: length1y, length2y, length3y, length4y, length1x, length2x, length3x, length4x
REAL(8) :: left_flag1, right_flag1, bottom_flag1, top_flag1   !LY REVISE, 2021-12-19
REAL(8) :: left_flag2, right_flag2, bottom_flag2, top_flag2   !LY REVISE, 2022-01-11
REAL(8) :: r0, h_size, sqrt_2, xc, yc, r1, r2 !LY REVISE, 2022-1-10
REAL(8), Allocatable :: narrow_size_left(:), narrow_size_right(:), narrow_size_bottom(:), narrow_size_top(:) !WSY ADD, 2022-10-19
Integer :: narrow_times, index, j
!=========LY modification, 2022-7-25=========
Integer :: in_index, out_index
Real(8) :: ax, ay, bx, by, cx, cy, dx, dy
Real(8) :: S123, S134, S234, S124, Sq, Sp, Spab, Spbc, Spcd, Spda, Spbd
Real(8) :: mark_abc, mark_bcd, mark_cda, mark_dab
Real(8) :: x, y
Real(8) :: cpx1, cpy1, cpx2, cpy2, cpx3, cpy3
Real(8) :: a1, b1, c1, a2, b2, c2
Real(8) :: xr, yr, radius
Real(8) :: vecACx, vecACy, vecAZx, vecAZy, vecANx, vecANy
Real(8) :: vecAN_len, vecAC_len, vecAZ_len
Real(8) :: vecAN_dot_vecAC, vecAN_dot_vecAZ
Integer, Dimension(:), Allocatable :: node_location
!=========LY modification, 2022-7-25=========


REAL(8), DIMENSION(4) :: temp
DATA temp /0.0,0.0,0.0,0.0/

n = 0
k = 0
r = 0.0
vertices = 0.0

row_number_of_T_partition = 4
number_of_element = SIZE(DGT, 2)
n_objects	=	SIZE(objects,1)
ALLOCATE(DGT_IFE_partition(row_number_of_T_partition+1, number_of_element))
DGT_IFE_partition(:,:) = -3

!=========LY modification, 2022-7-25=========
Allocate(node_location(Size(DGP,2)))
node_location = -1
!=========LY modification, 2022-7-25=========

n = 0
DO n = 1, number_of_element
  DGT_IFE_partition(1:row_number_of_T_partition, n) = DGT(1:row_number_of_T_partition, n)
END DO


OPEN(1,ACTION='READ',FILE='./INPUT/IDG_inf.inp')
READ(1, *)
READ(1, *)
READ(1, *) num_of_area

Do e = 1, num_of_area

READ(1, *) area_flag

!LY REVISE, 2021-12-19
  IF (area_flag == 1) THEN
  !----------judge refinement element: interface element----------
    !=========LY modification, 2022-7-25=========
    !=========Part 0 : Initilization=========
    in_index = -2
    out_index = -1
    count = 1
    !=========Part 0 : Initilization=========
    
    !=========Part 1 : Judge type and location of node=========
    Do i = 1, n_objects
      Do n = 1, number_of_element

        If (objects(i)%Shape==1 .AND. objects(i)%Axis==0) Then    !Intact Circle

          Do k = 1, 4
            If (node_location(DGT(k,n)) == -1) Then
              Call Adaptive_interface_function(DGP(1,DGT(k, n)), DGP(2,DGT(k,n)), objects(i), r)
              temp(k) = r

              If (temp(k) < -SmallValue) Then
                !The node locate inside circle.
                node_location(DGT(k,n)) = in_index
              Elseif (ABS(temp(k)) < SmallValue) Then
                !We consider the node on circle as inner node.
                node_location(DGT(k,n)) = in_index
              Else
                !The node locate outside circle.
                node_location(DGT(k,n)) = out_index
              End If
            End If
          End Do

        Elseif (objects(i)%Shape==3) Then   !Box or Rectangle

          vertices = DGP(1:2, DGT(1:4, n))

          Box_low_corner_x = objects(i)%Locations(1, 1)
          Box_low_corner_y = objects(i)%Locations(1, 2)
          Box_hig_corner_x = objects(i)%Locations(2, 1)
          Box_hig_corner_y = objects(i)%Locations(2, 2)

          Do k = 1, 4
            If (node_location(DGT(k,n)) == -1) Then
              x = DGP(1, DGT(k,n))
              y = DGP(2, DGT(k,n))

              If (x>(Box_low_corner_x-SmallValue) .AND. x<(Box_hig_corner_x+SmallValue) .AND. &
                  y>(Box_low_corner_y-SmallValue) .AND. y<(Box_hig_corner_y+SmallValue)) Then
                !That is 'x-object%Locations(1,1) > -SmallValue .AND. x-object%Locations(2,1) < SmallValue. AND.&'
                !        'y-object%Locations(1,2) > -SmallValue .AND. y-object%Locations(2,2) < SmallValue.
                !So, the node on Box line will be seen inner node.
                node_location(DGT(k,n)) = in_index
              Else
                !The node locate outside Box.
                node_location(DGT(k,n)) = out_index
              End If

            End If
          End Do

        Elseif (objects(i)%Shape==4) Then   !Triangle or Quardangle

          vertices = DGP(1:2, DGT(1:4, n))    !vertices(2,4)

          ax = objects(i)%Locations(1,1)
          ay = objects(i)%Locations(1,2)
          bx = objects(i)%Locations(2,1)
          by = objects(i)%Locations(2,2)
          cx = objects(i)%Locations(3,1)
          cy = objects(i)%Locations(3,2)
          dx = objects(i)%Locations(4,1)
          dy = objects(i)%Locations(4,2)

          !Firstly, we need to calculate some triangle area.
          !Note: In Concave quadrilateral(or ao quadrilateral), the inner node is 1 node, the bottom nodes is 2 and 4 node, the top node is 3 node.
          S123 = DABS(0.5*(ax*by + ay*cx + bx*cy - cx*by - cy*ax - bx*ay))  !left triangle of tu quadrilateral
          S134 = DABS(0.5*(ax*cy + ay*dx + cx*dy - dx*cy - dy*ax - cx*ay))  !right triangle of tu quadrilateral
          S234 = DABS(0.5*((cx - bx)*(dy - by) - (dx - bx)*(cy - by)))  !big triangle of ao quadrilateral
          S124 = DABS(0.5*((bx - ax)*(dy - ay) - (dx - ax)*(by - ay)))  !small triangle of ao quadrilateral(it will be minus)

          !the area of tu quadrilateral.
          Sq = S123 + S134  !Sq=S123+S134

          !Secondly, we need to judge the type of quadrilateral: tu or ao?
          mark_abc = ax*by + ay*cx + bx*cy - cx*by - cy*ax - bx*ay  !abc vector product
          mark_bcd = bx*cy + by*dx + cx*dy - dx*cy - dy*bx - cx*by  !bcd vector product
          mark_cda = cx*dy + cy*ax + dx*ay - ax*dy - ay*cx - dx*cy  !cda vector product
          mark_dab = dx*ay + dy*bx + ax*by - bx*ay - by*dx - ax*dy  !dab vector product

          !=========LY modification for Bool calculation, 2022-6-4=========
          Do k = 1, 4
            If (node_location(DGT(k,n)) == -1) Then
              x = vertices(1,k)
              y = vertices(2,k)

              !the triangle area between of mesh node and tu quadrilateral node.
              Spab = DABS(0.5*((ax - x)*(by - y) - (bx - x)*(ay - y)))  !P-ab triangle
              Spbc = DABS(0.5*((bx - x)*(cy - y) - (cx - x)*(by - y)))  !P-bc triangle
              Spcd = DABS(0.5*((cx - x)*(dy - y) - (dx - x)*(cy - y)))  !P-cd triangle
              Spda = DABS(0.5*((dx - x)*(ay - y) - (ax - x)*(dy - y)))  !P-da triangle

              !the triangle area between of mesh node and ao quadrilateral bottom node(2 and 4 node).
              Spbd = DABS(0.5*((bx - x)*(dy - y) - (dx - x)*(by - y)))  !P-bd triangle

              !the area of four triangle from mesh node and tu quadrilateral node.
              Sp = Spab + Spbc + Spcd + Spda

              !Finally, we can judge the type of quadrilateral and the node location.
              If ((mark_abc > -SmallValue .AND. mark_bcd > -SmallValue .AND. mark_cda > -SmallValue .AND. &
                mark_dab > -SmallValue) .OR. (mark_abc < SmallValue .AND. mark_bcd < SmallValue .AND. &
                mark_cda < SmallValue .AND. mark_dab < SmallValue)) Then   !tu quadrilateral
                If (DABS(Sp - Sq) < SmallValue) Then    !LY modification, 2022-6-4, we consider the objects endpoint and node on objects line as inner node.
                  node_location(DGT(k,n)) = in_index   !node locate inside tu quadrilateral
                Else
                  node_location(DGT(k,n)) = out_index   !node locate outside tu quadrilateral
                End If

              Else  !ao quadrilateral
                If (DABS(Spbd + Spbc + Spcd - S234)<SmallValue) Then  !node locate inside big triangle of ao quadrilateral.
                  !LY modification, 2022-6-4, we consider the object endpoint and node on objects line as inner node.
                  If ((DABS(Spbd + Spbc + Spcd - S234)<SmallValue) .AND. &
                    (DABS(Spda) > SmallValue) .AND. (DABS(Spab) > SmallValue)) Then
                    node_location(DGT(k,n)) = out_index  !node locate inside small triangle of ao quadrilateral.
                  Else
                    node_location(DGT(k,n)) = in_index  !node locate inside the triangle between big triangle minus small triangle of ao quadrilateral.
                  End If
                Else
                  node_location(DGT(k,n)) = out_index  !node locate outside big triangle of ao quadrilateral
                End If
              End If
            End If
          End Do
          !=========LY modification for Bool calculation, 2022-6-4=========

        Elseif (objects(i)%Shape==5) Then   !Ellipse

          Do k = 1, 4
            If (node_location(DGT(k,n)) == -1) Then
              Call Adaptive_interface_function(DGP(1,DGT(k, n)), DGP(2,DGT(k,n)), objects(i), r)
              temp(k) = r
              
              If (temp(k)-1.0 < -1.0D-8) Then   !LY modification for replace '-SmallValue' with '-1.0D-8'.
                !The node locate inside ellipse.
                node_location(DGT(k,n)) = in_index
              Elseif (ABS(temp(k)-1.0) < 1.0D-8) Then   !LY modification for replace 'SmallValue' with '1.0D-8'.
                !We consider the node on ellipse as inner node.
                node_location(DGT(k,n)) = in_index
              Else
                !The node locate outside ellipse.
                node_location(DGT(k,n)) = out_index
              End If
            End If
          End Do

        Elseif (objects(i)%Shape==6) Then   !Arch

          vertices = DGP(1:2, DGT(1:4, n))

          !=========LY modification, 2022-5-31=========
          cpx1=objects(i)%Locations(1,1)  !X coordinate of the first node  A
          cpy1=objects(i)%Locations(1,2)  !Y coordinate of the first node  A
          cpx2=objects(i)%Locations(2,1)  !X coordinate of the second node B
          cpy2=objects(i)%Locations(2,2)  !Y coordinate of the second node B
          cpx3=objects(i)%Locations(3,1)  !X coordinate of the third node  C
          cpy3=objects(i)%Locations(3,2)  !Y coordinate of the third node  C

          vecACx = cpx3 - cpx1    !vector AC x component
          vecACy = cpy3 - cpy1    !vector AC y component
          vecAZx = cpy3 - cpy1    !vector AZ x component
          vecAZy = cpx1 - cpx3    !vector AZ y component
          !=========LY modification, 2022-5-31=========

          a1=2*(cpx2-cpx1)
          b1=2*(cpy2-cpy1)
          c1=cpx2*cpx2+cpy2*cpy2-cpx1*cpx1-cpy1*cpy1
          a2=2*(cpx3-cpx2)
          b2=2*(cpy3-cpy2)
          c2=cpx3*cpx3+cpy3*cpy3-cpx2*cpx2-cpy2*cpy2

          xr=((c1*b2)-(c2*b1))/((a1*b2)-(a2*b1))    !Centre x of circle of the arch
          yr=((a1*c2)-(a2*c1))/((a1*b2)-(a2*b1))    !Centre y of circle of the arch
          radius=DSQRT((yr-cpy1)*(yr-cpy1)+(xr-cpx1)*(xr-cpx1))   !Radius of circle of the arch

          Do k = 1, 4
            x = vertices(1,k)
            y = vertices(2,k)

            If (node_location(DGT(k,n)) == -1) Then
              If (DSQRT((x-xr)*(x-xr)+(y-yr)*(y-yr))-radius < SmallValue) Then   !Node locate inside circle of the arch.

                !=========LY modification, 2022-5-31=========
                vecANx = x - cpx1   !vector AN x component
                vecANy = y - cpy1   !vector AN y component
                vecAN_len = DSQRT(vecANx**2 + vecANy**2)    !vector AN length
                vecAC_len = DSQRT(vecACx**2 + vecACy**2)    !vector AC length
                vecAZ_len = DSQRT(vecAZx**2 + vecAZy**2)    !vector AZ length
                vecAN_dot_vecAC = vecANx * vecACx + vecANy * vecACy   !dot product between vector AN and AC
                vecAN_dot_vecAZ = vecANx * vecAZx + vecANy * vecAZy   !dot product between vector AN and AZ

                !Node locate inside arch only cos¦È1 = (AN¡¤AC)/(|AN|*|AC|)>0 and cos¦È2 = (AN¡¤AZ)/(|AN|*|AZ|)>0.
                !At same time, we consider node locate inside arch when node is A, B and C.
                If (vecAN_len < SmallValue .OR. ABS(vecAN_dot_vecAZ) < SmallValue) Then
                  !|AN|=0, that is the node N is node A. We consider the node locate inside arch.
                  !AN¡¤AZ=0, that is the node N inside AC. We consider the node locate inside arch.
                  node_location(DGT(k,n)) = in_index
                Else
                  If ((vecAN_dot_vecAC/vecAN_len/vecAC_len > SmallValue) .AND. &
                    (vecAN_dot_vecAZ/vecAN_len/vecAZ_len > SmallValue)) Then   !Node locate inside arch.
                    node_location(DGT(k,n)) = in_index
                  Else    !Node locate outside arch.
                    node_location(DGT(k,n)) = out_index
                  End If
                End If
                !=========LY modification, 2022-5-31=========

              Else   !Node locate outside circle of the arch
                node_location(DGT(k,n)) = out_index
              End If
            End If
          End Do
        End If

      End Do
    End Do
    !=========Part 1 : Judge type and location of node=========
    
    !=========Part 2 : Judge refinement element=========
    Do n = 1, number_of_element
      If (node_location(DGT(1,n)) + node_location(DGT(2,n)) + node_location(DGT(3,n)) + node_location(DGT(4,n)) &
          == 4*in_index) Then
        !inside no-interface element.
        DGT_IFE_partition(row_number_of_T_partition+1, n) = max(-1, DGT_IFE_partition(row_number_of_T_partition+1, n))

      Elseif (node_location(DGT(1,n)) + node_location(DGT(2,n)) + node_location(DGT(3,n)) + node_location(DGT(4,n)) &
              == 4*out_index) Then
        !outside no-interface element.
        DGT_IFE_partition(row_number_of_T_partition+1, n) = max(-2, DGT_IFE_partition(row_number_of_T_partition+1, n))

      Else
        DGT_IFE_partition(row_number_of_T_partition+1, n) = count
        count = count + 1
      End If
    End Do
    !=========Part 2 : Judge refinement element=========
    
    !!=========Old Code=========
    !count = 1
    !DO i = 1, n_objects
    !  DO n = 1, number_of_element
    !
    !    IF (objects(i)%Shape==1 .AND. objects(i)%Axis==0) THEN
    !
    !      DO k = 1, 4
    !        CALL Adaptive_interface_function(DGP(1,DGT(k, n)), DGP(2,DGT(k,n)), objects(i), r)
    !        temp(k) = r
    !      END DO
    !
    !      IF (temp(1)<=0 .AND. temp(2)<=0 .AND. temp(3)<=0 .AND. temp(4)<=0) THEN
    !        DGT_IFE_partition(row_number_of_T_partition+1, n) = -1
    !      ELSE IF (temp(1)>=0 .AND. temp(2)>=0 .AND. temp(3)>=0 .AND. temp(4)>=0) THEN
    !        DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
    !      ELSE
    !        DGT_IFE_partition(row_number_of_T_partition+1, n) = count
    !        count = count + 1
    !      END IF
    !  
    !    ELSEIF (objects(i)%Shape==3) THEN
    !      
    !      vertices = DGP(1:2, DGT(1:4, n))
    !      Box_low_corner_x = objects(i)%Locations(1, 1)
    !      Box_low_corner_y = objects(i)%Locations(1, 2)
    !      Box_hig_corner_x = objects(i)%Locations(2, 1)
    !      Box_hig_corner_y = objects(i)%Locations(2, 2)
    !      local_1_x = vertices(1,1)
    !      local_1_y = vertices(2,1)
    !      local_2_x = vertices(1,2)
    !      local_2_y = vertices(2,2)
    !      local_3_x = vertices(1,3)
    !      local_3_y = vertices(2,3)
    !      local_4_x = vertices(1,4)
    !      local_4_y = vertices(2,4)
    !  
    !      length1y = local_1_y - Box_hig_corner_y
    !      length2y = local_3_y - Box_low_corner_y
    !      length3y = local_1_y - Box_low_corner_y
    !      length4y = local_3_y - Box_hig_corner_y
    !  
    !      length1x = local_1_x - Box_hig_corner_x
    !      length2x = local_2_x - Box_low_corner_x
    !      length3x = local_1_x - Box_low_corner_x
    !      length4x = local_2_x - Box_hig_corner_x
    !
    !      ! write(*,*) Box_low_corner_x, Box_hig_corner_x
    !      ! write(*,*) Box_low_corner_y, Box_hig_corner_y
    !      ! stop
    !      !IF (length1y > SmallValue) THEN
    !      !
    !      !  DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
    !      !
    !      !ELSEIF (length2y < -SmallValue) THEN
    !      !
    !      !  DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
    !      !
    !      !ELSEIF (length3y > SmallValue .AND. length4y < -SmallValue) THEN
    !      !
    !      !  IF (length1x > SmallValue) THEN
    !      !    DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
    !      !  ELSEIF (length2x < -SmallValue) THEN
    !      !    DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
    !      !  ELSEIF (length3x > SmallValue .AND. length4x < -SmallValue) THEN
    !      !    DGT_IFE_partition(row_number_of_T_partition+1, n) = -1
    !      !  ELSE
    !      !    DGT_IFE_partition(row_number_of_T_partition+1, n) = count
    !      !    count = count + 1
    !      !  END IF
    !      !
    !      !ELSE
    !      !
    !      !  IF (length1x > SmallValue) THEN
    !      !    DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
    !      !  ELSEIF (length2x < -SmallValue) THEN
    !      !    DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
    !      !  ELSE
    !      !    DGT_IFE_partition(row_number_of_T_partition+1, n) = count
    !      !    count = count + 1
    !      !  END IF
    !      !
    !      !END IF
    !      
    !      IF (length1y >= 0.0) THEN
    !  
    !        DGT_IFE_partition(row_number_of_T_partition+1, n) = max(-2, DGT_IFE_partition(row_number_of_T_partition+1, n))
    !    
    !      ELSEIF (length2y <= 0.0) THEN
    !  
    !        DGT_IFE_partition(row_number_of_T_partition+1, n) = max(-2, DGT_IFE_partition(row_number_of_T_partition+1, n))
    !    
    !      ELSEIF (length3y >= 0.0 .AND. length4y <= 0.0) THEN
    !  
    !        IF (length1x >= 0.0) THEN
    !          DGT_IFE_partition(row_number_of_T_partition+1, n) = max(-2, DGT_IFE_partition(row_number_of_T_partition+1, n))
    !        ELSEIF (length2x <= 0.0) THEN
    !          DGT_IFE_partition(row_number_of_T_partition+1, n) = max(-2, DGT_IFE_partition(row_number_of_T_partition+1, n))
    !        ELSEIF (length3x >= 0.0 .AND. length4x <= 0.0) THEN
    !          DGT_IFE_partition(row_number_of_T_partition+1, n) = max(-1, DGT_IFE_partition(row_number_of_T_partition+1, n))
    !        ELSE
    !          
    !          DGT_IFE_partition(row_number_of_T_partition+1, n) = count
    !          count = count + 1
    !        END IF
    !    
    !      ELSE
    !  
    !        IF (length1x >= 0.0) THEN
    !          DGT_IFE_partition(row_number_of_T_partition+1, n) = max(-2, DGT_IFE_partition(row_number_of_T_partition+1, n))
    !        ELSEIF (length2x <= 0.0) THEN
    !          DGT_IFE_partition(row_number_of_T_partition+1, n) = max(-2, DGT_IFE_partition(row_number_of_T_partition+1, n))
    !        ELSE
    !          
    !          DGT_IFE_partition(row_number_of_T_partition+1, n) = count
    !          count = count + 1
    !        END IF
    !    
    !      END IF
    !  
    !    END IF
    !
    !  END DO
    !END DO
    !!=========Old Code=========
  !----------judge refinement element: interface element----------

  ELSEIF (area_flag == 2) THEN
  !----------judge refinement element: designated area------------
    count = 1

      READ(1, *) left_flag1, left_flag2
      READ(1, *) right_flag1, right_flag2
      READ(1, *) bottom_flag1, bottom_flag2
      READ(1, *) top_flag1, top_flag2
      
      DO n = 1, number_of_element
        !IF ((left_flag1<=DGP(1,DGT(1,n)).AND.DGP(1,DGT(2,n))<=right_flag1) .OR. &
        !    (left_flag2<=DGP(1,DGT(1,n)).AND.DGP(1,DGT(2,n))<=right_flag2) .OR. &
        !    (bottom_flag1<=DGP(2,DGT(1,n)).AND.DGP(2,DGT(4,n))<=top_flag1) .OR. &
        !    (bottom_flag2<=DGP(2,DGT(1,n)).AND.DGP(2,DGT(4,n))<=top_flag2)) THEN
        !!IF (left_flag<=DGP(1,DGT(1,n)) .AND. DGP(1,DGT(2,n))<=right_flag .AND. bottom_flag<=DGP(2,DGT(1,n)) .AND. DGP(2,DGT(4,n))<=top_flag) THEN
        !  DGT_IFE_partition(row_number_of_T_partition+1, n) = count
        !  count = count + 1
        !ELSE
        !  DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
        !END IF
        
        !designated area--outside boundary refinement!
        !IF (left_flag1<=DGP(1,DGT(1,n)) .AND. DGP(1,DGT(2,n))<=right_flag1 .AND. &
        !    bottom_flag1<=DGP(2,DGT(1,n)) .AND. DGP(2,DGT(4,n))<=top_flag1) THEN
        !  DGT_IFE_partition(row_number_of_T_partition+1, n) = -1
        !ELSE
        !  DGT_IFE_partition(row_number_of_T_partition+1, n) = count
        !  count = count + 1
        !END IF
        
        !designated area--inside refinement!
        IF (left_flag1<=DGP(1,DGT(1,n)) .AND. DGP(1,DGT(2,n))<=right_flag1 .AND. &
            bottom_flag1<=DGP(2,DGT(1,n)) .AND. DGP(2,DGT(4,n))<=top_flag1) THEN
          DGT_IFE_partition(row_number_of_T_partition+1, n) = count
          count = count + 1
        ELSE
          DGT_IFE_partition(row_number_of_T_partition+1, n) = max(DGT_IFE_partition(row_number_of_T_partition+1, n) , -2)
        END IF
        
      END DO

    ! Do i = 1, number_of_element
    !   if (DGT_IFE_partition(row_number_of_T_partition+1, i) > 0) THEN
    !     write(*,*) i, DGT_IFE_partition(row_number_of_T_partition+1, i)
    !   end IF
    ! end DO 
    ! stop
  !----------judge refinement element: designated area------------

  ELSEIF (area_flag == 3) THEN
  !---------judge refinement element: SingleParticle for neighbor element----------
    count = 1
    OPEN(1,ACTION='READ',FILE='IDG_inf.inp')
      READ(1, *)
      READ(1, *)
      READ(1, *)
      READ(1, *) 
      READ(1, *) 
      READ(1, *) 
      READ(1, *) 
      READ(1, *) h_size
    CLOSE(1)
    DO i = 1, n_objects
      DO n = 1, number_of_element
        
        IF (objects(i)%Shape==1 .AND. objects(i)%Axis==0) THEN
          r0 = objects(i)%Dimensions(1)   !circle center.
          sqrt_2 = 1.414213562
          DO k = 1, 4
            CALL Adaptive_interface_function(DGP(1,DGT(k, n)), DGP(2,DGT(k,n)), objects(i), r)
            temp(k) = r
          END DO

          IF (temp(1)<=(2*(h_size**2)-2*sqrt_2*h_size*r0) .AND. temp(2)<=(2*(h_size**2)-2*sqrt_2*h_size*r0) .AND.&
              temp(3)<=(2*(h_size**2)-2*sqrt_2*h_size*r0) .AND. temp(4)<=(2*(h_size**2)-2*sqrt_2*h_size*r0)) THEN
            
            DGT_IFE_partition(row_number_of_T_partition+1, n) = -1
            
          ELSE IF (temp(1)>=(2*(h_size**2)+2*sqrt_2*h_size*r0) .AND. temp(2)>=(2*(h_size**2)+2*sqrt_2*h_size*r0) .AND.&
              temp(3)>=(2*(h_size**2)+2*sqrt_2*h_size*r0) .AND. temp(4)>=(2*(h_size**2)+2*sqrt_2*h_size*r0)) THEN

            DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
            
          ELSE
            DGT_IFE_partition(row_number_of_T_partition+1, n) = count
            count = count + 1
          END IF
          
        END IF
      END DO
    END DO
  !---------judge refinement element: SingleParticle for neighbor element----------
    
  ELSEIF (area_flag == 4) THEN
  !---------judge refinement element: when no objects, refine for FE element----------
    count = 1
    
    OPEN(1,ACTION='READ',FILE='IDG_inf.inp')
      READ(1, *)
      READ(1, *)
      READ(1, *)
      READ(1, *) 
      READ(1, *) 
      READ(1, *) 
      READ(1, *) 
      READ(1, *)
      READ(1, *) xc, yc, r1, r2
    CLOSE(1)
    
    DO n = 1, number_of_element
      
      DO k = 1, 4
        temp(k) = SQRT((DGP(1, DGT(k,n))-xc)**2 + (DGP(2, DGT(k,n))-yc)**2)
      END DO
      
      IF (temp(1)<=r1 .AND. temp(2)<=r1 .AND. temp(3)<=r1 .AND. temp(4)<=r1) THEN
        DGT_IFE_partition(row_number_of_T_partition+1, n) = -1
      ELSEIF (temp(1)>=r2 .AND. temp(2)>=r2 .AND. temp(3)>=r2 .AND. temp(4)>=r2) THEN
        DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
      ELSE
        DGT_IFE_partition(row_number_of_T_partition+1, n) = count
        count = count + 1
      END IF
      
    END DO
    
  !---------judge refinement element: when no objects, refine for FE element----------
    ELSE IF (area_flag == 5) THEN
    !----------Non-structure mesh, we probably need step by step refine mesh.----------
    !when we finish one refinement, we need change refine domain.
  
        count = 1
        OPEN(1,ACTION='READ',FILE='./INPUT/IDG_inf_step_refine.inp')
        READ(1, *) narrow_times
        READ(1, *) left_flag1, left_flag2
        READ(1, *) right_flag1, right_flag2
        READ(1, *) bottom_flag1, bottom_flag2
        READ(1, *) top_flag1, top_flag2
        Allocate(narrow_size_left(narrow_times+1))
        Allocate(narrow_size_right(narrow_times+1))
        Allocate(narrow_size_bottom(narrow_times+1))
        Allocate(narrow_size_top(narrow_times+1))
        READ(1, *) (narrow_size_left(j), j=1,narrow_times+1)
        READ(1, *) (narrow_size_right(j), j=1,narrow_times+1)
        READ(1, *) (narrow_size_bottom(j), j=1,narrow_times+1)
        READ(1, *) (narrow_size_top(j), j=1,narrow_times+1)
        READ(1, *) index
        CLOSE(1)
  
        DO n = 1, number_of_element
        IF (left_flag1<=DGP(1,DGT(1,n)) .AND. DGP(1,DGT(2,n))<=right_flag1 .AND. &
            bottom_flag1<=DGP(2,DGT(1,n)) .AND. DGP(2,DGT(4,n))<=top_flag1) THEN
            DGT_IFE_partition(row_number_of_T_partition+1, n) = count
            count = count + 1
        ELSE
            DGT_IFE_partition(row_number_of_T_partition+1, n) = -2
        END IF
        END DO
  
        left_flag1 = left_flag1 + narrow_size_left(index)
        right_flag1 = right_flag1 + narrow_size_right(index)
        bottom_flag1 = bottom_flag1 + narrow_size_bottom(index)
        top_flag1 = top_flag1 + narrow_size_top(index)
  
        OPEN(1,ACTION='WRITE',FILE='./INPUT/IDG_inf_step_refine.inp')
        WRITE(1, *) narrow_times
        WRITE(1, *) left_flag1, left_flag2
        WRITE(1, *) right_flag1, right_flag2
        WRITE(1, *) bottom_flag1, bottom_flag2
        WRITE(1, *) top_flag1, top_flag2
        WRITE(1, *) (narrow_size_left(j), j=1,narrow_times+1)
        WRITE(1, *) (narrow_size_right(j), j=1,narrow_times+1)
        WRITE(1, *) (narrow_size_bottom(j), j=1,narrow_times+1)
        WRITE(1, *) (narrow_size_top(j), j=1,narrow_times+1)
        WRITE(1, *) index + 1
        CLOSE(1)
  END IF

End Do
Close(1)

!=========LY modification, 2022-7-25=========
Deallocate(node_location)
!=========LY modification, 2022-7-25=========

! ------------------ wsy add for test(can delete) ------
!DGT_IFE_partition(5, :) = -2
!DO n = 19, 20
!  DGT_IFE_partition(5, n) = n
!END DO
!--------------------------------------------------------------------------

END SUBROUTINE