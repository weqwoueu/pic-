SUBROUTINE Update_Mesh_Objects_Intersection_2D  (	t_basic, p_basic, objects,		&
												element_index, information_1, information_2, node_index)
! Purpose:		Intersect an array of objects (generalized objects) with a Cartesian mesh.
! Last Update:	4/20/2004 06:01 PM

!USE PIC_MAIN_PARAM_2D
USE IFE_MAIN_PARAM
USE Object_Data_2D
USE Wall_2D
USE IFE_Interface, ONLY: El_Object_Intersection_info_2D_BL,El_Curve_Intersection_info_2D_BL,   &
                            Update_Locate_Node_2D

IMPLICIT NONE


INTEGER, DIMENSION(:,:), POINTER			::	t_basic
REAL(8), DIMENSION(:,:), POINTER			::	p_basic
!Yong Modified 2008-12-4 for PC Version
!TYPE(ObjectType), DIMENSION(:), ALLOCATABLE		::	objects
TYPE(ObjectType), DIMENSION(:), INTENT(IN)		::	objects
!INTEGER											nnxp
INTEGER, DIMENSION(:), POINTER				::	element_index

INTEGER, DIMENSION(:,:),POINTER							::  information_1
INTEGER, DIMENSION(:,:),POINTER							::  information_1_temp
INTEGER, DIMENSION(:),POINTER							::  information_1_el
REAL(8), DIMENSION(:,:),POINTER							::  information_2
REAL(8), DIMENSION(:,:),POINTER							::  information_2_temp
REAL(8), DIMENSION(:),POINTER							::  information_2_el
INTEGER, DIMENSION(:), INTENT(IN)						::	node_index

INTEGER		n_elements, n_nodes, ee, ie, n_int_elements,n_objects
INTEGER		el_type, el_region, i, j,n_wall		

REAL(8)		vert(2,4)

REAL(8), DIMENSION(:,:), ALLOCATABLE			::	ObjLowBound, ObjUpBound

INTEGER, DIMENSION(4)						::	node_index_el
LOGICAL											OutOfBound
REAL(8)											P_intrs(4,6)
!=========LY modification, 2022-6-13=========
Real(8) :: Fx1, Fy1, Fx2, Fy2, x_c, y_c
Real(8) :: Dx, Dy, Ex, Ey
!=========LY modification, 2022-6-13=========
!=========LY modification for Test Multi-Object-Fix-Potential, 2022-7-25=========
Integer :: outindex, inindex
!=========LY modification for Test Multi-Object-Fix-Potential, 2022-7-25=========

n_elements	=	SIZE(t_basic,2)
n_nodes		=	SIZE(p_basic,2)
n_objects	=	SIZE(objects,1)

ALLOCATE(ObjLowBound(n_objects,2), ObjUpBound(n_objects,2))

!ALLOCATE(nnode(n_objects))

!ALLOCATE(t_basic_int(5,n_elements))
! t_basic_int(6,:) will be used to store element orientation information (1-5)

ie = 0

n_int_elements	=	0.5 * n_elements ! This number is just an estimation

print *,'initial n_int_elements',n_int_elements

!ALLOCATE(element_index(n_elements))
ALLOCATE(information_1_temp(18,SIZE(t_basic,2)), information_2_temp(8,SIZE(t_basic,2)))

ALLOCATE(information_1_el(18), information_2_el(8))

information_1_temp = 0
information_2_temp = 0

DO i = 1, n_int_elements
	element_index(i)	=	0
END DO

! >>>>>>>>>>>>>>    	 VERY VERY VERY VERY IMPORTANT         <<<<<<<<<<<<<<
! Most exterior objects should come first, then next exterior objects come next.

DO i = 1, n_objects
	IF(objects(i)%Erosion == 0) THEN   !!! ²»¼ÆËã¸¯Ê´ bjw2019-5-30
		IF	(objects(i)%Shape==1 .AND. objects(i)%Axis==0) THEN	!Circle (whole)
			ObjLowBound(i,1)	= objects(i)%Locations(1,1)-objects(i)%Dimensions(1)
			ObjLowBound(i,2)	= objects(i)%Locations(1,2)-objects(i)%Dimensions(1)

			ObjUpBound(i,1)   = objects(i)%Locations(1,1)+objects(i)%Dimensions(1)
			ObjUpBound(i,2)   = objects(i)%Locations(1,2)+objects(i)%Dimensions(1)
		ELSE IF	(objects(i)%Shape==3) THEN	!box
			ObjLowBound(i,1)	= objects(i)%Locations(1,1)
			ObjLowBound(i,2)	= objects(i)%Locations(1,2)

			ObjUpBound(i,1)	  = objects(i)%Locations(2,1)
			ObjUpBound(i,2)	  = objects(i)%Locations(2,2)
        Elseif (objects(i)%Shape==4) Then  !LY modification, 2022-6-13, quadrangle or triangle.
          ObjLowBound(i,1)	= objects(i)%Locations(1,1)
			    ObjLowBound(i,2)	= objects(i)%Locations(1,2)

			    ObjUpBound(i,1)	  = objects(i)%Locations(3,1)
			    ObjUpBound(i,2)	  = objects(i)%Locations(3,2)
        Elseif (objects(i)%Shape==5) Then  !LY modification, 2022-6-13, ellipse.
          !=========LY modification, 2022-6-13=========
          !Note: In below code, we don't use the ObjLowBound and ObjUpBound to judge ellipse condition.
          !      It is similar triangle and arch, we use the node_index_el to judge.
          Fx1 = objects(i)%Locations(1,1)     !Focus1 x of the ellipse
          Fy1 = objects(i)%Locations(1,2)     !Focus1 y of the ellipse
          Fx2 = objects(i)%Locations(2,1)     !Focus2 x of the ellipse
          Fy2 = objects(i)%Locations(2,2)     !Focus2 y of the ellipse
          x_c = (Fx1+Fx2) / 2.0               !Centre x of the ellipse
          y_c = (Fy1+Fy2) / 2.0               !Centre y of the ellipse
      
          ObjLowBound(i,1)	= x_c - objects(i)%Dimensions(1)
			ObjLowBound(i,2)	= y_c - objects(i)%Dimensions(2)

			ObjUpBound(i,1)   = x_c + objects(i)%Dimensions(1)
			ObjUpBound(i,2)   = y_c + objects(i)%Dimensions(2)
      !=========LY modification, 2022-6-13=========
		ELSE IF	(objects(i)%Shape==11) THEN	!LINE
		    ObjLowBound(i,1)	= objects(i)%Locations(1,1)
			ObjLowBound(i,2)	= objects(i)%Locations(1,2)

			ObjUpBound(i,1)	= objects(i)%Locations(2,1)
			ObjUpBound(i,2)	= objects(i)%Locations(2,2)
		ENDIF
  ELSE
    !Not use Erosion, LY modification, 2022-6-13
    Write(*,*) 'Use Erosion, there are a error'
    Stop
		n_wall=objects(i)%Wall(1)%nxp
		IF	(objects(i)%Shape==1 .AND. objects(i)%Axis==0) THEN	!Circle (whole)
			IF (objects(i)%Dimensions(2)==0) THEN
				ObjLowBound(i,1)	= MINVAL(objects(i)%Wall(1)%node(1:n_wall,1))
				ObjLowBound(i,2)	= MINVAL(objects(i)%Wall(1)%node(1:n_wall,2))

				ObjUpBound(i,1)     = MAXVAL(objects(i)%Wall(1)%node(1:n_wall,1))
				ObjUpBound(i,2)	    = MAXVAL(objects(i)%Wall(1)%node(1:n_wall,2))
			ELSE
				ObjLowBound(i,1)	= objects(i)%Locations(1,1)-objects(i)%Dimensions(1)
				ObjLowBound(i,2)	= objects(i)%Locations(1,2)-objects(i)%Dimensions(1)

				ObjUpBound(i,1)     = objects(i)%Locations(1,1)+objects(i)%Dimensions(1)
				ObjUpBound(i,2)	    = objects(i)%Locations(1,2)+objects(i)%Dimensions(1)
			ENDIF
		ELSE IF	(objects(i)%Shape==3) THEN			!box
			IF (objects(i)%Wall(1)%Shape==1 .AND. objects(i)%Wall(1)%Channelwall==4) THEN
				ObjLowBound(i,1)	= objects(i)%Locations(1,1)
				ObjLowBound(i,2)	= objects(i)%Locations(1,2)

				ObjUpBound(i,1)	= objects(i)%Locations(2,1)     !MAXVAL(objects(i)%Wall(1)%node(1:nnode(i)-3,1))
				ObjUpBound(i,2)	= MAXVAL(objects(i)%Wall(1)%node(1:n_wall-3,2))			!objects(i)%Locations(2,2)
			ELSEIF (objects(i)%Wall(1)%Shape==1 .AND. objects(i)%Wall(1)%Channelwall==3) THEN
				ObjLowBound(i,1)	= objects(i)%Locations(1,1)
				ObjLowBound(i,2)	= MINVAL(objects(i)%Wall(1)%node(1:n_wall-3,2))			!MINVAL(ypx(:,j))		!objects(i)%Locations(1,2)

				ObjUpBound(i,1)	= objects(i)%Locations(2,1)         !MAXVAL(objects(i)%Wall(1)%node(1:nnode(i)-3,1))
				ObjUpBound(i,2)	= objects(i)%Locations(2,2)		!objects(i)%Locations(2,2)
			ELSEIF (objects(i)%Wall(1)%Shape==2 ) THEN
				ObjLowBound(i,1)	= MINVAL(objects(i)%Wall(1)%node(1:n_wall,1))
				ObjLowBound(i,2)	= MINVAL(objects(i)%Wall(1)%node(1:n_wall,2))

				ObjUpBound(i,1)	= MAXVAL(objects(i)%Wall(1)%node(1:n_wall,1))
				ObjUpBound(i,2)	= MAXVAL(objects(i)%Wall(1)%node(1:n_wall,2))			!objects(i)%Locations(2,2)
 
			ENDIF
		ENDIF
	ENDIF
!	ENDDO
END DO

!=========Old code, 2022-6-13=========
!DO ee=1,n_elements
!
!	vert = p_basic(:,t_basic(:,ee))
!	node_index_el = node_index(t_basic(:,ee))
!	information_1_el = 0
!	information_2_el = 0
!
!	el_region = 0
!
!	DO i = 1, n_objects
!		OutOfBound = .FALSE.
!		DO j = 1, 2
!			IF (	ALL(vert(j,:)<=ObjLowBound(i,j)) .OR.	&
!					ALL(vert(j,:)>=ObjUpBound(i,j)) ) THEN
!				OutOfBound = .TRUE.
!				EXIT 
!			END IF
!		END DO
!
!		IF (.NOT.OutOfBound) THEN
!
!			information_1_el(1:4) = t_basic(:,ee)
!			information_1_el(5)   = ee
!
!!print *,'pass11',ee,i,objects(i)%Erosion
!
!			IF (objects(i)%Erosion == 0)	THEN	
!				CALL El_Object_Intersection_info_2D_BL(vert, objects(i), P_intrs, el_type, el_region,	&
!													information_1_el, information_2_el, node_index_el)
!			ELSE
!			
!				CALL El_Curve_Intersection_info_2D_BL(vert, objects(i), P_intrs, el_type, el_region,	&
!													information_1_el, information_2_el, node_index_el)
!			    
!			    IF(ee==1 .AND. el_type/=0)THEN
!			        PRINT *,el_region,el_type, node_index_el,i,information_2_el(3:6)
!			    ENDIF
!			ENDIF
!!print *,'pass13',el_region,el_type
!			information_1_temp(:,ee) = information_1_el
!			information_2_temp(:,ee) = information_2_el
!
!		ELSEIF (OutOfBound) THEN
!
!			IF (el_region==0) THEN
!				el_type = 0
!				el_region = min(objects(i)%Regions(2),el_region)
!			ENDIF
!
!		END IF
!	END DO
!	
!   IF (el_type==0) THEN		! Non-interface element
!      element_index(ee)	= el_region
!   ELSE						! Interface element
!      ie = ie+1
!      element_index(ee)	= ie
!
!   END IF
!
!!	t_basic(4, ee) = Element_Index(ee)
!!print*, n_int_elements,ie
!
!END DO
!=========Old code, 2022-6-13=========

!=========LY modification for Multi-Objects-Bool Calculation, 2022-6-13=========
!=========LY modification for Test Multi-Object-Fix-Potential, 2022-7-25=========
outindex = -1
inindex = -2
!=========LY modification for Test Multi-Object-Fix-Potential, 2022-7-25=========
DO ee=1,n_elements

	vert = p_basic(:,t_basic(:,ee))
	node_index_el = node_index(t_basic(:,ee))
	information_1_el = 0
	information_2_el = 0

	el_region = 0

  DO i = 1, n_objects

    IF (node_index_el(1)+node_index_el(2)+node_index_el(3)+node_index_el(4)==4*outindex) THEN !must Non-interface element,outside
      OutOfBound = .TRUE.
      el_region=outindex !-1
      el_type = 0
    !ELSEIF(node_index_el(1)+node_index_el(2)+node_index_el(3)+node_index_el(4)==objects(i)%Regions(1)*4) THEN !must Non-interface element,inside
    ELSEIF(node_index_el(1)<=inindex .AND. node_index_el(2)<=inindex .AND. &
           node_index_el(3)<=inindex .AND. node_index_el(4)<=inindex) THEN   !LY modification for Test Multi-Object-Fix-Potential,2022-7-25
      OutOfBound = .TRUE.
      el_region=inindex !-2
      el_type = 0
    ELSE    !There denote at least one node have a different index.

      IF(i==1)THEN !first time,should judge at least once
        OutOfBound = .FALSE.
      ELSE
        IF (information_1_el(6)==0) THEN !have no intersection so far,can go to judge
          OutOfBound = .FALSE.
        ELSEIF( information_1_el(6)==1 .OR.information_1_el(6)==2) THEN !have intersection ,do not judge from now
          OutOfBound = .TRUE.
        ELSE  !error
          print *, 'error, please check information_1_le(6).'
          print *, 'information_1_el(6) = ', information_1_el(6)
          print *, 'element_index = ', ee
          STOP
        ENDIF
      ENDIF

      If (objects(i)%Shape==1 .OR. objects(i)%Shape==3) Then
        Do j = 1, 2
          If (ALL(vert(j,:)<=ObjLowBound(i,j)) .OR.	&
            ALL(vert(j,:)>=ObjUpBound(i,j)) ) Then
            OutOfBound = .TRUE.
            EXIT
          End If
        End Do
      End If
      
    ENDIF

    IF (.NOT.OutOfBound) THEN

      information_1_el(1:4) = t_basic(:,ee)
      information_1_el(5)   = ee

      IF (objects(i)%Erosion == 0)	THEN
        !Write(*,*) 'Program running here, the element_index = ', ee    !LY Debug using.
        CALL El_Object_Intersection_info_2D_BL(vert, objects(i), P_intrs, el_type, el_region,	&
                                               information_1_el, information_2_el, node_index_el)
      ELSE

        CALL El_Curve_Intersection_info_2D_BL(vert, objects(i), P_intrs, el_type, el_region,	&
                                              information_1_el, information_2_el, node_index_el)

        IF(ee==1 .AND. el_type/=0)THEN
          PRINT *,el_region,el_type, node_index_el,i,information_2_el(3:6)
        ENDIF
      ENDIF

      information_1_temp(:,ee) = information_1_el
      information_2_temp(:,ee) = information_2_el

    ELSEIF (OutOfBound .AND. information_1_el(6)==0) THEN

      IF (el_region==0) THEN
        el_type = 0
        !=========LY modification, 2022-5-6=========
        el_region = min(objects(i)%Regions(2),el_region)
        !=========LY modification, 2022-5-6=========
      ENDIF

    END IF
  END DO
	
  IF (el_type==0) THEN		! Non-interface element
    element_index(ee)	= el_region
  ELSE						! Interface element
    ie = ie+1
    element_index(ee)	= ie
  END IF

END DO
!=========LY modification for Multi-Objects-Bool Calculation, 2022-6-13=========

DO ee=1,n_elements
	IF(element_index(ee) < -1) THEN		!inside the object
		information_1_temp(:,ee) = 0
		information_2_temp(:,ee) = 0
	ENDIF
ENDDO

n_int_elements = ie
ALLOCATE(information_1(18,n_int_elements), information_2(8,n_int_elements))
information_1 = 0
information_2 = 0

i = 0

DO ee=1,n_elements
	IF(element_index(ee) > 0) THEN
		i = i+1
		information_1(:,i) = information_1_temp(:,ee)
		information_2(:,i) = information_2_temp(:,ee)
	ENDIF
ENDDO

!=========LY modification for error presentation in Multi-Objects-Bool calculation, 2022-6-13=========
Do i = 1, Size(information_2,2)
  vert = p_basic(1:2,t_basic(1:4,information_1(5,i)))
  
  Dx = information_2(3,i)
  Dy = information_2(4,i)
  Ex = information_2(5,i)
  Ey = information_2(6,i)
  
  If (Dx-MIN(vert(1,1),vert(1,2))>-SmallValue .AND. Dx-MAX(vert(1,1),vert(1,2))<SmallValue .AND. &
      Dy-MIN(vert(2,1),vert(2,4))>-SmallValue .AND. Dy-MAX(vert(2,1),vert(2,4))<SmallValue .AND. &
      Ex-MIN(vert(1,1),vert(1,2))>-SmallValue .AND. Ex-MAX(vert(1,1),vert(1,2))<SmallValue .AND. &
      Ey-MIN(vert(2,1),vert(2,4))>-SmallValue .AND. Ey-MAX(vert(2,1),vert(2,4))<SmallValue) Then

  Else
    Write(*,*) '=========error:Update_Mesh_Objects_Intersection_2D.f90========='
    Write(*,*) 'error interface element appearance Multi-Objects-Bool calculation!'
    Write(*,*) 'Please checking mesh and objects information and input limitation!'
    Write(*,*) 'error interface element information:'
    Write(*,*) 'element index = ', information_1(5,i)
    Write(*,*) 'the local node vertices = '
    Write(*,*) 'node 1: x = ', vert(1,1), 'y = ', vert(2,1)
    Write(*,*) 'node 2: x = ', vert(1,2), 'y = ', vert(2,2)
    Write(*,*) 'node 3: x = ', vert(1,3), 'y = ', vert(2,3)
    Write(*,*) 'node 4: x = ', vert(1,4), 'y = ', vert(2,4)
    Write(*,*) 'intersection node 1: Dx =', Dx, 'Dy =', Dy
    Write(*,*) 'intersection node 2: Ex =', Ex, 'Ey =', Ey
    Write(*,*) '=========error:Update_Mesh_Objects_Intersection_2D.f90========='
    Stop
  End IF
End Do
!=========LY modification for error presentation in Multi-Objects-Bool calculation, 2022-6-13=========

DEALLOCATE(information_1_temp,information_2_temp)
DEALLOCATE(information_1_el,information_2_el)

DEALLOCATE(ObjLowBound, ObjUpBound)
!DEALLOCATE(nnode)

print *,'end n_int_elements',n_int_elements

END