SUBROUTINE Mesh_Objects_Intersection_info_2D_BL_cell_volume  (	t_basic, p_basic, 	&
												objects, Int_El_Frac,	&
												element_index, information_1, information_2, node_index)

USE Object_Data_2D
USE IFE_INTERFACE, ONLY: El_Object_Intersection_info_2D_BL
IMPLICIT NONE


REAL(8),				 INTENT(IN)						::	Int_El_Frac
INTEGER, DIMENSION(:,:), INTENT(IN)						::	t_basic
REAL(8), DIMENSION(:,:), INTENT(IN)						::	p_basic
TYPE(ObjectType), DIMENSION(:), INTENT(IN)				::	objects
INTEGER, DIMENSION(:), POINTER							::	element_index
INTEGER, DIMENSION(:,:),POINTER							::  information_1
INTEGER, DIMENSION(:,:),POINTER							::  information_1_temp
INTEGER, DIMENSION(:),POINTER							::  information_1_el
REAL(8), DIMENSION(:,:),POINTER							::  information_2
REAL(8), DIMENSION(:,:),POINTER							::  information_2_temp
REAL(8), DIMENSION(:),POINTER							::  information_2_el
INTEGER, DIMENSION(:), INTENT(IN)						::	node_index



INTEGER		n_elements, n_nodes, e, ie, n_int_elements,				&
			el_type, el_region,										&
			i, j, n_objects

REAL(8)		vert(2,4)
REAL(8)		temp1(4), temp2(4)
REAL(8), PARAMETER	::	SmallValue = 1.0D-5

REAL(8), DIMENSION(:,:), POINTER						::	ObjLowBound, ObjUpBound
INTEGER, DIMENSION(4)									::	node_index_el
LOGICAL											OutOfBound
REAL(8)											P_intrs(4,6)


n_elements	=	SIZE(t_basic,2)
n_nodes		=	SIZE(p_basic,2)
n_objects   =   SIZE(objects,1)

ALLOCATE(ObjLowBound(n_objects,2), ObjUpBound(n_objects,2))

ie = 0

n_int_elements	=	Int_El_Frac * n_elements

ALLOCATE(element_index(n_elements))
ALLOCATE(information_1_temp(18,SIZE(t_basic,2)), information_2_temp(8,SIZE(t_basic,2)))

ALLOCATE(information_1_el(18), information_2_el(8))

information_1_temp = 0
information_2_temp = 0


DO i = 1, n_int_elements
	element_index(i)	=	0
END DO


DO i = 1, n_objects
		IF	(objects(i)%Shape==1 .AND. objects(i)%Axis==0) THEN	!Circle (whole)
			ObjLowBound(i,1)	= objects(i)%Locations(1,1)-objects(i)%Dimensions(1)
			ObjLowBound(i,2)	= objects(i)%Locations(1,2)-objects(i)%Dimensions(1)

			ObjUpBound(i,1)     = objects(i)%Locations(1,1)+objects(i)%Dimensions(1)
			ObjUpBound(i,2)	    = objects(i)%Locations(1,2)+objects(i)%Dimensions(1)
		ELSE IF	(objects(i)%Shape==3) THEN	!box
			ObjLowBound(i,1)	= objects(i)%Locations(1,1)
			ObjLowBound(i,2)	= objects(i)%Locations(1,2)

			ObjUpBound(i,1)	= objects(i)%Locations(2,1)
			ObjUpBound(i,2)	= objects(i)%Locations(2,2)
		ELSEIF(objects(i)%Shape==4) THEN   !Triangluar	
			ObjLowBound(i,1)	= MIN(objects(i)%Dimensions(1), objects(i)%Locations(1,1), objects(i)%Locations(2,1))
			ObjLowBound(i,2)	= MIN(objects(i)%Dimensions(2), objects(i)%Locations(1,2), objects(i)%Locations(2,2))

			ObjUpBound(i,1)	= MAX(objects(i)%Dimensions(1), objects(i)%Locations(1,1), objects(i)%Locations(2,1))
			ObjUpBound(i,2)	= MAX(objects(i)%Dimensions(2), objects(i)%Locations(1,2), objects(i)%Locations(2,2))
		ENDIF
	
END DO



DO e=1,n_elements
	vert = p_basic(:,t_basic(:,e))
	node_index_el = node_index(t_basic(:,e))
	information_1_el = 0
	information_2_el = 0
	el_region = 0

	DO i = 1, n_objects
		OutOfBound = .FALSE.
		DO j = 1, 2
			temp1 = vert(j,:) - ObjLowBound(i,j)
			temp2 = vert(j,:) - ObjUpBound(i,j)
			IF (	ALL(temp1 < SmallValue) .OR.	&
					ALL(temp2 > -SmallValue)) THEN
				OutOfBound = .TRUE.
				EXIT 
			END IF
		END DO
		IF (.NOT.OutOfBound) THEN
				information_1_el(1:4) = t_basic(:,e)
				information_1_el(5)   = e
				CALL El_Object_Intersection_info_2D_BL(vert, objects(i), P_intrs, el_type, el_region,	&
														information_1_el, information_2_el, node_index_el)
				information_1_temp(:,e) = information_1_el
				information_2_temp(:,e) = information_2_el
		ELSEIF (OutOfBound) THEN
			IF (el_region==0) THEN
				el_type = 0
				el_region = objects(i)%Regions(2)
			ENDIF
		END IF
	END DO

	
   IF (el_type==0) THEN		! Non-interface element
      element_index(e)	= el_region
   ELSE						    ! Interface element
      ie = ie+1
      element_index(e)	= ie
	  IF(ie>n_int_elements) THEN
		WRITE(6,*) ie,'>', 'Int_El_Frac * n_elements=',n_int_elements
		WRITE(6,*) 'Increase Int_El_Frac',Int_El_Frac 
		STOP
	  ENDIF
   END IF

END DO

DO e=1,n_elements
	IF(element_index(e) < -1) THEN		!inside the object
		information_1_temp(:,e) = 0
		information_2_temp(:,e) = 0
	ENDIF
ENDDO

n_int_elements = ie
ALLOCATE(information_1(18,n_int_elements), information_2(8,n_int_elements))
information_1 = 0
information_2 = 0

i = 0
DO e=1,n_elements
	IF(element_index(e) > 0) THEN
		i = i+1
		information_1(:,i) = information_1_temp(:,e)
		information_2(:,i) = information_2_temp(:,e)
	ENDIF
ENDDO

DEALLOCATE(information_1_temp,information_2_temp)
DEALLOCATE(information_1_el,information_2_el)



DEALLOCATE(ObjLowBound, ObjUpBound)

END